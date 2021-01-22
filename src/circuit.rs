// QuiZX - Rust library for quantum circuit rewriting and optimisation
//         using the ZX-calculus
// Copyright (C) 2021 - Aleks Kissinger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use std::fmt;
use num::Rational;
use num::traits::{Zero,One};
use rustc_hash::FxHashMap;
use regex::Regex;

#[derive(PartialEq,Eq,Clone,Copy,Debug)]
pub enum GType {
    XPhase,
    NOT,
    ZPhase,
    Z,
    S,
    T,
    Sdg,
    Tdg,
    CNOT,
    CZ,
    ParityPhase,
    XCX,
    SWAP,
    HAD,
    Tofolli,
    CCZ,
    InitAncilla,
    PostSelect,
    UnknownGate,
}

#[derive(PartialEq,Eq,Clone,Debug)]
pub struct Gate {
    pub t: GType,
    pub qs: Vec<usize>,
    pub phase: Rational,
}

use GType::*;

impl GType {
    pub fn from_qasm_name(s: &str) -> GType {
        match s {
            "rz"   => ZPhase,
            "rx"   => XPhase,
            "x"    => NOT,
            "z"    => Z,
            "s"    => S,
            "t"    => T,
            "sdg"  => Sdg,
            "tdg"  => Tdg,
            "h"    => HAD,
            "cx"   => CNOT,
            "CX"   => CNOT,
            "cz"   => CZ,
            "ccx"  => Tofolli,
            "ccz"  => CCZ,
            "swap" => SWAP,
            // n.b. these are pyzx-specific gates
            "pp"       => ParityPhase,
            "xcx"      => XCX,
            "init_anc" => InitAncilla,
            "post_sel" => PostSelect,
            _          => UnknownGate,
        }
    }

    pub fn qasm_name(&self) -> &'static str {
        match self {
            ZPhase => "rz",
            NOT => "x",
            XPhase => "rx",
            Z => "z",
            S => "s",
            T => "t",
            Sdg => "sdg",
            Tdg => "tdg",
            HAD => "h",
            CNOT => "cx",
            CZ => "cz",
            Tofolli => "ccx",
            CCZ => "ccz",
            SWAP => "swap",
            // n.b. these are pyzx-specific gates
            ParityPhase => "pp",
            XCX => "xcx",
            InitAncilla => "init_anc",
            PostSelect => "post_sel",
            UnknownGate => "UNKNOWN",
        }
    }
}

/// A type for quantum circuits
#[derive(PartialEq,Eq,Clone,Debug)]
pub struct Circuit {
    nqubits: usize,
    gates: Vec<Gate>
}

impl Gate {
    pub fn from_qasm_name(s: &str) -> Gate {
        Gate {
            t: GType::from_qasm_name(s),
            qs: vec![],
            phase: Rational::zero()
        }
    }

    pub fn qasm_name(&self) -> &'static str { self.t.qasm_name() }

    pub fn to_qasm(&self) -> String {
        let mut s = String::from(self.qasm_name());

        if let ZPhase | XPhase = self.t {
            s += &format!("({}*pi)", self.phase);
        }

        s += " ";
        let qs: Vec<String> = self.qs.iter()
            .map(|i| format!("q[{}]", i)).collect();
        s += &qs.join(", ");

        s
    }

    pub fn adjoint(&mut self) {
        match self.t {
            ZPhase | XPhase | ParityPhase => {
                self.phase *= -1;
            },
            S => { self.t = Sdg },
            T => { self.t = Tdg },
            Sdg => { self.t = S },
            Tdg => { self.t = T },
            _ => {}, // everything else is self-adjoint
        }
    }
}

impl Circuit {
    pub fn new(nqubits: usize) -> Circuit {
        Circuit {
            gates: vec![],
            nqubits
        }
    }

    pub fn num_qubits(&self) -> usize { self.nqubits }

    pub fn push(&mut self, g: Gate) {
        self.gates.push(g);
    }

    pub fn add_gate_with_phase(&mut self, name: &str,
                               qs: Vec<usize>, phase: Rational)
    {
        self.push(Gate {
            t: GType::from_qasm_name(name),
            qs,
            phase});
    }

    pub fn add_gate(&mut self, name: &str, qs: Vec<usize>) {
        self.add_gate_with_phase(name, qs, Rational::zero());
    }

    pub fn adjoint(&mut self) {
        self.gates.reverse();
        for g in &mut self.gates {
            g.adjoint();
        }
    }

    pub fn to_adjoint(&self) -> Circuit {
        let mut c = self.clone();
        c.adjoint();
        c
    }

    pub fn to_qasm(&self) -> String {
        String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n") +
            &self.to_string()
    }

    pub fn parse_qasm(source: &str) -> Result<Vec<qasm::AstNode>, String> {
        // bypass the defalt preprocessor, and just strip comments and includes
        let re = Regex::new(r#"//.*|include\s*"(.*)";"#).unwrap();
        let processed_source = re.replace_all(&source, "");
        let mut tokens = qasm::lex(&processed_source);
        qasm::parse(&mut tokens).map_err(|e| format!("QASM parsing error: {}", e))
    }

    pub fn parse_phase(p: &str) -> Option<Rational> {
        let starts_pi = Regex::new(r#"^\s*(-)\s*pi\s*"#).unwrap();
        let has_pi = Regex::new(r#"\s*\*?\s*pi\s*"#).unwrap();

        if has_pi.is_match(p) {
            let p1 = starts_pi.replace(p, "${1}1");
            let p1 = has_pi.replace(&*p1, "");
            // println!("p1 = '{}'", p1);
            if p1.is_empty() { Some(Rational::one()) }
            else if let Ok(r) = p1.parse::<Rational>() { Some(r) }
            else if let Ok(f) = p1.parse::<f32>() { Rational::approximate_float(f) }
            else { None }
        } else {
            if let Ok(f) = p.parse::<f32>() {
                let f1: f32 = f / std::f32::consts::PI;
                Rational::approximate_float(f1)
            } else { None }
        }
    }

    pub fn from_qasm(source: &str) -> Result<Circuit, String> {
        use qasm::AstNode::{ApplyGate,QReg};
        use qasm::Argument::Qubit;
        let ast = Circuit::parse_qasm(source)?;

        let mut c = Circuit::new(0);
        // maintain a mapping from named registers to qubit offsets
        let mut reg: FxHashMap<String,usize> = FxHashMap::default();

        for nd in ast {
            match nd {
                QReg(rname, sz) => {
                    if reg.contains_key(&rname) {
                        return Err(format!("Re-declaration of qreg: {}", &rname));
                    }

                    reg.insert(rname, c.nqubits);
                    c.nqubits += sz as usize;
                },
                ApplyGate(name, pos, args) => {
                    let t = GType::from_qasm_name(&name);
                    if t == UnknownGate {
                        return Err(format!("Unknown gate: {}", name));
                    }

                    let mut qs: Vec<usize> = Vec::with_capacity(pos.len());
                    for p in pos {
                        if let Qubit(rname, q) = p {
                            if let Some(&offset) = reg.get(&rname) {
                                qs.push(offset + q as usize);
                            } else {
                                return Err(format!("Undeclared register: {}", rname));
                            }
                        } else {
                            return Err(format!("Unsupported: applying gate to entire register"));
                        }
                    }

                    let phase =
                        if args.is_empty() {
                            Rational::zero()
                        } else {
                            if let Some(p) = Circuit::parse_phase(&args[0]) { p }
                            else { return Err(format!("Bad phase: {}", &args[0])); }
                        };
                    c.push(Gate { t, qs, phase });
                },
                // quietly ignore other QASM constructions, rather than giving an "Unsupported" error
                _ => {}, 
            }
        }

        Ok(c)
    }
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "qreg q[{}];\n", self.num_qubits())?;

        for g in &self.gates {
            write!(f, "{};\n", g.to_qasm())?;
        }

        Ok(())
    }
}


impl std::ops::Add<&Circuit> for Circuit {
    type Output = Circuit;
    fn add(mut self, rhs: &Circuit) -> Self::Output {
        if self.num_qubits() != rhs.num_qubits() {
            panic!("Cannot append circuits with different numbers of qubits");
        }
        self.gates.extend_from_slice(&rhs.gates);
        self
    }
}

impl std::ops::Add<Circuit> for Circuit {
    type Output = Circuit;
    fn add(self, rhs: Circuit) -> Self::Output { self.add(&rhs) } }

impl std::ops::Add<Circuit> for &Circuit {
    type Output = Circuit;
    fn add(self, rhs: Circuit) -> Self::Output { self.clone().add(&rhs) } }

