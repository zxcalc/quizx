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
use num::traits::Zero;
use rustc_hash::FxHashMap;
use regex::Regex;
use std::fs::File;
use std::io::Read;
use crate::scalar::Mod2;

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

pub enum QASM {
    QReg(String, usize),
    ApplyGate(String, Vec<(String, usize)>, Vec<String>),
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
    pub gates: Vec<Gate>
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

    pub fn num_gates(&self) -> usize { self.gates.len() }

    pub fn push(&mut self, g: Gate) {
        self.gates.push(g);
    }

    pub fn add_gate_with_phase(&mut self, name: &str,
                               qs: Vec<usize>, phase: Rational)
    {
        self.push(Gate {
            t: GType::from_qasm_name(name),
            qs,
            phase: phase.mod2() });
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

    fn parse_phase(p: &str) -> Option<Rational> {
        let spc = Regex::new(r#"\s*"#).unwrap();
        let starts_pi = Regex::new(r#"^(-?)pi"#).unwrap();
        let has_pi = Regex::new(r#"\*?pi"#).unwrap();

        // strip whitespace
        let p1 = spc.replace_all(p, "");

        if has_pi.is_match(p) {
            // replace (-)pi with (-)1 at the beginning of the string
            let p1 = starts_pi.replace(&p1, "${1}1");
            // remove any other occurance of (*)pi
            let p1 = has_pi.replace(&p1, "");

            // println!("p1 = '{}'", p1);

            if let Ok(r) = p1.parse::<Rational>() { Some(r) }
            else if let Ok(f) = p1.parse::<f32>() { Rational::approximate_float(f) }
            else { None }
        } else {
            if let Ok(f) = p1.parse::<f32>() {
                let f1: f32 = f / std::f32::consts::PI;
                Rational::approximate_float(f1)
            } else { None }
        }
    }

    fn from_qasm(source: &str) -> Result<Circuit, String> {
        // strip comments and includes
        let strip = Regex::new(r#"OPENQASM .*;|//.*|include\s*".*";"#).unwrap();
        let source = strip.replace_all(&source, "");

        // register declaration
        let qreg = Regex::new(r#"^qreg\s+([a-zA-Z0-9_]+)\s*\[\s*([0-9]+)\s*\]$"#).unwrap();
        // gate application
        let gate = Regex::new(r#"^([a-zA-Z0-9_]+)\s*(\(([^)]*)\))?\s+(.*)$"#).unwrap();
        // a qubut location "REG[x]"
        let loc = Regex::new(r#"^\s*([a-zA-Z0-9_]+)\s*\s*\[\s*([0-9]+)\s*\]\s*$"#).unwrap();

        // the circuit we are building
        let mut c = Circuit::new(0);
        // a mapping from named registers to qubit offsets
        let mut reg: FxHashMap<String,usize> = FxHashMap::default();

        for line in source.split(";") {
            let line = line.trim();
            if line.is_empty() { continue; }
            if let Some(caps) = qreg.captures(line) {
                let rname = caps[1].to_owned();
                let sz = caps[2].parse::<usize>().unwrap();
                // println!("rname = {}, sz = {}", rname, sz);
                if reg.contains_key(&rname) {
                    return Err(format!("Re-declaration of qreg: {}", &rname));
                }

                reg.insert(rname, c.nqubits);
                c.nqubits += sz as usize;
            } else if let Some(caps) = gate.captures(line) {
                let name = caps[1].to_owned();

                let phase =
                    if let Some(m) = caps.get(3) {
                        if let Some(p) = Circuit::parse_phase(m.as_str()) { p.mod2() }
                        else { return Err(format!("Bad phase: {}", m.as_str())); }
                    } else {
                        Rational::zero()
                    };

                // println!("name = {}, phase = {}, rest = {}", name, phase, caps[4].to_owned());

                let t = GType::from_qasm_name(&name);
                if t == UnknownGate {
                    return Err(format!("Unknown gate: {}", name));
                }

                let mut qs: Vec<usize> = Vec::new();
                for loc_str in caps[4].split(",") {
                    if let Some(cap1) = loc.captures(loc_str) {
                        let rname = cap1[1].to_owned();
                        let q = cap1[2].parse::<usize>().unwrap();
                        if let Some(&offset) = reg.get(&rname) {
                            qs.push(offset + q as usize);
                        } else {
                            return Err(format!("Undeclared register: {}", rname));
                        }
                    } else {
                        return Err(format!("Bad qubit location: {}", loc_str));
                    }
                }

                c.push(Gate { t, qs, phase });
            }
        }
        // let mut tokens = qasm::lex(&processed_source);
        // qasm::parse(&mut tokens).map_err(|e| format!("QASM parsing error: {}", e))
        Ok(c)
    }

    pub fn from_file(name: &str) -> Result<Circuit, String> {
        let mut f = File::open(name).map_err(|e| e.to_string())?;
        let mut source = String::new();
        f.read_to_string(&mut source).map_err(|e| e.to_string())?;
        Circuit::from_qasm(&source)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mk_circuit() {
        let mut c = Circuit::new(3);
        c.add_gate("cz", vec![0, 1]);
        c.add_gate("z", vec![1]);
        c.add_gate("cx", vec![1, 2]);
        c.add_gate("h", vec![0]);
        assert_eq!(c.num_qubits(), 3);
        assert_eq!(c.num_gates(), 4);

        let qasm = r#"
            OPENQASM 2.0;
            include "qelib1.inc";
            qreg q[3];
            cz q[0], q[1];
            z q[1];
            cx q[1], q[2];
            h q[0];
        "#;

        let c1 = Circuit::from_qasm(qasm);
        assert_eq!(c1, Ok(c));
    }

    #[test]
    fn mk_circuit_with_phase() {
        let mut c = Circuit::new(1);
        c.add_gate_with_phase("rz", vec![0], Rational::new(1,1));
        c.add_gate_with_phase("rz", vec![0], Rational::new(1,1));
        c.add_gate_with_phase("rz", vec![0], Rational::new(1,3));
        c.add_gate_with_phase("rz", vec![0], Rational::new(1,3));
        c.add_gate_with_phase("rz", vec![0], Rational::new(2,3));
        c.add_gate_with_phase("rz", vec![0], Rational::new(2,3));
        c.add_gate_with_phase("rz", vec![0], Rational::new(-1,3));
        c.add_gate_with_phase("rz", vec![0], Rational::new(-1,3));
        c.add_gate_with_phase("rz", vec![0], Rational::new(-1,3));
        c.add_gate_with_phase("rz", vec![0], Rational::new(1,1));
        c.add_gate_with_phase("rz", vec![0], Rational::new(-1,2));

        let qasm = r#"
            OPENQASM 2.0;
            include "qelib1.inc";
            qreg q[1];
            rz(pi) q[0];
            rz(-pi) q[0];
            rz(pi/3) q[0];
            rz(1/3 * pi) q[0];
            rz(2pi/3) q[0];
            rz(2/3 * pi) q[0];
            rz(-pi/3) q[0];
            rz(-1/3 * pi) q[0];
            rz(-0.333333333 * pi) q[0];
            rz(3.14159265359) q[0];
            rz(-1.57079632679) q[0];
        "#;

        let c1 = Circuit::from_qasm(qasm);
        assert_eq!(c1, Ok(c));
    }

    #[test]
    fn mk_circuit_2reg() {
        let mut c = Circuit::new(5);
        c.add_gate("cx", vec![0, 1]);
        c.add_gate("cx", vec![1, 2]);
        c.add_gate("cx", vec![2, 3]);
        c.add_gate("cx", vec![3, 4]);

        let qasm = r#"
            OPENQASM 2.0;
            include "qelib1.inc";
            qreg q[2];
            qreg r[3];
            cx q[0], q[1];
            cx q[1], r[0];
            cx r[0], r[1];
            cx r[1], r[2];
        "#;

        let c1 = Circuit::from_qasm(qasm);
        assert_eq!(c1, Ok(c));
    }
}
