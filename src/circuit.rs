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
use std::str;
use num::{Rational,Zero};
use regex::Regex;
use std::fs::File;
use std::io::prelude::*;
use crate::scalar::Mod2;
use crate::gate::*;

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

    pub fn from_qasm(source: &str) -> Result<Circuit, String> {
        let lines = source.split(';');

        // pattern matching a qreg declaration
        let qreg: Regex = Regex::new(r#"^qreg\s+([a-zA-Z0-9_]+)\s*\[\s*([0-9]+)\s*\]"#).unwrap();

        // the circuit we are building
        let mut c = Circuit::new(0);

        // a mapping from named registers to qubit offset and size. Note a vec
        // with linear lookups seems better than hashing for ~15 or fewer
        // registers.
        let mut reg: Vec<(String,usize,usize)> = Vec::new();

        let mut first = true;

        for line in lines {
            let line = line.trim_start();
            if line.is_empty() { continue; }
            if first && line.starts_with("OPENQASM") { first = false; continue; }

            if line.starts_with("include") { continue; }

            if line.starts_with("qreg") {
                if let Some(caps) = qreg.captures(&line) {
                    let rname = caps[1].to_owned();
                    let sz = caps[2].parse::<usize>().unwrap();
                    // println!("rname = {}, sz = {}", rname, sz);
                    if reg.iter().any(|(s,_,_)| s == &rname) {
                        return Err(format!("Re-declaration of qreg: {}", line));
                    }

                    reg.push((rname, c.nqubits, sz));
                    c.nqubits += sz as usize;
                }
            } else {
                // strip comments and look for a phase
                let mut parts = line.splitn(2, "//").next().unwrap().splitn(2, '(');

                // if a phase is found, this first part of the split is the
                // gate name. If no phase is found, this is the whole command.
                let mut name = parts.next().unwrap().trim_end();

                // continue if this line only contains a comment
                if name.is_empty() { continue; }

                // parsed from argument gate has an argument, otherwise
                // set to 0.
                let phase;

                // save the rest of the command which isn't gate name or arg
                let rest;

                if let Some(arg) = parts.next() {
                    let mut parts = arg.splitn(2, ')');
                    let arg = parts.next().unwrap();
                    if let Some(p) = Circuit::parse_phase(arg) {
                        if let Some(r) = parts.next() {
                            rest = r;
                        } else {
                            return Err(format!("Bad gate application: {}", line));
                        }
                        phase = p.mod2()
                    } else {
                        return Err(format!("Bad phase: {}", line));
                    }
                } else {
                    let mut parts = name.splitn(2, |c| c==' ' || c == '\t');
                    name = parts.next().unwrap();
                    if let Some(r) = parts.next() {
                        rest = r;
                    } else {
                        return Err(format!("Bad gate application: {}", line));
                    }
                    phase = Rational::zero();
                }

                let t = GType::from_qasm_name(&name);
                if t == UnknownGate {
                    return Err(format!("Unknown gate: {}", line));
                }

                let mut qs: Vec<usize> = Vec::new();

                for loc_str in rest.split(",") {
                    let parts = loc_str.trim().splitn(2, "[").collect::<Vec<_>>();

                    if parts.len() == 2 || parts[1].ends_with("]") {
                        let rname = parts[0].to_owned();
                        let q_str = &parts[1][0..parts[1].len()-1];
                        if let Ok(q) = q_str.parse::<usize>() {
                            if let Some(&(_,offset,sz)) = reg.iter().find(|(s,_,_)| s == &rname) {
                                if q < sz {
                                    qs.push(offset + q as usize);
                                } else { return Err(format!("Index out of bounds: {}", line)); }
                            } else { return Err(format!("Undeclared register: {}", line)); }
                        } else { return Err(format!("Expected numeric qubit index: {}", line)); }
                    } else { return Err(format!("Bad qubit location: {}", line)); }
                }

                if let Some(numq) = t.num_qubits() {
                    if numq != qs.len() {
                        return Err(format!("Wrong number of qubits for gate: {}", line));
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
        // let r = BufReader::new(f);
        // let it = r.lines().map(|ln| ln.expect("IO error reading file"));
        // let it = r.split(b';').map(|x|
        //         String::from_utf8(
        //             x.expect("IO error reading file.")
        //         ).expect("Error converting to UTF8."));
        // Circuit::from_lines(source.split(';'))
    }

    /// returns a copy of the circuit, decomposed into 1- and 2-qubit Clifford +
    /// phase gates.
    pub fn to_basic_gates(&self) -> Circuit {
        // calculate the space needed in advance
        let sz = self.gates.iter().map(|g| g.num_basic_gates()).sum();
        let mut gs: Vec<Gate> = Vec::with_capacity(sz);
        for g in &self.gates { g.push_basic_gates(&mut gs); }

        Circuit { gates: gs, nqubits: self.nqubits }
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
