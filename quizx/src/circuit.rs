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
use std::collections::VecDeque;
use crate::scalar::Mod2;
use crate::gate::*;
use crate::graph::*;
use crate::linalg::RowOps;

/// A type for quantum circuits
#[derive(PartialEq,Eq,Clone,Debug)]
pub struct Circuit {
    nqubits: usize,
    pub gates: VecDeque<Gate>
}

#[derive(PartialEq,Eq,Clone,Copy,Debug)]
pub struct CircuitStats {
    pub qubits: usize,
    pub total: usize,
    pub oneq: usize,
    pub twoq: usize,
    pub moreq: usize,
    pub cliff: usize,
    pub non_cliff: usize,
}

impl CircuitStats {
    pub fn make(c: &Circuit) -> Self {
        let mut s = CircuitStats { qubits: c.num_qubits(), total: c.num_gates(),
          oneq: 0, twoq: 0, moreq: 0, cliff: 0, non_cliff: 0, };
        for g in &c.gates {
            match g.qs.len() {
                1 => { s.oneq += 1; },
                2 => { s.twoq += 1; },
                _ => { s.moreq += 1; },
            }

            match g.t {
                NOT | Z | S | Sdg | CNOT | CZ | SWAP | HAD => { s.cliff += 1; },
                ZPhase | XPhase => {
                    if *g.phase.denom() == 1 ||
                       *g.phase.denom() == 2
                    { s.cliff += 1; }
                    else { s.non_cliff += 1; }
                },
                _ => { s.non_cliff += 1; }
            }
        }
        s
    }

    pub fn into_array(self) -> [usize; 7] {
        [self.qubits, self.total, self.oneq, self.twoq, self.moreq, self.cliff, self.non_cliff]
    }
}

impl fmt::Display for CircuitStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Circuit with {} qubits, {} gates\n  1-qubit: {}\n  2-qubit: {}\n  n-qubit: {}\n  clifford: {}\n  non-clifford: {}", self.qubits, self.total, self.oneq, self.twoq, self.moreq, self.cliff, self.non_cliff)
    }
}

impl Circuit {
    pub fn new(nqubits: usize) -> Circuit {
        Circuit {
            gates: VecDeque::new(),
            nqubits
        }
    }

    pub fn num_qubits(&self) -> usize { self.nqubits }

    pub fn num_gates(&self) -> usize { self.gates.len() }

    pub fn num_gates_of_type(&self, t: GType) -> usize {
        let mut n = 0;
        for g in &self.gates {
            if g.t == t { n += 1; }
        }
        n
    }

    pub fn push(&mut self, g: Gate) {
        self.gates.push_back(g);
    }

    pub fn push_back(&mut self, g: Gate) {
        self.gates.push_back(g);
    }

    pub fn push_front(&mut self, g: Gate) {
        self.gates.push_front(g);
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

    pub fn reverse(&mut self) {
        self.gates.make_contiguous().reverse();
    }

    pub fn adjoint(&mut self) {
        self.reverse();
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
        let mut c = Circuit { gates: VecDeque::with_capacity(sz), nqubits: self.nqubits };
        for g in &self.gates { g.push_basic_gates(&mut c); }

        c
    }

    pub fn to_graph<G: GraphLike>(&self) -> G {
        let mut graph = G::new();
        let mut qs = Vec::with_capacity(self.nqubits);
        let mut inputs = Vec::with_capacity(self.nqubits);

        // we start counting rows from 1, to allow coordinate
        // (0,0) to mean "no coordinate"
        for i in 0..self.nqubits {
            let v = graph.add_vertex_with_data(VData {
                ty: VType::B,
                phase: Rational::zero(),
                qubit: i as i32,
                row: 1
            });
            qs.push(Some(v));
            inputs.push(v);
        }

        graph.set_inputs(inputs);

        for g in &self.gates {
            g.add_to_graph(&mut graph, &mut qs);
        }

        let last_row = qs.iter()
            .map(|&v| match v { Some(v1) => graph.row(v1), None => 0 })
            .max()
            .unwrap_or(0);

        let mut outputs = Vec::with_capacity(self.nqubits);
        for (i,&q) in qs.iter().enumerate() {
            if let Some(v0) = q {
                let v = graph.add_vertex_with_data(VData {
                    ty: VType::B,
                    phase: Rational::zero(),
                    qubit: i as i32,
                    row: last_row + 1
                });
                graph.add_edge(v0, v);
                outputs.push(v);
            }
        }

        graph.set_outputs(outputs);
        graph
    }

    pub fn stats(&self) -> CircuitStats {
        CircuitStats::make(self)
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


impl std::ops::Add<Circuit> for Circuit {
    type Output = Circuit;
    fn add(mut self, mut rhs: Circuit) -> Self::Output {
        if self.num_qubits() != rhs.num_qubits() {
            panic!("Cannot append circuits with different numbers of qubits");
        }
        self.gates.append(&mut rhs.gates);
        self
    }
}

impl std::ops::Add<&Circuit> for Circuit {
    type Output = Circuit;
    fn add(mut self, rhs: &Circuit) -> Self::Output {
        if self.num_qubits() != rhs.num_qubits() {
            panic!("Cannot append circuits with different numbers of qubits");
        }
        self.gates.extend(rhs.gates.iter().cloned());
        self
    }
}

impl std::ops::Add<Circuit> for &Circuit {
    type Output = Circuit;
    fn add(self, rhs: Circuit) -> Self::Output { self.clone().add(rhs) } }

impl std::ops::Add<&Circuit> for &Circuit {
    type Output = Circuit;
    fn add(self, rhs: &Circuit) -> Self::Output { self.clone().add(rhs) } }

/// A circuit can pretend to be a matrix, where row/column operations correspond
/// to appending or prepending CNOT gates.
///
/// For example, we can synthesise a CNOT circuit corresponding to the parity
/// matrix `m` as follows:
///
/// ```
/// use quizx::circuit::Circuit;
/// use quizx::linalg::*;
/// let mut c = Circuit::new(3); // c|b> = |id * b>
/// let mut m = Mat2::new(vec![vec![1,1,1], vec![0,1,1], vec![0,0,1]]);
/// m.gauss_x(true, 1, &mut c);  // c|b> = |m^-1 * b>
/// c.reverse();                 // c|b> = |m * b>
///
/// ```
impl RowOps for Circuit {
    fn row_add(&mut self, r0: usize, r1: usize) {
        self.push_back(Gate::new(GType::CNOT, vec![r1, r0]));
    }

    fn row_swap(&mut self, r0: usize, r1: usize) {
        self.push_back(Gate::new(GType::SWAP, vec![r0, r1]));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::*;
    use crate::vec_graph::Graph;

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

    #[test]
    fn tograph_cz() {
        let c = Circuit::from_qasm(r#"
        qreg q[2];
        cz q[0], q[1];
        "#).unwrap();

        let h: Graph = c.to_graph();
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_edge(0,2);
        g.add_edge(1,3);
        g.add_edge_with_type(2,3,EType::H);
        g.add_edge(2,4);
        g.add_edge(3,5);
        g.set_inputs(vec![0,1]);
        g.set_outputs(vec![4,5]);
        g.scalar_mut().mul_sqrt2_pow(1);

        println!("g = {}\n\nh = {}\n\n", g.to_dot(), h.to_dot());

        assert_eq!(c.to_tensor4(), Tensor::cphase(Rational::new(1,1), 2));
        assert_eq!(g.to_tensor4(), Tensor::cphase(Rational::new(1,1), 2));
    }

    #[test]
    fn tograph_3cnot() {
        let c = Circuit::from_qasm(r#"
            qreg q[2];
            cx q[0], q[1];
            cx q[1], q[0];
            cx q[0], q[1];
        "#).unwrap();

        let g: Graph = c.to_graph();

        assert_eq!(g.num_vertices(), 10);
        assert_eq!(g.num_edges(), 11);
        assert_eq!(c.to_tensor4(), g.to_tensor4());
    }

    #[test]
    fn tograph_more() {
        let c = Circuit::from_qasm(r#"
            qreg q[3];
            ccz q[0], q[1], q[2];
        "#).unwrap();

        let g: Graph = c.to_graph();
        assert_eq!(c.to_tensor4(), g.to_tensor4());

        let c = Circuit::from_qasm(r#"
            qreg q[4];
            h q[1];
            ccx q[0], q[1], q[2];
            ccx q[1], q[2], q[3];
            h q[2];
            z q[0];
        "#).unwrap();

        let g: Graph = c.to_graph();
        assert_eq!(c.to_tensor4(), g.to_tensor4());
    }
}
