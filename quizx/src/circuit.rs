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

use crate::gate::*;
use crate::graph::*;
use crate::linalg::RowOps;
use crate::phase::Phase;
use num::{Rational64, Zero};
use openqasm::{ast::Symbol, translate::Value, GenericError, ProgramVisitor};
use std::collections::VecDeque;
use std::fmt;
use std::str;

/// A type for quantum circuits
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Circuit {
    nqubits: usize,
    pub gates: VecDeque<Gate>,
}

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
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
        let mut s = CircuitStats {
            qubits: c.num_qubits(),
            total: c.num_gates(),
            oneq: 0,
            twoq: 0,
            moreq: 0,
            cliff: 0,
            non_cliff: 0,
        };
        for g in &c.gates {
            match g.qs.len() {
                1 => {
                    s.oneq += 1;
                }
                2 => {
                    s.twoq += 1;
                }
                _ => {
                    s.moreq += 1;
                }
            }

            match g.t {
                NOT | Z | S | Sdg | CNOT | CZ | SWAP | HAD => {
                    s.cliff += 1;
                }
                ZPhase | XPhase => {
                    if g.phase.is_clifford() {
                        s.cliff += 1;
                    } else {
                        s.non_cliff += 1;
                    }
                }
                _ => {
                    s.non_cliff += 1;
                }
            }
        }
        s
    }

    pub fn into_array(self) -> [usize; 7] {
        [
            self.qubits,
            self.total,
            self.oneq,
            self.twoq,
            self.moreq,
            self.cliff,
            self.non_cliff,
        ]
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
            nqubits,
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.nqubits
    }

    pub fn num_gates(&self) -> usize {
        self.gates.len()
    }

    pub fn num_gates_of_type(&self, t: GType) -> usize {
        let mut n = 0;
        for g in &self.gates {
            if g.t == t {
                n += 1;
            }
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

    pub fn add_gate_with_phase(&mut self, name: &str, qs: Vec<usize>, phase: impl Into<Phase>) {
        self.push(Gate {
            t: GType::from_qasm_name(name),
            qs,
            phase: phase.into(),
        });
    }

    pub fn add_gate(&mut self, name: &str, qs: Vec<usize>) {
        self.add_gate_with_phase(name, qs, Rational64::zero());
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
        String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n") + &self.to_string()
    }

    fn from_qasm_parser(read: impl FnOnce(&mut openqasm::Parser)) -> Result<Circuit, String> {
        let mut cache = openqasm::SourceCache::new();
        let mut parser = openqasm::Parser::new(&mut cache)
            .with_file_policy(openqasm::parser::FilePolicy::Ignore);
        read(&mut parser);
        parser.parse_source::<String>(
            "
            opaque rz(phase) q;
            opaque rx(phase) q;
            opaque x q;
            opaque z q;
            opaque s q;
            opaque t q;
            opaque sdg q;
            opaque tdg q;
            opaque h q;
            opaque cx a, b;
            opaque cz a, b;
            opaque ccx a, b, c;
            opaque ccz a, b, c;
            opaque swap a, b;
            opaque xcx a, b;
            opaque init_anc a;
            opaque post_sel a;
        "
            .to_string(),
            None,
        );

        let program = parser.done().to_errors().map_err(|e| e.to_string())?;
        program
            .type_check()
            .to_errors()
            .map_err(|e| e.to_string())?;

        let mut writer = CircuitWriter {
            circuit: Circuit::new(0),
        };
        let mut linearize = openqasm::Linearize::new(&mut writer, usize::MAX);
        linearize
            .visit_program(&program)
            .to_errors()
            .map_err(|e| e.to_string())?;

        Ok(writer.circuit)
    }

    pub fn from_qasm(source: &str) -> Result<Circuit, String> {
        Circuit::from_qasm_parser(|parser| parser.parse_source::<String>(source.to_string(), None))
    }

    pub fn from_file(name: &str) -> Result<Circuit, String> {
        Circuit::from_qasm_parser(|parser| parser.parse_file(name))
    }

    /// returns a copy of the circuit, decomposed into 1- and 2-qubit Clifford +
    /// phase gates.
    pub fn to_basic_gates(&self) -> Circuit {
        // calculate the space needed in advance
        let sz = self.gates.iter().map(|g| g.num_basic_gates()).sum();
        let mut c = Circuit {
            gates: VecDeque::with_capacity(sz),
            nqubits: self.nqubits,
        };
        for g in &self.gates {
            g.push_basic_gates(&mut c);
        }

        c
    }

    pub fn to_graph_with_options<G: GraphLike>(&self, postselect: bool) -> G {
        let mut graph = G::new();
        let mut qs = Vec::with_capacity(self.nqubits);
        let mut inputs = Vec::with_capacity(self.nqubits);

        // we start counting rows from 1, to allow coordinate
        // (0,0) to mean "no coordinate"
        for i in 0..self.nqubits {
            let v = graph.add_vertex_with_data(VData {
                ty: VType::B,
                phase: Phase::zero(),
                qubit: i as f64,
                row: 1.0,
            });
            qs.push(Some(v));
            inputs.push(v);
        }

        graph.set_inputs(inputs);

        for g in &self.gates {
            g.add_to_graph(&mut graph, &mut qs, postselect);
        }

        let last_row = qs
            .iter()
            .fold(None, |r, &v| match (r, v.map(|v1| graph.row(v1))) {
                (Some(r), Some(r1)) => Some(if r < r1 { r1 } else { r }),
                (Some(r), None) => Some(r),
                (None, Some(r1)) => Some(r1),
                (None, None) => None,
            })
            .unwrap_or(0.0);

        let mut outputs = Vec::with_capacity(self.nqubits);
        for (i, &q) in qs.iter().enumerate() {
            if let Some(v0) = q {
                let v = graph.add_vertex_with_data(VData {
                    ty: VType::B,
                    phase: Phase::zero(),
                    qubit: i as f64,
                    row: last_row + 1.0,
                });
                graph.add_edge(v0, v);
                outputs.push(v);
            }
        }

        graph.set_outputs(outputs);
        graph
    }

    pub fn to_graph<G: GraphLike>(&self) -> G {
        self.to_graph_with_options(false)
    }

    pub fn stats(&self) -> CircuitStats {
        CircuitStats::make(self)
    }
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "qreg q[{}];", self.num_qubits())?;

        for g in &self.gates {
            writeln!(f, "{};", g.to_qasm())?;
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
    fn add(self, rhs: Circuit) -> Self::Output {
        self.clone().add(rhs)
    }
}

impl std::ops::Add<&Circuit> for &Circuit {
    type Output = Circuit;
    fn add(self, rhs: &Circuit) -> Self::Output {
        self.clone().add(rhs)
    }
}

impl std::ops::AddAssign<&Circuit> for Circuit {
    fn add_assign(&mut self, rhs: &Self) {
        self.gates.extend(rhs.gates.iter().cloned());
    }
}

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

struct CircuitWriter {
    circuit: Circuit,
}

#[derive(Debug)]
#[allow(clippy::enum_variant_names)]
enum CircuitWriterError {
    UnitaryNotSupported,
    BarrierNotSupported,
    ResetNotSupported,
    MeasureNotSupported,
    ConditionalNotSupported,
}

impl std::fmt::Display for CircuitWriterError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CircuitWriterError::UnitaryNotSupported => {
                write!(f, "arbitrary unitaries are not supported")
            }
            CircuitWriterError::BarrierNotSupported => write!(f, "barriers are not supported"),
            CircuitWriterError::ResetNotSupported => write!(f, "resets are not supported"),
            CircuitWriterError::MeasureNotSupported => write!(f, "measurements are not supported"),
            CircuitWriterError::ConditionalNotSupported => {
                write!(f, "conditionals are not supported")
            }
        }
    }
}

impl std::error::Error for CircuitWriterError {}

impl openqasm::GateWriter for &mut CircuitWriter {
    type Error = CircuitWriterError;

    fn initialize(&mut self, qubits: &[Symbol], _: &[Symbol]) -> Result<(), Self::Error> {
        self.circuit = Circuit::new(qubits.len());
        Ok(())
    }

    fn write_cx(&mut self, a: usize, b: usize) -> Result<(), Self::Error> {
        self.circuit.push(Gate::new(GType::CNOT, vec![a, b]));
        Ok(())
    }

    fn write_opaque(
        &mut self,
        name: &Symbol,
        params: &[Value],
        regs: &[usize],
    ) -> Result<(), Self::Error> {
        fn param_to_phase(value: Value) -> Phase {
            if value.a.is_zero() {
                Rational64::new(*value.b.numer(), *value.b.denom()).into()
            } else {
                let a = *value.a.numer() as f32 / *value.a.denom() as f32;
                let mut r =
                    Rational64::approximate_float(a / std::f32::consts::PI).unwrap_or(0.into());
                r += Rational64::new(*value.b.numer(), *value.b.denom());
                Phase::new(r)
            }
        }

        let mut g = Gate::from_qasm_name(name.as_str());
        g.qs.extend_from_slice(regs);
        if !params.is_empty() {
            g.phase = param_to_phase(params[0]);
        }

        self.circuit.push(g);

        Ok(())
    }

    fn write_u(&mut self, _: Value, _: Value, _: Value, _: usize) -> Result<(), Self::Error> {
        Err(CircuitWriterError::UnitaryNotSupported)
    }

    fn write_barrier(&mut self, _: &[usize]) -> Result<(), Self::Error> {
        Err(CircuitWriterError::BarrierNotSupported)
    }

    fn write_reset(&mut self, _: usize) -> Result<(), Self::Error> {
        Err(CircuitWriterError::ResetNotSupported)
    }

    fn write_measure(&mut self, _: usize, _: usize) -> Result<(), Self::Error> {
        Err(CircuitWriterError::MeasureNotSupported)
    }

    fn start_conditional(&mut self, _: usize, _: usize, _: u64) -> Result<(), Self::Error> {
        Err(CircuitWriterError::ConditionalNotSupported)
    }

    fn end_conditional(&mut self) -> Result<(), Self::Error> {
        Err(CircuitWriterError::ConditionalNotSupported)
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
        c.add_gate_with_phase("rz", vec![0], Rational64::new(1, 1));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(1, 1));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(1, 3));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(1, 3));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(2, 3));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(2, 3));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(-1, 3));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(-1, 3));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(-1, 3));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(1, 1));
        c.add_gate_with_phase("rz", vec![0], Rational64::new(-1, 2));

        let qasm = r#"
            OPENQASM 2.0;
            include "qelib1.inc";
            qreg q[1];
            rz(pi) q[0];
            rz(-pi) q[0];
            rz(pi/3) q[0];
            rz(1/3 * pi) q[0];
            rz(2*pi/3) q[0];
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
        let c = Circuit::from_qasm(
            r#"
        qreg q[2];
        cz q[0], q[1];
        "#,
        )
        .unwrap();

        let h: Graph = c.to_graph();
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_edge(0, 2);
        g.add_edge(1, 3);
        g.add_edge_with_type(2, 3, EType::H);
        g.add_edge(2, 4);
        g.add_edge(3, 5);
        g.set_inputs(vec![0, 1]);
        g.set_outputs(vec![4, 5]);
        g.scalar_mut().mul_sqrt2_pow(1);

        println!("g = {}\n\nh = {}\n\n", g.to_dot(), h.to_dot());

        assert_eq!(c.to_tensor4(), Tensor::cphase(Rational64::new(1, 1), 2));
        assert_eq!(g.to_tensor4(), Tensor::cphase(Rational64::new(1, 1), 2));
    }

    #[test]
    fn tograph_3cnot() {
        let c = Circuit::from_qasm(
            r#"
            qreg q[2];
            cx q[0], q[1];
            cx q[1], q[0];
            cx q[0], q[1];
        "#,
        )
        .unwrap();

        let g: Graph = c.to_graph();

        assert_eq!(g.num_vertices(), 10);
        assert_eq!(g.num_edges(), 11);
        assert_eq!(c.to_tensor4(), g.to_tensor4());
    }

    #[test]
    fn tograph_more() {
        let c = Circuit::from_qasm(
            r#"
            qreg q[3];
            ccz q[0], q[1], q[2];
        "#,
        )
        .unwrap();

        let g: Graph = c.to_graph();
        assert_eq!(c.to_tensor4(), g.to_tensor4());

        let c = Circuit::from_qasm(
            r#"
            qreg q[4];
            h q[1];
            ccx q[0], q[1], q[2];
            ccx q[1], q[2], q[3];
            h q[2];
            z q[0];
        "#,
        )
        .unwrap();

        let g: Graph = c.to_graph();
        assert_eq!(c.to_tensor4(), g.to_tensor4());
    }

    #[test]
    fn tograph_postsel() {
        let c = Circuit::from_qasm(
            r#"
            qreg q[3];
            ccz q[0], q[1], q[2];
        "#,
        )
        .unwrap();
        let g: Graph = c.to_graph_with_options(true);
        assert_eq!(c.to_tensor4(), g.to_tensor4());

        let c = Circuit::from_qasm(
            r#"
            qreg q[3];
            ccx q[0], q[1], q[2];
        "#,
        )
        .unwrap();
        let g: Graph = c.to_graph_with_options(true);
        assert_eq!(c.to_tensor4(), g.to_tensor4());
    }
}
