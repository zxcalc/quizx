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

use crate::circuit::Circuit;
use crate::graph::*;
use crate::params::{Parity, Var};
use crate::phase::Phase;
use crate::scalar::Scalar4;
use num::{Rational64, Zero};

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
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
    TOFF,
    CCZ,
    InitAncilla,
    PostSelect,
    Measure,
    MeasureReset,
    UnknownGate,
}

use rustc_hash::FxHashMap;
pub use GType::*;

impl GType {
    pub fn from_qasm_name(s: &str) -> GType {
        match s {
            "rz" => ZPhase,
            "rx" => XPhase,
            "x" => NOT,
            "z" => Z,
            "s" => S,
            "t" => T,
            "sdg" => Sdg,
            "tdg" => Tdg,
            "h" => HAD,
            "cx" => CNOT,
            "CX" => CNOT,
            "cz" => CZ,
            "ccx" => TOFF,
            "ccz" => CCZ,
            "swap" => SWAP,
            // n.b. these are pyzx-specific gates
            "pp" => ParityPhase,
            "xcx" => XCX,
            "init_anc" => InitAncilla,
            "post_sel" => PostSelect,
            "measure_d" => Measure,
            "measure_r" => MeasureReset,
            _ => UnknownGate,
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
            TOFF => "ccx",
            CCZ => "ccz",
            SWAP => "swap",
            // n.b. these are pyzx-specific gates
            ParityPhase => "pp",
            XCX => "xcx",
            InitAncilla => "init_anc",
            PostSelect => "post_sel",
            Measure => "measure_d",
            MeasureReset => "measure_r",
            UnknownGate => "UNKNOWN",
        }
    }

    /// number of qubits the gate acts on
    ///
    /// If the gate type requires a fixed number of qubits, return it,
    /// otherwise None.
    pub fn num_qubits(&self) -> Option<usize> {
        match self {
            CNOT | CZ | XCX | SWAP => Some(2),
            TOFF | CCZ => Some(3),
            ParityPhase | UnknownGate => None,
            _ => Some(1),
        }
    }
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Gate {
    pub t: GType,
    pub qs: Vec<usize>,
    pub phase: Phase,
    pub vars: Parity,
}

impl Default for Gate {
    fn default() -> Self {
        Gate {
            t: UnknownGate,
            qs: vec![],
            phase: Phase::zero(),
            vars: Parity::zero(),
        }
    }
}

impl Gate {
    pub fn from_qasm_name(s: &str) -> Gate {
        Gate {
            t: GType::from_qasm_name(s),
            qs: vec![],
            phase: Phase::zero(),
            vars: Parity::zero(),
        }
    }

    pub fn qasm_name(&self) -> &'static str {
        self.t.qasm_name()
    }

    pub fn to_qasm(&self) -> String {
        let mut s = String::from(self.qasm_name());

        if let ZPhase | XPhase = self.t {
            s += &format!("({}*pi)", self.phase.to_f64());
        }

        s += " ";
        let qs: Vec<String> = self.qs.iter().map(|i| format!("q[{i}]")).collect();
        s += &qs.join(", ");

        s
    }

    pub fn adjoint(&mut self) {
        match self.t {
            ZPhase | XPhase | ParityPhase => {
                self.phase *= -1;
            }
            S => self.t = Sdg,
            T => self.t = Tdg,
            Sdg => self.t = S,
            Tdg => self.t = T,
            _ => {} // everything else is self-adjoint
        }
    }
}

impl Gate {
    pub fn new(t: GType, qs: Vec<usize>) -> Gate {
        Gate {
            t,
            qs,
            ..Default::default()
        }
    }

    pub fn new_with_phase(t: GType, qs: Vec<usize>, phase: impl Into<Phase>) -> Gate {
        Gate {
            t,
            qs,
            phase: phase.into(),
            ..Default::default()
        }
    }

    pub fn new_with_phase_and_vars(
        t: GType,
        qs: Vec<usize>,
        phase: impl Into<Phase>,
        vars: impl Into<Parity>,
    ) -> Gate {
        Gate {
            t,
            qs,
            phase: phase.into(),
            vars: vars.into(),
        }
    }

    fn push_ccz_decomp(circ: &mut Circuit, qs: &[usize]) {
        circ.push(Gate::new(CNOT, vec![qs[1], qs[2]]));
        circ.push(Gate::new(Tdg, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[2]]));
        circ.push(Gate::new(T, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[1], qs[2]]));
        circ.push(Gate::new(Tdg, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[2]]));
        circ.push(Gate::new(T, vec![qs[1]]));
        circ.push(Gate::new(T, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[1]]));
        circ.push(Gate::new(T, vec![qs[0]]));
        circ.push(Gate::new(Tdg, vec![qs[1]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[1]]));
    }

    /// number of 1- and 2-qubit Clifford + phase gates needed to realise this gate
    pub fn num_basic_gates(&self) -> usize {
        match self.t {
            CCZ => 13,
            TOFF => 15,
            ParityPhase => {
                if self.qs.is_empty() {
                    0
                } else {
                    self.qs.len() * 2 - 1
                }
            }
            _ => 1,
        }
    }

    /// decompose as 1 and 2 qubit Clifford + phase gates and push on to given vec
    ///
    /// If a gate is already basic, push a copy of itself.
    pub fn push_basic_gates(&self, circ: &mut Circuit) {
        match self.t {
            CCZ => {
                Gate::push_ccz_decomp(circ, &self.qs);
            }
            TOFF => {
                circ.push(Gate::new(HAD, vec![self.qs[2]]));
                Gate::push_ccz_decomp(circ, &self.qs);
                circ.push(Gate::new(HAD, vec![self.qs[2]]));
            }
            ParityPhase => {
                if let Some(&t) = self.qs.last() {
                    let sz = self.qs.len();
                    for &c in self.qs[0..sz - 1].iter() {
                        circ.push(Gate::new(CNOT, vec![c, t]));
                    }
                    circ.push(Gate::new_with_phase(ZPhase, vec![t], self.phase));
                    for &c in self.qs[0..sz - 1].iter().rev() {
                        circ.push(Gate::new(CNOT, vec![c, t]));
                    }
                }
            }
            _ => circ.push(self.clone()),
        }
    }

    fn add_spider<G: GraphLike>(
        graph: &mut G,
        qs: &mut FxHashMap<usize, usize>,
        qubit: usize,
        ty: VType,
        et: EType,
        phase: impl Into<Phase>,
    ) -> Option<usize> {
        if let Some(&i) = qs.get(&qubit) {
            let v0 = graph.outputs()[i];
            let qubit = graph.qubit(v0);
            let row = graph.row(v0);
            graph.set_vertex_type(v0, ty);
            graph.set_phase(v0, phase.into());
            let outp = graph.add_vertex_with_data(VData {
                ty: VType::B,
                qubit,
                row: row + 1.0,
                ..Default::default()
            });

            graph.add_edge(v0, outp);
            graph.outputs_mut()[i] = outp;

            if et == EType::H {
                let v1_opt = graph.neighbors(v0).next();
                if let Some(v1) = v1_opt {
                    graph.toggle_edge_type(v0, v1);
                }
            }

            Some(v0)
        } else {
            None
        }
    }

    /// A postselected ZX implementation of CCZ with 4 T-like phases
    ///
    /// Based on the circuit construction of Cody Jones (Phys Rev A 022328, 2013). Note this is intended
    /// only for applications where the circuit doesn't need to be re-extracted (e.g. classical simulation).
    fn add_ccz_postselected<G: GraphLike>(
        graph: &mut G,
        qs: &mut FxHashMap<usize, usize>,
        qubits: &[usize],
    ) {
        if qs.get(&qubits[0]).is_some()
            && qs.get(&qubits[1]).is_some()
            && qs.get(&qubits[2]).is_some()
        {
            let v0 =
                Gate::add_spider(graph, qs, qubits[0], VType::Z, EType::N, Phase::zero()).unwrap();
            let v1 =
                Gate::add_spider(graph, qs, qubits[1], VType::Z, EType::N, Phase::zero()).unwrap();
            let v2 = Gate::add_spider(
                graph,
                qs,
                qubits[2],
                VType::Z,
                EType::N,
                Rational64::new(-1, 2),
            )
            .unwrap();
            // add spiders, 3 in "circuit-like" positions, and one extra
            let s = graph.add_vertex(VType::Z);
            graph.set_phase(s, Rational64::new(-1, 4));
            graph.add_edge_with_type(s, v2, EType::H);

            // add 3 phase gadgets
            let g0 = [
                graph.add_vertex(VType::Z),
                graph.add_vertex(VType::Z),
                graph.add_vertex(VType::Z),
            ];
            let g1 = [
                graph.add_vertex(VType::Z),
                graph.add_vertex(VType::Z),
                graph.add_vertex(VType::Z),
            ];
            graph.set_phase(g1[0], Rational64::new(-1, 4));
            graph.set_phase(g1[1], Rational64::new(-1, 4));
            graph.set_phase(g1[2], Rational64::new(1, 4));
            for i in 0..3 {
                graph.add_edge_with_type(g1[i], g0[i], EType::H);
            }

            // connect gadgets to v0, v1, and s
            graph.add_edge_with_type(g0[0], v0, EType::H);
            graph.add_edge_with_type(g0[0], s, EType::H);
            graph.add_edge_with_type(g0[1], v1, EType::H);
            graph.add_edge_with_type(g0[1], s, EType::H);
            graph.add_edge_with_type(g0[2], v0, EType::H);
            graph.add_edge_with_type(g0[2], v1, EType::H);
            graph.add_edge_with_type(g0[2], s, EType::H);

            // fix scalar
            *graph.scalar_mut() *= Scalar4::new([0, 1, 0, 0], 2);
        }
    }

    /// add the gate to the given graph using spiders
    ///
    /// This method takes mutable parameters for the graph being built, and a vec `qs` mapping qubit
    /// number to the corresponding index of graph.outputs(), which could be different if measurements
    /// or post-selections have happened.
    ///
    /// For basic gates, this returns a set of vertices that have been modified, which can be used to
    /// guide the simplifier. For compound gates, this returns an empty vec.
    ///
    /// TODO: return vertices modified for compound gates.
    pub fn add_to_graph(
        &self,
        fresh_var: &mut Var,
        graph: &mut impl GraphLike,
        qs: &mut FxHashMap<usize, usize>,
        postselect: bool,
    ) -> Vec<V> {
        match self.t {
            ZPhase => Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, self.phase)
                .into_iter()
                .collect(),
            Z => Gate::add_spider(
                graph,
                qs,
                self.qs[0],
                VType::Z,
                EType::N,
                Rational64::new(1, 1),
            )
            .into_iter()
            .collect(),
            S => Gate::add_spider(
                graph,
                qs,
                self.qs[0],
                VType::Z,
                EType::N,
                Rational64::new(1, 2),
            )
            .into_iter()
            .collect(),
            Sdg => Gate::add_spider(
                graph,
                qs,
                self.qs[0],
                VType::Z,
                EType::N,
                Rational64::new(-1, 2),
            )
            .into_iter()
            .collect(),
            T => Gate::add_spider(
                graph,
                qs,
                self.qs[0],
                VType::Z,
                EType::N,
                Rational64::new(1, 4),
            )
            .into_iter()
            .collect(),
            Tdg => Gate::add_spider(
                graph,
                qs,
                self.qs[0],
                VType::Z,
                EType::N,
                Rational64::new(-1, 4),
            )
            .into_iter()
            .collect(),
            XPhase => Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, self.phase)
                .into_iter()
                .collect(),
            NOT => Gate::add_spider(
                graph,
                qs,
                self.qs[0],
                VType::X,
                EType::N,
                Rational64::new(1, 1),
            )
            .into_iter()
            .collect(),
            HAD => Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::H, Phase::zero())
                .into_iter()
                .collect(),
            CNOT => {
                if let (Some(v1), Some(v2)) = (
                    Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Phase::zero()),
                    Gate::add_spider(graph, qs, self.qs[1], VType::X, EType::N, Phase::zero()),
                ) {
                    let o1 = graph.outputs()[*qs.get(&self.qs[0]).unwrap()];
                    let o2 = graph.outputs()[*qs.get(&self.qs[1]).unwrap()];
                    let r1 = graph.row(o1);
                    let r2 = graph.row(o2);
                    let row = if r1 < r2 { r2 } else { r1 };
                    graph.set_row(v1, row - 1.0);
                    graph.set_row(v2, row - 1.0);
                    graph.set_row(o1, row);
                    graph.set_row(o2, row);

                    graph.add_edge(v1, v2);
                    graph.scalar_mut().mul_sqrt2_pow(1);
                    vec![v1, v2]
                } else {
                    vec![]
                }
            }
            CZ => {
                if let (Some(v1), Some(v2)) = (
                    Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Phase::zero()),
                    Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Phase::zero()),
                ) {
                    let o1 = graph.outputs()[*qs.get(&self.qs[0]).unwrap()];
                    let o2 = graph.outputs()[*qs.get(&self.qs[1]).unwrap()];
                    let r1 = graph.row(o1);
                    let r2 = graph.row(o2);
                    let row = if r1 < r2 { r2 } else { r1 };
                    graph.set_row(v1, row - 1.0);
                    graph.set_row(v2, row - 1.0);
                    graph.set_row(o1, row);
                    graph.set_row(o2, row);

                    graph.add_edge_with_type(v1, v2, EType::H);
                    graph.scalar_mut().mul_sqrt2_pow(1);
                    vec![v1, v2]
                } else {
                    vec![]
                }
            }
            XCX => {
                if let (Some(v1), Some(v2)) = (
                    Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, Phase::zero()),
                    Gate::add_spider(graph, qs, self.qs[1], VType::X, EType::N, Phase::zero()),
                ) {
                    let o1 = graph.outputs()[*qs.get(&self.qs[0]).unwrap()];
                    let o2 = graph.outputs()[*qs.get(&self.qs[1]).unwrap()];
                    let r1 = graph.row(o1);
                    let r2 = graph.row(o2);
                    let row = if r1 < r2 { r2 } else { r1 };
                    graph.set_row(v1, row - 1.0);
                    graph.set_row(v2, row - 1.0);
                    graph.set_row(o1, row);
                    graph.set_row(o2, row);

                    graph.add_edge_with_type(v1, v2, EType::H);
                    graph.scalar_mut().mul_sqrt2_pow(1);
                    vec![v1, v2]
                } else {
                    vec![]
                }
            }
            SWAP => {
                if let (Some(&i0), Some(&i1)) = (qs.get(&self.qs[0]), qs.get(&self.qs[1])) {
                    qs.insert(self.qs[0], i1);
                    qs.insert(self.qs[1], i0);
                }
                vec![]
            }
            InitAncilla => {
                if let Some(&i) = qs.get(&self.qs[0]) {
                    let outp = graph.outputs()[i];
                    let inp_opt = graph.neighbors(outp).next();
                    if let Some(inp) = inp_opt {
                        if graph.vertex_type(inp) == VType::B {
                            // this is a noop if a gate has already been applied to this qubit
                            let inputs: Vec<_> = graph
                                .inputs()
                                .iter()
                                .copied()
                                .filter(|&w| w != inp)
                                .collect();
                            graph.set_inputs(inputs);
                            graph.set_vertex_type(inp, VType::X);
                            graph.scalar_mut().mul_sqrt2_pow(-1);
                            return vec![inp];
                        }
                    }
                }
                vec![]
            }
            PostSelect => {
                if let Some(&i) = qs.get(&self.qs[0]) {
                    let outp = graph.outputs()[i];
                    if graph.vertex_type(outp) == VType::B {
                        graph.set_vertex_type(outp, VType::X);

                        // all later gates involving this qubit are quietly ignored
                        graph.outputs_mut().remove(i);
                        qs.remove(&self.qs[0]);

                        // adjust qubit indices to account for the missing output
                        for (_, v1) in qs.iter_mut() {
                            if *v1 > i {
                                *v1 -= 1;
                            }
                        }
                        graph.scalar_mut().mul_sqrt2_pow(-1);
                        return vec![outp];
                    }
                }
                vec![]
            }
            Measure => {
                if let Some(&i) = qs.get(&self.qs[0]) {
                    let outp = graph.outputs()[i];
                    if graph.vertex_type(outp) == VType::B {
                        graph.set_vertex_type(outp, VType::X);

                        if !self.vars.is_empty() {
                            graph.set_vars(outp, self.vars.clone());
                        } else {
                            graph.set_vars(outp, Parity::single(*fresh_var));
                            *fresh_var += 1;
                        }

                        graph.scalar_mut().mul_sqrt2_pow(-1);

                        // all later gates involving this qubit are quietly ignored
                        graph.outputs_mut().remove(i);
                        qs.remove(&self.qs[0]);

                        // adjust qubit indices to account for the missing output
                        for (_, v1) in qs.iter_mut() {
                            if *v1 > i {
                                *v1 -= 1;
                            }
                        }

                        return vec![outp];
                    }
                }
                vec![]
            }
            MeasureReset => {
                if let Some(&i) = qs.get(&self.qs[0]) {
                    let v = graph.outputs()[i];
                    if graph.vertex_type(v) == VType::B {
                        graph.set_vertex_type(v, VType::X);

                        if !self.vars.is_empty() {
                            graph.set_vars(v, self.vars.clone());
                        } else {
                            graph.set_vars(v, Parity::single(*fresh_var));
                            *fresh_var += 1;
                        }

                        let qubit = graph.qubit(v);
                        let row = graph.row(v);

                        let v1 = graph.add_vertex_with_data(VData {
                            ty: VType::X,
                            qubit,
                            row: row + 1.0,
                            ..Default::default()
                        });

                        let outp = graph.add_vertex_with_data(VData {
                            ty: VType::B,
                            qubit,
                            row: row + 2.0,
                            ..Default::default()
                        });

                        graph.add_edge(v1, outp);
                        graph.outputs_mut()[i] = outp;

                        graph.scalar_mut().mul_sqrt2_pow(-2);
                        return vec![v, v1];
                    }
                }
                vec![]
            }
            CCZ => {
                if postselect {
                    Gate::add_ccz_postselected(graph, qs, &self.qs);
                } else {
                    let mut c = Circuit::new(0);
                    self.push_basic_gates(&mut c);
                    for g in c.gates {
                        g.add_to_graph(fresh_var, graph, qs, postselect);
                    }
                }
                vec![]
            }
            TOFF => {
                if postselect {
                    Gate::add_spider(graph, qs, self.qs[2], VType::Z, EType::H, Phase::zero());
                    Gate::add_ccz_postselected(graph, qs, &self.qs);
                    Gate::add_spider(graph, qs, self.qs[2], VType::Z, EType::H, Phase::zero());
                } else {
                    let mut c = Circuit::new(0);
                    self.push_basic_gates(&mut c);
                    for g in c.gates {
                        g.add_to_graph(fresh_var, graph, qs, postselect);
                    }
                }
                vec![]
            }
            ParityPhase => {
                // TODO add directly as phase gadget?
                let mut c = Circuit::new(0);
                self.push_basic_gates(&mut c);
                for g in c.gates {
                    g.add_to_graph(fresh_var, graph, qs, postselect);
                }
                vec![]
            }
            UnknownGate => {
                vec![]
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gtype_from_qasm_name() {
        assert_eq!(GType::from_qasm_name("rz"), ZPhase);
        assert_eq!(GType::from_qasm_name("rx"), XPhase);
        assert_eq!(GType::from_qasm_name("x"), NOT);
        assert_eq!(GType::from_qasm_name("z"), Z);
        assert_eq!(GType::from_qasm_name("s"), S);
        assert_eq!(GType::from_qasm_name("t"), T);
        assert_eq!(GType::from_qasm_name("sdg"), Sdg);
        assert_eq!(GType::from_qasm_name("tdg"), Tdg);
        assert_eq!(GType::from_qasm_name("h"), HAD);
        assert_eq!(GType::from_qasm_name("cx"), CNOT);
        assert_eq!(GType::from_qasm_name("CX"), CNOT);
        assert_eq!(GType::from_qasm_name("cz"), CZ);
        assert_eq!(GType::from_qasm_name("ccx"), TOFF);
        assert_eq!(GType::from_qasm_name("ccz"), CCZ);
        assert_eq!(GType::from_qasm_name("swap"), SWAP);
        assert_eq!(GType::from_qasm_name("unknown"), UnknownGate);
    }

    #[test]
    fn test_gtype_qasm_name_roundtrip() {
        let gate_types = [
            ZPhase, XPhase, NOT, Z, S, T, Sdg, Tdg, HAD, CNOT, CZ, TOFF, CCZ, SWAP, ParityPhase,
            XCX, InitAncilla, PostSelect, Measure, MeasureReset,
        ];

        for gt in gate_types {
            let name = gt.qasm_name();
            // Most gates should roundtrip (except for aliases like CX).
            if gt != CNOT {
                assert_eq!(GType::from_qasm_name(name), gt, "Failed for {:?}", gt);
            }
        }
    }

    #[test]
    fn test_gtype_num_qubits() {
        assert_eq!(ZPhase.num_qubits(), Some(1));
        assert_eq!(XPhase.num_qubits(), Some(1));
        assert_eq!(HAD.num_qubits(), Some(1));
        assert_eq!(T.num_qubits(), Some(1));
        assert_eq!(S.num_qubits(), Some(1));
        assert_eq!(CNOT.num_qubits(), Some(2));
        assert_eq!(CZ.num_qubits(), Some(2));
        assert_eq!(SWAP.num_qubits(), Some(2));
        assert_eq!(XCX.num_qubits(), Some(2));
        assert_eq!(TOFF.num_qubits(), Some(3));
        assert_eq!(CCZ.num_qubits(), Some(3));
        assert_eq!(ParityPhase.num_qubits(), None);
        assert_eq!(UnknownGate.num_qubits(), None);
    }

    #[test]
    fn test_gate_new() {
        let g = Gate::new(CNOT, vec![0, 1]);
        assert_eq!(g.t, CNOT);
        assert_eq!(g.qs, vec![0, 1]);
        assert!(g.phase.is_zero());
    }

    #[test]
    fn test_gate_new_with_phase() {
        let g = Gate::new_with_phase(ZPhase, vec![0], (1, 4));
        assert_eq!(g.t, ZPhase);
        assert_eq!(g.qs, vec![0]);
        assert_eq!(g.phase, Phase::from((1, 4)));
    }

    #[test]
    fn test_gate_default() {
        let g = Gate::default();
        assert_eq!(g.t, UnknownGate);
        assert!(g.qs.is_empty());
        assert!(g.phase.is_zero());
    }

    #[test]
    fn test_gate_to_qasm() {
        let g = Gate::new(CNOT, vec![0, 1]);
        assert_eq!(g.to_qasm(), "cx q[0], q[1]");

        let g = Gate::new(HAD, vec![2]);
        assert_eq!(g.to_qasm(), "h q[2]");

        let g = Gate::new_with_phase(ZPhase, vec![0], (1, 2));
        assert!(g.to_qasm().starts_with("rz("));
        assert!(g.to_qasm().contains("q[0]"));
    }

    #[test]
    fn test_gate_adjoint() {
        // T -> Tdg
        let mut g = Gate::new(T, vec![0]);
        g.adjoint();
        assert_eq!(g.t, Tdg);

        // Tdg -> T
        let mut g = Gate::new(Tdg, vec![0]);
        g.adjoint();
        assert_eq!(g.t, T);

        // S -> Sdg
        let mut g = Gate::new(S, vec![0]);
        g.adjoint();
        assert_eq!(g.t, Sdg);

        // Sdg -> S
        let mut g = Gate::new(Sdg, vec![0]);
        g.adjoint();
        assert_eq!(g.t, S);

        // ZPhase negates phase.
        let mut g = Gate::new_with_phase(ZPhase, vec![0], (1, 4));
        g.adjoint();
        assert_eq!(g.phase, Phase::from((-1, 4)));

        // Self-adjoint gates unchanged
        let mut g = Gate::new(HAD, vec![0]);
        g.adjoint();
        assert_eq!(g.t, HAD);

        let mut g = Gate::new(CNOT, vec![0, 1]);
        g.adjoint();
        assert_eq!(g.t, CNOT);
    }

    #[test]
    fn test_gate_num_basic_gates() {
        assert_eq!(Gate::new(HAD, vec![0]).num_basic_gates(), 1);
        assert_eq!(Gate::new(CNOT, vec![0, 1]).num_basic_gates(), 1);
        assert_eq!(Gate::new(CCZ, vec![0, 1, 2]).num_basic_gates(), 13);
        assert_eq!(Gate::new(TOFF, vec![0, 1, 2]).num_basic_gates(), 15);

        // ParityPhase: 2*n - 1 gates for n qubits
        let g = Gate::new(ParityPhase, vec![0, 1, 2]);
        assert_eq!(g.num_basic_gates(), 5); // 2*3 - 1

        let g = Gate::new(ParityPhase, vec![]);
        assert_eq!(g.num_basic_gates(), 0);
    }

    #[test]
    fn test_gate_push_basic_gates() {
        let mut c = Circuit::new(1);
        let g = Gate::new(HAD, vec![0]);
        g.push_basic_gates(&mut c);
        assert_eq!(c.num_gates(), 1);
        assert_eq!(c.gates[0].t, HAD);

        let mut c = Circuit::new(3);
        let g = Gate::new(CCZ, vec![0, 1, 2]);
        g.push_basic_gates(&mut c);
        assert_eq!(c.num_gates(), 13);

        // TOFF = H on target, then CCZ, then H on target.
        let mut c = Circuit::new(3);
        let g = Gate::new(TOFF, vec![0, 1, 2]);
        g.push_basic_gates(&mut c);
        assert_eq!(c.num_gates(), 15);
        assert_eq!(c.gates[0].t, HAD);
        assert_eq!(c.gates[14].t, HAD);
    }

    #[test]
    fn test_gate_parity_phase_decomposition() {
        // ParityPhase decomposes into CNOT ladder to accumulate parity on last qubit,
        // a ZPhase on that qubit, then the reverse CNOT ladder to restore.
        let mut c = Circuit::new(3);
        let g = Gate::new_with_phase(ParityPhase, vec![0, 1, 2], (1, 4));
        g.push_basic_gates(&mut c);

        assert_eq!(c.num_gates(), 5);
        assert_eq!((c.gates[0].t, &c.gates[0].qs), (CNOT, &vec![0, 2]));
        assert_eq!((c.gates[1].t, &c.gates[1].qs), (CNOT, &vec![1, 2]));
        assert_eq!((c.gates[2].t, &c.gates[2].qs), (ZPhase, &vec![2]));
        assert_eq!((c.gates[3].t, &c.gates[3].qs), (CNOT, &vec![1, 2]));
        assert_eq!((c.gates[4].t, &c.gates[4].qs), (CNOT, &vec![0, 2]));
        assert_eq!(c.gates[2].phase, Phase::from((1, 4)));
    }

    #[test]
    fn test_gate_from_qasm_name() {
        let g = Gate::from_qasm_name("h");
        assert_eq!(g.t, HAD);
        assert!(g.qs.is_empty());

        let g = Gate::from_qasm_name("unknown_gate");
        assert_eq!(g.t, UnknownGate);
    }
}
