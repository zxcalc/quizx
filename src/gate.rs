use std::cmp::max;
use num::{Rational,Zero};
use crate::graph::*;

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
    TOFF,
    CCZ,
    InitAncilla,
    PostSelect,
    UnknownGate,
}

pub use GType::*;

#[derive(PartialEq,Eq,Clone,Debug)]
pub struct Gate {
    pub t: GType,
    pub qs: Vec<usize>,
    pub phase: Rational,
}

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
            "ccx"  => TOFF,
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
            TOFF => "ccx",
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

impl Gate {
    pub fn new(t: GType, qs: Vec<usize>) -> Gate {
        Gate { t, qs, phase: Rational::zero() }
    }

    pub fn new_with_phase(t: GType, qs: Vec<usize>, phase: Rational) -> Gate {
        Gate { t, qs, phase }
    }

    fn push_ccz_decomp(gs: &mut Vec<Gate>, qs: &Vec<usize>) {
        gs.push(Gate::new(CNOT, vec![qs[1], qs[2]]));
        gs.push(Gate::new(Tdg, vec![qs[2]]));
        gs.push(Gate::new(CNOT, vec![qs[0], qs[2]]));
        gs.push(Gate::new(T, vec![qs[2]]));
        gs.push(Gate::new(CNOT, vec![qs[1], qs[2]]));
        gs.push(Gate::new(Tdg, vec![qs[2]]));
        gs.push(Gate::new(CNOT, vec![qs[0], qs[2]]));
        gs.push(Gate::new(T, vec![qs[1]]));
        gs.push(Gate::new(T, vec![qs[2]]));
        gs.push(Gate::new(CNOT, vec![qs[0], qs[1]]));
        gs.push(Gate::new(T, vec![qs[0]]));
        gs.push(Gate::new(Tdg, vec![qs[1]]));
        gs.push(Gate::new(CNOT, vec![qs[0], qs[1]]));
    }

    /// number of 1- and 2-qubit Clifford + phase gates needed to realise this gate
    pub fn num_basic_gates(&self) -> usize {
        match self.t {
            CCZ => 13,
            TOFF => 15,
            ParityPhase => if self.qs.is_empty() { 0 } else { self.qs.len() * 2 - 1 },
            _ => 1,
        }
    }

    /// decompose as 1 and 2 qubit Clifford + phase gates and push on to given vec
    ///
    /// If a gate is already basic, push a copy of itself.
    pub fn push_basic_gates(&self, gs: &mut Vec<Gate>) {
        match self.t {
            CCZ => {
                Gate::push_ccz_decomp(gs, &self.qs);
            },
            TOFF => {
                gs.push(Gate::new(HAD, vec![self.qs[2]]));
                Gate::push_ccz_decomp(gs, &self.qs);
                gs.push(Gate::new(HAD, vec![self.qs[2]]));
            },
            ParityPhase => {
                if let Some(&t) = self.qs.last() {
                    let sz = self.qs.len();
                    for &c in self.qs[0..sz-1].iter() {
                        gs.push(Gate::new(CNOT, vec![c, t]));
                    }
                    gs.push(Gate::new_with_phase(ZPhase, vec![t], self.phase));
                    for &c in self.qs[0..sz-1].iter().rev() {
                        gs.push(Gate::new(CNOT, vec![c, t]));
                    }
                }
            }
            _ => gs.push(self.clone()),
        }
    }

    pub fn to_basic_gates(&self) -> Vec<Gate> {
        let mut gates = Vec::with_capacity(self.num_basic_gates());
        self.push_basic_gates(&mut gates);
        gates
    }

    fn add_spider<G: GraphLike>(graph: &mut G, qs: &mut Vec<Option<usize>>, qubit: usize,
                  ty: VType, et: EType, phase: Rational) -> Option<usize>
    {
        if let Some(v0) = qs[qubit] {
            let row = graph.row(v0) + 1;
            let v = graph.add_vertex_with_data(VData { ty, phase, qubit: (qubit as i32), row });
            graph.add_edge_with_type(v0, v, et);
            qs[qubit] = Some(v);
            Some(v)
        } else {
            None
        }
    }

    /// add the gate to the given graph using spiders
    ///
    /// This method takes mutable parameters for the graph being built, and a vec `qs` mapping qubit
    /// number to the most recent vertex in that spot.
    pub fn add_to_graph(&self, graph: &mut impl GraphLike, qs: &mut Vec<Option<usize>>) {
        match self.t {
            ZPhase => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, self.phase); },
            Z      => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(1,1)); },
            S      => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(1,2)); },
            Sdg    => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(-1,2)); },
            T      => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(1,4)); },
            Tdg    => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(-1,4)); },
            XPhase => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, self.phase); },
            NOT    => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(1,1)); },
            HAD    => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::H, Rational::zero()); },
            CNOT => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::X, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    graph.add_edge(v1, v2);
                    graph.scalar().mul_sqrt2_pow(1);
                }
            },
            CZ => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    graph.add_edge_with_type(v1, v2, EType::H);
                    graph.scalar().mul_sqrt2_pow(1);
                }
            },
            XCX => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::X, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    graph.add_edge_with_type(v1, v2, EType::H);
                    graph.scalar().mul_sqrt2_pow(1);
                }
            },
            SWAP => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    qs.swap(self.qs[0], self.qs[1]);
                    Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero());
                    Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Rational::zero());
                }
            },
            InitAncilla => {
                if let Some(v) = qs[self.qs[0]] {
                    // this is a noop if a gate has already been applied to this qubit
                    if graph.vertex_type(v) == VType::B {
                        graph.set_vertex_type(v, VType::X);
                        graph.scalar().mul_sqrt2_pow(-1);
                    }
                }
            },
            PostSelect => {
                Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, Rational::zero());
                graph.scalar().mul_sqrt2_pow(-1);

                // all later gates involving this qubit are quietly ignored
                qs[self.qs[0]] = None;
            },
            CCZ | TOFF | ParityPhase => {
                for g in self.to_basic_gates() {
                    g.add_to_graph(graph, qs);
                }
            }
            UnknownGate => {},
        };
    }
}
