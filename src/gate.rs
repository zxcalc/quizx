use num::{Rational,Zero};

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

    fn ccz_decomp(gs: &mut Vec<Gate>, qs: &Vec<usize>) {
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

    /// decompose as 1 and 2 qubit Clifford + phase gates
    ///
    /// If a gate is already basic, return a singleton containing a copy of
    /// self.
    pub fn to_basic_gates(&self) -> Vec<Gate> {
        match self.t {
            CCZ => {
                let mut gs = Vec::with_capacity(13);
                Gate::ccz_decomp(&mut gs, &self.qs);
                gs
            },
            TOFF => {
                let mut gs = Vec::with_capacity(15);
                gs.push(Gate::new(HAD, vec![self.qs[2]]));
                Gate::ccz_decomp(&mut gs, &self.qs);
                gs.push(Gate::new(HAD, vec![self.qs[2]]));
                gs
            },
            ParityPhase => {
                if let Some(&t) = self.qs.last() {
                    let sz = self.qs.len();
                    let mut gs = Vec::with_capacity(sz * 2 - 1);
                    for &c in self.qs[0..sz-1].iter() {
                        gs.push(Gate::new(CNOT, vec![c, t]));
                    }
                    gs.push(Gate::new_with_phase(ZPhase, vec![t], self.phase));
                    for &c in self.qs[0..sz-1].iter().rev() {
                        gs.push(Gate::new(CNOT, vec![c, t]));
                    }
                    gs
                } else {
                    vec![]
                }
            }
            _ => vec![self.clone()],
        }
    }
}
