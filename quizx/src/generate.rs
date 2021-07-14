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

use std::mem;

use crate::circuit::*;
use crate::gate::*;
use rand::{SeedableRng, Rng};
use rand::rngs::StdRng;
use num::Rational;

pub struct RandomCircuitBuilder {
    pub rng: StdRng,
    pub qubits: usize,
    pub depth: usize,
    pub p_cnot: f32,
    pub p_cz: f32,
    pub p_h: f32,
    pub p_s: f32,
    pub p_t: f32,
}

pub struct RandomHiddenShiftCircuitBuilder {
    pub rng: StdRng,
    pub qubits: usize,
    pub clifford_depth: usize,
    pub n_ccz: usize,
}

pub struct RandomPauliGadgetCircuitBuilder {
    pub rng: StdRng,
    pub qubits: usize,
    pub depth: usize,
    pub min_weight: usize,
    pub max_weight: usize,
    pub phase_denom: usize,
}

impl Circuit {
    pub fn random() -> RandomCircuitBuilder {
        RandomCircuitBuilder {
            rng:    StdRng::from_entropy(),
            qubits: 0,
            depth:  0,
            p_cnot: 0.0,
            p_cz:   0.0,
            p_h:    0.0,
            p_s:    0.0,
            p_t:    0.0,
        }
    }

    pub fn random_hidden_shift() -> RandomHiddenShiftCircuitBuilder {
        RandomHiddenShiftCircuitBuilder {
            rng: StdRng::from_entropy(),
            qubits: 40,
            clifford_depth: 200,
            n_ccz: 5,
        }
    }

    pub fn random_pauli_gadget() -> RandomPauliGadgetCircuitBuilder {
        RandomPauliGadgetCircuitBuilder {
            rng:         StdRng::from_entropy(),
            qubits:      40,
            depth:       10,
            min_weight:  2,
            max_weight:  4,
            phase_denom: 4,
        }
    }
}

impl RandomCircuitBuilder {
    pub fn seed(&mut self, seed: u64) -> &mut Self { self.rng = StdRng::seed_from_u64(seed); self }
    pub fn qubits(&mut self, qubits: usize) -> &mut Self { self.qubits = qubits; self }
    pub fn depth(&mut self, depth: usize) -> &mut Self { self.depth = depth; self }
    pub fn p_cnot(&mut self, p_cnot: f32) -> &mut Self { self.p_cnot = p_cnot; self }
    pub fn p_cz(&mut self, p_cz: f32) -> &mut Self { self.p_cz = p_cz; self }
    pub fn p_h(&mut self, p_h: f32) -> &mut Self { self.p_h = p_h; self }
    pub fn p_s(&mut self, p_s: f32) -> &mut Self { self.p_s = p_s; self }
    pub fn p_t(&mut self, p_t: f32) -> &mut Self { self.p_t = p_t; self }

    /// Distribute the remaining probability evenly among Clifford (CNOT, H, S) gates
    pub fn with_cliffords(&mut self) -> &mut Self {
        let p = (1.0 - self.p_t - self.p_cz) / 3.0;
        self.p_cnot = p;
        self.p_h = p;
        self.p_s = p;
        self
    }

    /// Convenience method for generating Clifford+T circuits
    ///
    /// Takes a probability of T gates, then distributes the rest
    /// of the probability evenly among CNOT, H, and S.
    pub fn clifford_t(&mut self, p_t: f32) -> &mut Self {
        self.p_t(p_t).with_cliffords()
    }

    pub fn uniform(&mut self) -> &mut Self {
        self.p_cnot = 0.2;
        self.p_cz = 0.2;
        self.p_h = 0.2;
        self.p_s = 0.2;
        self.p_t = 0.2;
        self
    }

    pub fn build(&mut self) -> Circuit {
        let mut c = Circuit::new(self.qubits);

        for _ in 0..self.depth {
            let mut p0 = 0.0;
            let p: f32 = self.rng.gen();
            let q0 = self.rng.gen_range(0..self.qubits);
            let mut q1 = self.rng.gen_range(0..self.qubits-1);
            if q1 >= q0 { q1 += 1; }

            p0 += self.p_cnot;
            if p < p0 { c.push(Gate::new(CNOT, vec![q0,q1])); continue; }

            p0 += self.p_cz;
            if p < p0 { c.push(Gate::new(CZ, vec![q0,q1])); continue; }

            p0 += self.p_h;
            if p < p0 { c.push(Gate::new(HAD, vec![q0])); continue; }

            p0 += self.p_s;
            if p < p0 { c.push(Gate::new(S, vec![q0])); continue; }

            p0 += self.p_t;
            if p < p0 { c.push(Gate::new(T, vec![q0])); continue; }

        }

        c
    }
}

impl RandomHiddenShiftCircuitBuilder {
    pub fn seed(&mut self, seed: u64) -> &mut Self { self.rng = StdRng::seed_from_u64(seed); self }
    pub fn qubits(&mut self, qubits: usize) -> &mut Self { self.qubits = qubits; self }
    pub fn clifford_depth(&mut self, clifford_depth: usize) -> &mut Self { self.clifford_depth = clifford_depth; self }
    pub fn n_ccz(&mut self, n_ccz: usize) -> &mut Self { self.n_ccz = n_ccz; self }

    fn random_clifford_layer(&mut self, c: &mut Circuit) {
        let qs = self.qubits/2;
        for _ in 0..self.clifford_depth {
            let q0 = self.rng.gen_range(0..qs);

            if self.rng.gen_bool(0.5) {
                c.push(Gate::new(Z, vec![q0]));
            } else {
                let mut q1 = self.rng.gen_range(0..qs-1);
                if q1 >= q0 { q1 += 1; }
                c.push(Gate::new(CZ, vec![q0,q1]));
            }
        }
    }

    fn random_ccz(&mut self, c: &mut Circuit) {
        let qs = self.qubits/2;
        let mut q0 = self.rng.gen_range(0..qs);
        let mut q1 = self.rng.gen_range(0..qs-1);
        let mut q2 = self.rng.gen_range(0..qs-2);
        if q1 >= q0 { q1 += 1; } else { mem::swap(&mut q0, &mut q1); }
        if q2 >= q0 { q2 += 1; }
        if q2 >= q1 { q2 += 1; }
        c.push(Gate::new(CCZ, vec![q0,q1,q2]));
    }

    pub fn build(&mut self) -> (Circuit, Vec<u8>) {
        if self.qubits < 6 || self.qubits % 2 != 0 {
            panic!("Random hidden shift circuits must have an even number of qubits >= 6.");
        }
        let mut oraclef = Circuit::new(self.qubits);
        for _ in 0..self.n_ccz {
            self.random_clifford_layer(&mut oraclef);
            self.random_ccz(&mut oraclef);
        }
        self.random_clifford_layer(&mut oraclef);

        let mut oracleg = oraclef.clone();
        oracleg.gates.iter_mut().for_each(|g|
            g.qs.iter_mut().for_each(|q| *q += self.qubits/2)
        );

        for q in 0..self.qubits/2 {
            oraclef.push(Gate::new(CZ, vec![q,q+(self.qubits/2)]));
            oracleg.push(Gate::new(CZ, vec![q,q+(self.qubits/2)]));
        }

        let mut shift = vec![];
        let mut shift_c = Circuit::new(self.qubits);
        for q in 0..self.qubits {
            if self.rng.gen_bool(0.5) {
                shift.push(1);
                shift_c.push(Gate::new(Z, vec![q]));
            } else {
                shift.push(0);
            }
        }

        let mut hs = Circuit::new(self.qubits);
        for q in 0..self.qubits { hs.push(Gate::new(HAD, vec![q])); }

        let mut c = Circuit::new(self.qubits);
        c += &hs;
        c += &oraclef;
        c += &hs;
        c += &shift_c;
        c += &oracleg;
        c += &hs;
        (c, shift)
    }
}

impl RandomPauliGadgetCircuitBuilder {
    pub fn seed(&mut self, seed: u64) -> &mut Self { self.rng = StdRng::seed_from_u64(seed); self }
    pub fn qubits(&mut self, qubits: usize) -> &mut Self { self.qubits = qubits; self }
    pub fn depth(&mut self, depth: usize) -> &mut Self { self.depth = depth; self }
    pub fn min_weight(&mut self, min_weight: usize) -> &mut Self { self.min_weight = min_weight; self }
    pub fn max_weight(&mut self, max_weight: usize) -> &mut Self { self.max_weight = max_weight; self }
    pub fn weight(&mut self, weight: usize) -> &mut Self { self.min_weight = weight; self.max_weight = weight; self }
    pub fn phase_denom(&mut self, phase_denom: usize) -> &mut Self { self.phase_denom = phase_denom; self }

    pub fn build(&mut self) -> Circuit {
        let mut c = Circuit::new(self.qubits);

        // add "depth" pauli gadgets
        for _ in 0..self.depth {
            // pick a random weight in the provided range, inclusive
            let w = self.rng.gen_range(self.min_weight..=self.max_weight);
            if w > self.qubits { panic!("Weight larger than total qubits"); }

            // randomly pick w qubits
            let mut all_qs: Vec<_> = (0..self.qubits).collect();
            let mut qs = vec![];
            for _ in 0..w {
                let q = all_qs.swap_remove(self.rng.gen_range(0..all_qs.len()));
                qs.push(q);
            }
            qs.sort();

            // local clifford circuit randomly applies I, H, or HSH to each chosen qubit to make
            // phase gadget into random pauli gadget
            let mut lc = Circuit::new(self.qubits);
            for &q in &qs {
                match self.rng.gen_range(0..=2) {
                    1 => { lc.push(Gate::new(HAD, vec![q])) },
                    2 => {
                        let mut g = Gate::new(XPhase, vec![q]);
                        g.phase = Rational::new(1,2);
                        lc.push(g);
                    },
                    _ => {},
                }
            }

            // choose random non-clifford (if possible) phase with a fixed denominator
            let phase_num =
                if self.phase_denom >= 4 && self.phase_denom % 2 == 0 {
                    // if the denominator is even, avoid picking rationals equal to 1/2, 1, or 3/2
                    let mut p = self.rng.gen_range(1..(2*self.phase_denom)-3);
                    if p >= self.phase_denom/2 { p += 1; }
                    if p >= self.phase_denom { p += 1; }
                    if p >= (3*self.phase_denom)/2 { p += 1; }

                    p
                } else {
                    // if the denominator is odd or too small, just choose any non-trivial phase
                    self.rng.gen_range(1..2*self.phase_denom)
                };

            let mut g = Gate::new(ParityPhase, qs);
            g.phase = Rational::new(phase_num as isize, self.phase_denom as isize);

            // add pauli gadget to the circuit
            c += &lc;
            c.push(g);
            lc.adjoint();
            c += &lc;
        }

        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_circ() {
        let c = Circuit::random()
            .qubits(5)
            .depth(20)
            .p_cz(0.25)
            .p_cnot(0.25)
            .p_h(0.3)
            .p_s(0.1)
            .p_t(0.1)
            .build();

        assert_eq!(c.num_qubits(), 5);
        assert_eq!(c.num_gates(), 20);
    }

    #[test]
    fn random_seeds() {
        let mut builder = Circuit::random();
        builder.qubits(5)
            .depth(20)
            .uniform();

        builder.seed(1337);
        let c1 = builder.build();

        builder.seed(1337);
        let c2 = builder.build();

        builder.seed(1338);
        let c3 = builder.build();

        assert_eq!(c1, c2);
        assert_ne!(c1, c3);
    }

    #[test]
    fn random_all_gates() {
        // this could fail with some (small) probablity, so try some fixed seeds
        for &seed in &[1337, 800, 40104] {
            let c = Circuit::random().seed(seed)
                .qubits(10)
                .depth(100)
                .with_cliffords()
                .build();
            assert_ne!(c.num_gates_of_type(CNOT), 0);
            assert_eq!(c.num_gates_of_type(CZ), 0);
            assert_ne!(c.num_gates_of_type(HAD), 0);
            assert_ne!(c.num_gates_of_type(S), 0);
            assert_eq!(c.num_gates_of_type(T), 0);

            let c = Circuit::random().seed(seed)
                .qubits(10)
                .depth(100)
                .p_t(0.3)
                .with_cliffords()
                .build();
            assert_ne!(c.num_gates_of_type(CNOT), 0);
            assert_eq!(c.num_gates_of_type(CZ), 0);
            assert_ne!(c.num_gates_of_type(HAD), 0);
            assert_ne!(c.num_gates_of_type(S), 0);
            assert_ne!(c.num_gates_of_type(T), 0);

            let c = Circuit::random().seed(seed)
                .qubits(10)
                .depth(100)
                .uniform()
                .build();
            assert_ne!(c.num_gates_of_type(CNOT), 0);
            assert_ne!(c.num_gates_of_type(CZ), 0);
            assert_ne!(c.num_gates_of_type(HAD), 0);
            assert_ne!(c.num_gates_of_type(S), 0);
            assert_ne!(c.num_gates_of_type(T), 0);
        }
    }

    #[test]
    fn random_hidden_shift() {
        // this could fail with some (small) probablity, so try some fixed seeds
        for &seed in &[1337, 800, 40104] {
            let q = 40;
            let c = Circuit::random_hidden_shift().seed(seed)
                .qubits(q)
                .clifford_depth(200)
                .n_ccz(6)
                .build().0;
            assert_ne!(c.num_gates_of_type(Z), 0);
            assert_ne!(c.num_gates_of_type(CZ), 0);
            assert_eq!(c.num_gates_of_type(HAD), q*3);
            assert_eq!(c.num_gates_of_type(CCZ), 12);
        }
    }

    #[test]
    fn random_pauli_gadget() {
        for &seed in &[1337, 800, 40104] {
            let depth = 100;
            let c = Circuit::random_pauli_gadget().seed(seed)
                .qubits(50)
                .min_weight(4)
                .max_weight(10)
                .depth(depth)
                .phase_denom(8)
                .build();

            assert_eq!(c.num_gates_of_type(ParityPhase), depth);
            assert_ne!(c.num_gates_of_type(HAD), 0);
            assert_ne!(c.num_gates_of_type(XPhase), 0);

            let c = c.to_basic_gates();
            assert!(c.num_gates_of_type(CNOT) >= depth * 4 * 2);
            assert_eq!(c.num_gates_of_type(ZPhase), depth);
        }
    }
}
