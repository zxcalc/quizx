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
            qubits: 0,
            clifford_depth: 200,
            n_ccz: 5,
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

    pub fn build(&mut self) -> Circuit {
        let mut oracle = Circuit::new(self.qubits);
        for _ in 0..self.n_ccz {
            self.random_clifford_layer(&mut oracle);
            self.random_ccz(&mut oracle);
        }
        self.random_clifford_layer(&mut oracle);
        for q in 0..self.qubits/2 {
            oracle.push(Gate::new(CZ, vec![q,q+(self.qubits/2)]))
        }

        let mut shift = Circuit::new(self.qubits);
        for q in 0..self.qubits {
            if self.rng.gen_bool(0.5) { shift.push(Gate::new(NOT, vec![q])); }
        }

        let mut hs = Circuit::new(self.qubits);
        for q in 0..self.qubits { hs.push(Gate::new(HAD, vec![q])); }

        let mut c = Circuit::new(self.qubits);
        c += &hs;
        c += &oracle;
        c += &hs;
        c += &shift;
        c += &oracle;
        c += &shift;
        c += &hs;
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
}
