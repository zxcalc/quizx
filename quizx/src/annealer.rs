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

// use crate::circuit::*;
// use crate::gate::*;
use crate::graph::*;
use crate::extract::*;

use rand::{SeedableRng, Rng};
use rand::rngs::StdRng;

pub struct Annealer<G: GraphLike> {
    pub g: G,
    rng: StdRng,
    scoref: fn (&G) -> usize,
    actions: Vec<fn (&mut StdRng, &mut G)>,
    temp: f64,
    cool: f64,
    iters: usize,
}

impl <G: GraphLike> Annealer<G> {
    pub fn extract_2q_score(g: &G) -> usize {
        let c = g.to_circuit().unwrap();
        c.stats().twoq
    }

    pub fn new(g: G) -> Self {
        Annealer {
            g,
            rng: StdRng::from_entropy(),
            scoref: Annealer::extract_2q_score,
            actions: Vec::new(),
            temp: 25.0,
            cool: 0.005,
            iters: 1000,
        }
    }


    pub fn seed(&mut self, seed: u64) -> &mut Self { self.rng = StdRng::seed_from_u64(seed); self }
    pub fn scoref(&mut self, scoref: fn (&G) -> usize) -> &mut Self
    { self.scoref = scoref; self }

    pub fn temp(&mut self, temp: f64) -> &mut Self { self.temp = temp; self }
    pub fn cool(&mut self, cool: f64) -> &mut Self { self.cool = cool; self }
    pub fn iters(&mut self, iters: usize) -> &mut Self { self.iters = iters; self }

    pub fn anneal(&mut self) {
        if self.actions.is_empty() { return; }
        let mut temp = self.temp;
        let mut current_score = (self.scoref)(&self.g);

        for _ in 0..self.iters {
            // select and action uniformly at random
            let i = self.rng.gen_range(0..self.actions.len());
            let mut g = self.g.clone();
            self.actions[i](&mut self.rng, &mut g);
            let new_score = (self.scoref)(&self.g);
            if new_score < current_score || 
                (temp != 0.0 &&
                 self.rng.gen_bool(f64::min(1.0, ((current_score - new_score) as f64 / temp).exp())))
            {
                self.g = g;
                current_score = new_score;
            }

            temp *= 1.0 - self.cool;
        }
    }
}
