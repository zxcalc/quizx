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

use std::time::Instant;
use std::env;
use std::io::{self,Write};
use itertools::Itertools;
use quizx::circuit::*;
use quizx::graph::*;
use quizx::scalar::*;
use quizx::vec_graph::Graph;
use quizx::decompose::Decomposer;
use rand::rngs::StdRng;
use rand::{SeedableRng, Rng};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let debug = true;
    let args: Vec<_> = env::args().collect();
    let (qs, depth, min_weight, max_weight, seed) =
        if args.len() >= 6 {
            (args[1].parse().unwrap(),
             args[2].parse().unwrap(),
             args[3].parse().unwrap(),
             args[4].parse().unwrap(),
             args[5].parse().unwrap())
        } else { (50, 40, 2, 4, 1337) };
    if debug { println!("qubits: {}, depth: {}, min_weight: {}, max_weight: {}, seed: {}",
                        qs, depth, min_weight, max_weight, seed); }
    let qs = 50;
    let seed = 1337;
    let c = Circuit::random_pauli_gadget()
        .qubits(qs)
        .depth(depth)
        .seed(seed)
        .min_weight(min_weight)
        .max_weight(max_weight)
        .build();

    let mut rng = StdRng::seed_from_u64(seed * 37);

    let mut g: Graph = c.to_graph();
    if debug { println!("g has T-count: {}", g.tcount()); }
    quizx::simplify::full_simp(&mut g);
    if debug { println!("g has reduced T-count: {}", g.tcount()); }

    g.plug_inputs(&vec![BasisElem::Z0; qs]);

    let mut renorm = Scalar::one();
    let mut prob = Scalar::one();
    let mut meas = vec![];

    for i in 0..qs {
        let mut h = g.clone();
        h.plug_output(0, BasisElem::Z1);
        h.plug(&h.to_adjoint());

        quizx::simplify::full_simp(&mut h);
        if debug {
            print!("Q{} ({}T):\t", i, h.tcount());
            io::stdout().flush().unwrap();
        }

        let time = Instant::now();
        let mut d = Decomposer::new(&h);
        d.with_full_simp();
        // let max = d.max_terms();

        let d = d.decomp_parallel(3);
        // d.decomp_all();
        // println!("Finished in {:.2?}", time.elapsed());

        // println!("Got {} terms for T-count {} (naive {} terms)", d.nterms, h.tcount(), max);
        // println!("{:?}", d.scalar);

        prob = d.scalar;
        let mut p = prob.float_value().re / renorm.float_value().re;
        // println!("prob(1): {}", p);
        // println!("prior: {}", renorm);
        // if p < 0.0 {
        //     println!("WARNING: p < 0 qubit {}", i);
        //     p = 0.0;
        // } else if p > 1.0 {
        //     println!("WARNING: p > 0 qubit {}", i);
        //     p = 1.0;
        // }

        let outcome =
            if rng.gen_bool(p) {
                renorm = prob.clone();
                g.plug_output(0, BasisElem::Z1);
                1
            } else {
                p = 1.0 - p;
                prob = renorm + Scalar::minus_one() * prob;
                renorm = prob.clone();
                g.plug_output(0, BasisElem::Z0);
                0
            };

        meas.push(outcome);

        if debug { println!("{} (p: {}, terms: {}, time: {:.2?})", outcome, p, d.nterms, time.elapsed()); }
    }

    println!("Got: {} (P: {})", meas.iter().format(""), prob);

    Ok(())
}
