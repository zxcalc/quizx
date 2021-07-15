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
        } else { (50, 50, 4, 4, 1337) };
    if debug { println!("qubits: {}, depth: {}, min_weight: {}, max_weight: {}, seed: {}",
                        qs, depth, min_weight, max_weight, seed); }
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

        // let |h> = (<1| ⊗ I)|g>
        h.plug_output(0, BasisElem::Z1);

        // form <h|h>
        h.plug(&h.to_adjoint());

        quizx::simplify::full_simp(&mut h);
        if debug {
            print!("Q{} ({}T):\t", i, h.tcount());
            io::stdout().flush().unwrap();
        }

        let time = Instant::now();
        let mut d = Decomposer::new(&h);
        d.with_full_simp();

        let d = d.decomp_parallel(3);

        // compute <h|h> by stabiliser decomposition
        prob = d.scalar;
        // if debug { println!("{} / {}", prob, renorm); }

        // n.b. |g> is sub-normalised in general, so let p = <h|h>/<g|g>
        let mut p = prob.float_value().re / renorm.float_value().re;

        // should not happen (unless there are some rounding errors)
        if p < 0.0 {
            println!("WARNING: p < 0 qubit {}. p = {}", i, p);
            p = 0.0;
        } else if p > 1.0 {
            println!("WARNING: p > 1 qubit {}. p = {}", i, p);
            p = 1.0;
        }

        // randomly pick an outcome according to p
        let outcome =
            if rng.gen_bool(p) {
                // outcome 1: let |g> = |h> = (<1| ⊗ I)|g>
                g.plug_output(0, BasisElem::Z1);
                // and save <g|g> = <h|h>
                renorm = prob.clone();
                1
            } else {
                // outcome 0: for |h'> = (<0| ⊗ I)|g>
                //   we have <g|g> = <h'|h'> + <h|h>, so <h'|h'> = <g|g> - <h|h>
                
                // let |g> = |h'>
                g.plug_output(0, BasisElem::Z0);

                // and <g|g> = <h'|h'>
                prob = renorm + Scalar::minus_one() * prob;
                renorm = prob.clone();

                p = 1.0 - p; // complement probability for output below

                0
            };

        meas.push(outcome);

        if debug { println!("{} (p: {}, terms: {}, time: {:.2?})", outcome, p, d.nterms, time.elapsed()); }
    }

    println!("Got: {} (P: {})", meas.iter().format(""), prob);

    Ok(())
}
