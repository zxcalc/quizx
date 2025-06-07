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

use std::env;
use std::time::{Duration, Instant};
// use std::fs;
use num::Complex;
use std::io::{self, Write};
// use itertools::Itertools;
use quizx::circuit::*;
use quizx::graph::*;
// use quizx::scalar::*;
use quizx::decompose::{terms_for_tcount, Decomposer};
use quizx::vec_graph::Graph;
// use rand::rngs::StdRng;
// use rand::{SeedableRng, Rng};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let debug = true;
    let args: Vec<_> = env::args().collect();
    let (qs, depth, min_weight, max_weight, nsamples, seed) = if args.len() >= 6 {
        (
            args[1].parse().unwrap(),
            args[2].parse().unwrap(),
            args[3].parse().unwrap(),
            args[4].parse().unwrap(),
            args[5].parse().unwrap(),
            args[6].parse().unwrap(),
        )
    } else {
        (50, 50, 2, 4, 1, 1339)
        // (13, 15, 2, 4, 3, 1337)
    };
    if debug {
        println!(
            "qubits: {}, depth: {}, min_weight: {}, max_weight: {}, seed: {}",
            qs, depth, min_weight, max_weight, seed
        );
    }
    let c = Circuit::random_pauli_gadget()
        .qubits(qs)
        .depth(depth)
        .seed(seed)
        .min_weight(min_weight)
        .max_weight(max_weight)
        .build();

    // let mut rng = StdRng::seed_from_u64(seed * 37);

    let mut g: Graph = c.to_graph();
    let tcount = g.tcount();
    if debug {
        println!("g has T-count: {}", g.tcount());
    }

    let time_all = Instant::now();
    quizx::simplify::full_simp(&mut g);
    if debug {
        println!("g has reduced T-count: {}", g.tcount());
    }

    let mut tcounts = vec![];
    let mut terms = 0;
    let mut time = Duration::from_millis(0);

    for s in 1..=nsamples {
        println!("Sample {} of {}", s, nsamples);
        g = c.to_graph();
        g.plug_inputs(&vec![BasisElem::Z0; qs]);
        let mut prob;

        for i in 0..qs {
            let mut h = g.clone();

            // let |h> = (<1|_i ⊗ I)|g>
            h.plug_output(i, BasisElem::Z1);

            // form <h|h>
            h.plug(&h.to_adjoint());

            quizx::simplify::full_simp(&mut h);
            if debug {
                print!("Q{} ({}T):\t", i, h.tcount());
                io::stdout().flush().unwrap();
            }

            tcounts.push(h.tcount());

            let time_single = Instant::now();
            let mut d = Decomposer::new(&h);
            d.with_full_simp();

            let d = d.decompose_parallel();

            // compute <h|h> by stabiliser decomposition
            prob = d.scalar();

            if debug {
                let prob_c: Complex<f64> = prob.into();
                println!(
                    "P: {} (re(P) ~ {}, terms: {}, time: {:.2?})",
                    prob,
                    prob_c.re,
                    d.nterms,
                    time_single.elapsed()
                );
            }
            terms += d.nterms;
        }

        // println!("Got: {} (P: {}, re(P) ~ {})", meas.iter().format(""), prob, prob.float_value().re);
        time += time_all.elapsed();

        if debug {
            println!(
                "Circuit with {} qubits and T-count {} simulated in {:.2?}",
                qs, tcount, time
            );
        }
    }

    let naive: f64 = (nsamples as f64) * (qs as f64) * terms_for_tcount(2 * tcount);
    let no_simp: f64 = tcounts.iter().map(|&t| terms_for_tcount(t)).sum();
    println!(
        "Got {} terms across all samples ({:+e} naive, {:+e} no-simp)",
        terms, naive, no_simp
    );

    // let data = format!("\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{:+e}\",\"{:+e}\"\n",
    //                    qs, depth, tcount, min_weight, max_weight, nsamples, seed, terms, time.as_millis(), tcounts.iter().format(","), no_simp, naive);
    // if success {
    //     print!("OK {}", data);
    //     fs::write(&format!("pauli_gadget_{}_{}_{}_{}_{}_{}", qs, depth, min_weight, max_weight, nsamples, seed), data).expect("Unable to write file");
    // } else {
    //     print!("FAILED {}", data);
    // }

    Ok(())
}
