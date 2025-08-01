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

use itertools::Itertools;
use quizx::circuit::*;
use quizx::decompose::{terms_for_tcount, Decomposer};
use quizx::graph::*;
use quizx::tensor::*;
use quizx::vec_graph::Graph;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::env;
use std::fs;
use std::io::{self, Write};
use std::time::{Duration, Instant};

fn meas_str(e: &[BasisElem]) -> String {
    format!(
        "{}",
        e.iter()
            .map(|&b| if b == BasisElem::Z0 { 0 } else { 1 })
            .format("")
    )
}

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
        (50, 70, 4, 4, 1, 1337)
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

    let mut rng = StdRng::seed_from_u64(seed * 37);

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
    let mut success = true;
    let mut time = Duration::from_millis(0);

    for s in 1..=nsamples {
        println!("Amplitude {} of {}", s, nsamples);
        g = c.to_graph();
        g.plug_inputs(&vec![BasisElem::Z0; qs]);
        quizx::simplify::full_simp(&mut g);

        // let time_single = Instant::now();
        // let mut terms_single = 0;

        let effect: Vec<_> = (0..qs)
            .map(|_| {
                if rng.random_bool(0.5) {
                    BasisElem::Z0
                } else {
                    BasisElem::Z1
                }
            })
            .collect();

        let mut h = g.clone();
        h.plug_outputs(&effect);
        quizx::simplify::full_simp(&mut h);
        tcounts.push(h.tcount());
        // print!("{} ", h.tcount());

        if debug {
            println!("Plugging {} ({}T)", meas_str(&effect), h.tcount());
            // io::stdout().flush().unwrap();
        }

        let mut d = Decomposer::new(&h);
        d.with_full_simp();

        let d = d.decompose_standard();
        let prob = d.scalar() * d.scalar().conj();
        terms += d.nterms;

        if debug {
            println!(
                "terms = {}, P = {}, re(P) ~ {}",
                d.nterms,
                prob,
                prob.complex_value().re
            );
        }

        time += time_all.elapsed();

        // for small numbers of qubits, it is feasible to check the final probablility
        success = success
            && if qs <= 15 {
                print!("Checking tensor...");
                io::stdout().flush().unwrap();
                let mut check: Graph = c.to_graph();
                check.plug_inputs(&vec![BasisElem::Z0; qs]);
                check.plug_outputs(&effect);
                let amp = check.to_tensor4()[[]];
                let check_prob = amp * amp.conj();
                if check_prob == prob {
                    println!("OK");
                    true
                } else {
                    println!("FAILED {} != {}", check_prob, prob);
                    false
                }
            } else {
                println!("Skipping tensor check (too big).");
                true
            };

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
        "Got {} terms across all samples ({:+e} naive)",
        terms, naive
    );

    let data = format!("\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{:+e}\",\"{:+e}\"\n",
                       qs, depth, tcount, min_weight, max_weight, nsamples, seed, terms, time.as_millis(), tcounts.iter().format(","), no_simp, naive);
    if success {
        print!("OK {}", data);
        fs::write(
            format!(
                "pauli_gadget_{}_{}_{}_{}_{}_{}",
                qs, depth, min_weight, max_weight, nsamples, seed
            ),
            data,
        )
        .expect("Unable to write file");
    } else {
        print!("FAILED {}", data);
    }

    Ok(())
}
