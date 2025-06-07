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
use quizx::fscalar::*;
use quizx::graph::*;
use quizx::random_graph::*;
use quizx::tensor::*;
use quizx::vec_graph::Graph;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::env;
use std::fs;
use std::io::{self, Write};
use std::time::{Duration, Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let debug = true;
    let args: Vec<_> = env::args().collect();
    let (qs, depth, min_weight, max_weight, nsamples, norm_est_log_samples, seed) =
        if args.len() >= 6 {
            (
                args[1].parse().unwrap(),
                args[2].parse().unwrap(),
                args[3].parse().unwrap(),
                args[4].parse().unwrap(),
                args[5].parse().unwrap(),
                args[6].parse().unwrap(),
                args[7].parse().unwrap(),
            )
        } else {
            (50, 30, 10, 10, 1, 3, 1337)
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

    let mut random_state_builder = EquatorialStabilizerStateBuilder::new();
    random_state_builder.seed(seed * 51);

    for s in 1..=nsamples {
        println!("Sample {} of {}", s, nsamples);
        g = c.to_graph();
        g.plug_inputs(&vec![BasisElem::Z0; qs]);
        quizx::simplify::full_simp(&mut g);

        let mut renorm = FScalar::one();
        let mut prob = FScalar::one();
        let mut meas = vec![];

        for i in 0..qs {
            if debug {
                print!("Q{}:\t", i);
                io::stdout().flush().unwrap();
            }

            // form <psi|h> for |psi> a random equatorial stabilizer state
            random_state_builder.qubits(g.outputs().len() - 1);

            let time_single = Instant::now();
            let mut mean = FScalar::zero();
            let mut terms_single = 0;

            for _ in 0..2usize.pow(norm_est_log_samples) {
                let mut h = g.clone();
                // let |h> = (<1| ⊗ I)|g>
                h.plug_output(0, BasisElem::Z1);
                let mut psi: Graph = random_state_builder.build();
                psi.adjoint();
                h.plug(&psi);
                quizx::simplify::full_simp(&mut h);
                print!("{} ", h.tcount());
                io::stdout().flush().unwrap();

                tcounts.push(h.tcount());

                let mut d = Decomposer::new(&h);
                d.with_full_simp();

                let d = d.decompose_parallel();
                mean += d.scalar() * d.scalar().conj();
                terms_single += d.nterms;
            }

            let power = 2 * (qs as i32 - norm_est_log_samples as i32);
            mean.mul_sqrt2_pow(power);

            prob = mean;
            // println!("\nprob = {}", prob);
            // if debug { println!("{} / {}", prob, renorm); }

            // n.b. |g> is sub-normalised in general, so let p = <h|h>/<g|g>
            let mut p = prob.complex_value().re / renorm.complex_value().re;

            // should not happen (unless there are some rounding errors)
            if p < 0.0 {
                println!("\nWARNING: p < 0 qubit {}. p = {}", i, p);
                p = 0.0;
            } else if p > 1.0 {
                println!("\nWARNING: p > 1 qubit {}. p = {}", i, p);
                p = 1.0;
            }

            // randomly pick an outcome according to p
            let outcome = if rng.gen_bool(p) {
                // outcome 1: let |g> = |h> = (<1| ⊗ I)|g>
                g.plug_output(0, BasisElem::Z1);
                // and save <g|g> = <h|h>
                renorm = prob;
                1
            } else {
                // outcome 0: for |h'> = (<0| ⊗ I)|g>
                //   we have <g|g> = <h'|h'> + <h|h>, so <h'|h'> = <g|g> - <h|h>

                // let |g> = |h'>
                g.plug_output(0, BasisElem::Z0);

                // and <g|g> = <h'|h'>
                prob = renorm + FScalar::minus_one() * prob;
                renorm = prob;

                p = 1.0 - p; // complement probability for output below

                0
            };

            meas.push(outcome);

            if debug {
                println!(
                    "{} (p: {}, terms: {}, time: {:.2?})",
                    outcome,
                    p,
                    terms_single,
                    time_single.elapsed()
                );
            }
            terms += terms_single;
        }

        println!(
            "Got: {} (P: {}, re(P) ~ {})",
            meas.iter().format(""),
            prob,
            prob.complex_value().re
        );
        time += time_all.elapsed();

        // for small numbers of qubits, it is feasible to check the final probablility
        success = success
            && if qs <= 15 {
                print!("Checking tensor...");
                io::stdout().flush().unwrap();
                let mut check: Graph = c.to_graph();
                let effect: Vec<_> = meas
                    .iter()
                    .map(|&b| if b == 0 { BasisElem::Z0 } else { BasisElem::Z1 })
                    .collect();
                check.plug_inputs(&vec![BasisElem::Z0; qs]);
                check.plug_outputs(&effect);
                let amp = check.to_tensorf()[[]];
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
