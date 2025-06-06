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
use quizx::decompose::Decomposer;
use quizx::fscalar::*;
use quizx::graph::*;
use quizx::vec_graph::Graph;
use std::env;
use std::fs;
use std::io::{self, Write};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let debug = true;
    let args: Vec<_> = env::args().collect();
    let (qs, n_ccz, seed) = if args.len() >= 4 {
        (
            args[1].parse().unwrap(),
            args[2].parse().unwrap(),
            args[3].parse().unwrap(),
        )
    } else {
        (50, 30, 1337)
    };
    if debug {
        println!("qubits: {}, # ccz: {}, seed: {}", qs, n_ccz, seed);
    }

    // generate hidden shift circuit as in Bravyi-Gosset 2016
    let (c, shift) = Circuit::random_hidden_shift()
        .qubits(qs)
        .n_ccz(n_ccz) // T = CCZ * 2 * 7
        .seed((seed * qs * n_ccz) as u64)
        .build();

    // compute T-count and theoretical max terms
    let mut g: Graph = c.to_graph();
    let tcount = g.tcount();
    g.plug(&g.to_adjoint());
    let mut d = Decomposer::new(&g);
    let naive = d.max_terms();

    let time_all = Instant::now();
    let mut shift_m = vec![];
    let mut terms = 0;
    let mut tcounts = vec![];

    // Hidden shift is deterministic ==> only need 1-qubit marginals

    // compute marginals P(qi = 1)
    for i in 0..qs {
        g = c.to_graph();

        // |00..0> as input
        g.plug_inputs(&vec![BasisElem::Z0; qs]);

        // <1|_qi ⊗ I as output
        g.plug_output(i, BasisElem::Z1);

        // compute norm as <psi|psi>. Doubles T-count!
        g.plug(&g.to_adjoint());

        quizx::simplify::full_simp(&mut g);
        tcounts.push(g.tcount());

        if debug {
            print!("Q{} ({}T):\t", i, g.tcount());
        }
        io::stdout().flush().unwrap();

        let time = Instant::now();

        // do the decomposition, with full_simp called eagerly
        d = Decomposer::new(&g);
        d.with_driver(quizx::decompose::Driver::BssWithCats(false));
        d.with_full_simp();
        let d = d.decompose_parallel();
        terms += d.nterms;

        // record the measurement outcome. Since hidden shift is deterministic, we
        // only need to check if the marginal P(q_i = 1) is zero for each i.
        let outcome = if d.scalar().is_zero() { 0 } else { 1 };
        shift_m.push(outcome);

        if debug {
            println!(
                "{} (terms: {}, time: {:.2?})",
                outcome,
                d.nterms,
                time.elapsed()
            );
        }
    }

    let time = time_all.elapsed();
    if debug {
        println!("Shift: {}", shift.iter().format(""));
        println!("Simul: {}", shift_m.iter().format(""));
        println!("Check: {}", shift == shift_m);
        println!(
            "Circuit with {} qubits and T-count {} simulated in {:.2?}",
            qs, tcount, time
        );
        println!("Got {} terms ({:+e} naive)", terms, (qs as f64) * naive);
    }

    let data = format!(
        "\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\"\n",
        qs,
        n_ccz,
        seed,
        terms,
        time.as_millis(),
        tcounts.iter().format(",")
    );
    if shift == shift_m {
        print!("OK {}", data);
        fs::write(format!("hidden_shift_{}_{}_{}", qs, n_ccz, seed), data)
            .expect("Unable to write file");
    } else {
        print!("FAILED {}", data);
    }
    Ok(())
}
