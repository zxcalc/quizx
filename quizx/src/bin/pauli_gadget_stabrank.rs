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
use quizx::circuit::*;
use quizx::graph::*;
use quizx::vec_graph::Graph;
use quizx::decompose::Decomposer;
use rand::rngs::StdRng;
use rand::{SeedableRng, Rng};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let qs = 50;
    let seed = 1337;
    let c = Circuit::random_pauli_gadget()
        .qubits(qs)
        .depth(40)
        .seed(seed)
        .min_weight(2)
        .max_weight(4)
        .build();

    let mut rng = StdRng::seed_from_u64(seed * 37);

    let mut g: Graph = c.to_graph();
    println!("g has T-count: {}", g.tcount());
    quizx::simplify::full_simp(&mut g);
    println!("g has reduced T-count: {}", g.tcount());

    g.plug_inputs(&vec![BasisElem::Z0; qs]);

    for i in 0..qs {
        println!("plugging qubit {}", i);

        let mut h = g.clone();
        h.plug_output(0, BasisElem::Z1);
        h.plug(&h.to_adjoint());

        quizx::simplify::full_simp(&mut h);
        println!("h = g^dag o g has reduced T-count: {}", h.tcount());

        let time = Instant::now();
        let mut d = Decomposer::new(&h);
        d.with_full_simp();
        let max = d.max_terms();
        println!("Naive: {} terms", max);

        println!("Decomposing h...");
        let d = d.decomp_parallel(3);
        // d.decomp_all();
        println!("Finished in {:.2?}", time.elapsed());

        println!("Got {} terms for T-count {} (naive {} terms)", d.nterms, h.tcount(), max);
        println!("{:?}", d.scalar);

        let mut p = d.scalar.float_value().re;
        if p < 0.0 {
            println!("WARNING: p < 0 qubit {}", i);
            p = 0.0;
        } else if p > 1.0 {
            println!("WARNING: p > 0 qubit {}", i);
            p = 1.0;
        }

        if rng.gen_bool(p) {
            g.plug_output(0, BasisElem::Z1);
        } else {
            g.plug_output(0, BasisElem::Z0);
        }
    }

    Ok(())
}
