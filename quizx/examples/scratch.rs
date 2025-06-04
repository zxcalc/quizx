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

use quizx::circuit::*;
use quizx::graph::*;
use std::time::Instant;
// use quizx::tensor::*;
// use quizx::scalar::*;
use quizx::decompose::Decomposer;
use quizx::vec_graph::Graph;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let qs = 40;
    let c = Circuit::random()
        .qubits(qs)
        .depth(400)
        .seed(1337)
        .p_t(0.4)
        .with_cliffords()
        .build();
    println!("Circuit:\n{}", c.stats());
    let mut g: Graph = c.to_graph();
    g.plug_inputs(&vec![BasisElem::Z0; qs]);
    g.plug_outputs(&vec![BasisElem::Z0; qs]);
    quizx::simplify::full_simp(&mut g);

    println!("Decomposing g with (reduced) T-count: {}", g.tcount());

    let time = Instant::now();
    let mut d = Decomposer::new(&g);
    d.with_full_simp()
        .with_driver(quizx::decompose::Driver::BssTOnly(true));
    let mut max = d.max_terms();
    let mut best_d = d.clone();

    for _ in 0..100 {
        let mut d1 = d.clone();
        d1.decomp_until_depth(1)
            .with_driver(quizx::decompose::Driver::BssTOnly(false))
            .decomp_until_depth(3);
        if d1.max_terms() < max {
            max = d1.max_terms();
            println!("lower max: {}", max);
            best_d = d1;
        }
    }

    d = best_d;
    d.decompose();
    println!("Finished in {:.2?}", time.elapsed());
    println!(
        "got {} terms for T-count {} (naive {} terms)",
        d.nterms,
        g.tcount(),
        max
    );

    // let t = g.to_tensor4();
    // println!("{:?}", t.first().unwrap());
    // println!("{:?}", d.scalar);
    // let s = ScalarN::from_scalar(&t.first().unwrap());

    // assert_eq!(s, d.scalar);

    Ok(())
}
