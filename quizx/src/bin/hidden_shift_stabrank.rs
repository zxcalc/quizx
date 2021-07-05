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
// use quizx::tensor::*;
use quizx::scalar::*;
use quizx::vec_graph::Graph;
use quizx::decompose::Decomposer;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let qs = 50;

    let (c,shift) = Circuit::random_hidden_shift()
        .qubits(qs)
        .n_ccz(11)
        .seed(1337)
        .build();
    let mut shift_m = vec![];

    for i in 0..qs {
        let mut g: Graph = c.to_graph_with_options(false);
        g.plug_inputs(&vec![BasisElem::Z0; qs]);
        g.plug_output(i, BasisElem::Z1);
        println!("g has T-count: {}", g.tcount());

        let h = g.to_adjoint();
        g.plug(&h);

        println!("g^dag o g has T-count: {}", g.tcount());
        quizx::simplify::full_simp(&mut g);
        println!("g^dag o g has reduced T-count: {}", g.tcount());

        let time = Instant::now();
        let mut d = Decomposer::new(&g);
        d.with_full_simp();
        let max = d.max_terms();
        println!("Naive: {} terms", max);

        println!("Decomposing g...");
        let d = d.decomp_parallel(3);
        println!("Finished in {:.2?}", time.elapsed());

        println!("Got {} terms for T-count {} (naive {} terms)", d.nterms, g.tcount(), max);
        println!("{:?}", d.scalar);
        if d.scalar.is_zero() { shift_m.push(0); }
        else { shift_m.push(1); }
    }

    println!("       Hidden shift: {:?}", shift);
    println!("Measurement outcome: {:?}", shift_m);

    Ok(())
}
