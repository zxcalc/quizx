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

// use std::time::Instant;
// use quizx::circuit::*;
// use std::fs;
// use quizx::scalar::*;
use quizx::graph::*;
use quizx::tensor::*;
use quizx::vec_graph::Graph;
use quizx::decompose::Decomposer;
use num::Rational;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut g = Graph::new();
    let mut outs = vec![];
    for _ in 0..6 {
        let v = g.add_vertex_with_phase(VType::Z, Rational::new(1,4));
        let w = g.add_vertex(VType::B);
        outs.push(w);
        g.add_edge(v, w);
    }
    g.set_outputs(outs);

    let mut d = Decomposer::new(&g);
    d.decomp_top();

    let t = g.to_tensor4();
    let mut tsum: Tensor4 = d.stack[0].to_tensor4();
    for i in 1..d.stack.len() {
        tsum = tsum + d.stack[i].to_tensor4();
    }

    if t == tsum {
        println!("tensors are equal!");
    }

    Ok(())
}
