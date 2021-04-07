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

use quizx::graph::*;
use quizx::vec_graph::Graph;
// use quizx::hash_graph::Graph;
use quizx::basic_rules::*;
use std::time::Instant;

fn main() {
    let sz = 100_000;
    println!("Building Z-spider chain of size: {}...", sz);
    let time = Instant::now();
    let mut g = Graph::new();
    g.add_vertex(VType::Z);

    for i in 1..sz {
        g.add_vertex(VType::Z);
        g.add_edge(i-1, i);
    }

    println!("Done in {:.2?}", time.elapsed());

    assert_eq!(g.num_vertices(), sz);

    println!("Fusing all spiders...");
    let time = Instant::now();

    loop {
        match g.find_edge(|v0,v1,_| check_spider_fusion(&g, v0, v1)) {
            Some((v0,v1,_)) => spider_fusion_unchecked(&mut g, v0, v1),
            None => break,
        };
    }

    println!("Done in {:.2?}", time.elapsed());
    assert_eq!(g.num_vertices(), 1);
}
