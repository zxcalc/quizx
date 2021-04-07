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


// use quizx::graph::*;
// use quizx::vec_graph::Graph;
// use quizx::tensor::ToTensor;
// use quizx::basic_rules::*;
// use std::time::Instant;
use ndarray::prelude::*;
use ndarray::{Axis,stack};

fn main() {
    // let a1 = array![[1, 0],[0, 1]];
    let a = array![[1, 2],[3, 4]];
    let b = stack![Axis(0),a,a];
    // let sh = Vec::from(b.shape());
    println!("{}", b);
    println!("{:?}", b.shape());
}
