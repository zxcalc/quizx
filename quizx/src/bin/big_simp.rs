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
use quizx::circuit::*;
use quizx::vec_graph::*;
use quizx::simplify::*;
use quizx::extract::*;
use std::time::Instant;
// use quizx::tensor::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Circuit::from_file("../circuits/large/hwb10.qasm")?
                    .to_basic_gates();
    println!("stats before: {}", c.stats());
    let mut g: Graph = c.to_graph();

    println!("simplifying...");
    let time = Instant::now();
    clifford_simp(&mut g);
    println!("Done in {:.2?}", time.elapsed());

    let time = Instant::now();
    println!("extracting...");

    let result = g.extractor()
        .gflow()
        .up_to_perm()
        .extract();

    match result {
        Ok(c1) => {
            println!("Done in {:.2?}", time.elapsed());
            println!("extracted ok");
            println!("stats after: {}", c1.stats());
        },
        Err(ExtractError(msg, _c, _g)) => {
            println!("extract failed: {}", msg);
            // println!("{}\n\n{}\n\n{}", msg, _c, _g.to_dot());
        },
    }
    Ok(())
}
