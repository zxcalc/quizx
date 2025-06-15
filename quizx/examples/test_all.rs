// QuiZX - Rust library for quantum circuit rewriting and optimisation
//         using the ZX-calculus
// Copyright (C) 2021 - Vlad-mihai Moldoveanu, Aleks Kissinger
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
use quizx::extract::*;
use quizx::simplify::*;
use quizx::tensor::*;
use quizx::vec_graph::*;
use std::fs;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for e in fs::read_dir("circuits/small")? {
        if let Some(f) = e?.path().to_str() {
            let time = Instant::now();
            println!("{}", f);
            let c = Circuit::from_file(f)
                .unwrap_or_else(|e| panic!("circuit failed to parse: {}. {}", f, e));
            println!("...done reading in {:.2?}", time.elapsed());
            // if c.num_qubits() > 10 { continue; }

            println!("Simplifying circuit...");
            let time = Instant::now();
            let mut g: Graph = c.to_graph();
            full_simp(&mut g);
            println!("Done simplifying in {:.2?}", time.elapsed());

            println!("Extracting circuit...");
            let time = Instant::now();
            let result = g
                .extractor()
                .gflow()
                // .up_to_perm()
                .extract();

            match result {
                Ok(c1) => {
                    println!("Done in {:.2?}", time.elapsed());
                    println!("extracted ok");

                    if c1.num_qubits() > 5 {
                        println!("Circuit too big, not comparing tensors.");
                        continue;
                    }

                    println!("Comparing tensors...");
                    let time = Instant::now();
                    if !Tensor4::scalar_compare(&c, &c1) {
                        println!("Tensors differ!");
                    } else {
                        println!("Checked successfully in {:.2?}", time.elapsed());
                    }
                }
                Err(ExtractError(msg, _c, _g)) => {
                    println!("extract failed: {}", msg);
                    // println!("{}\n\n{}\n\n{}", msg, _c, _g.to_dot());
                }
            }
        }
    }

    Ok(())
}
