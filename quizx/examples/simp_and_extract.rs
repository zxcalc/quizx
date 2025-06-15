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
use quizx::extract::*;
use quizx::simplify::*;
use quizx::tensor::*;
use quizx::vec_graph::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Circuit::from_file("circuits/small/mod5_4.qasm")?;
    let c0 = c.to_basic_gates();
    println!("Before:\n{}", c0.stats());
    // let c = c.to_basic_gates();
    // let c1 = c.to_basic_gates();
    // if !Tensor4::scalar_compare(&c, &c1) {
    //     panic!("Tensors don't match: c, c1");
    // }

    // let c = Circuit::random()
    //     .seed(1337)
    //     .qubits(5)
    //     .depth(30)
    //     .p_t(0.2)
    //     .with_cliffords()
    //     .build();
    let mut g: Graph = c.to_graph();
    full_simp(&mut g);
    // println!("{}", g.to_dot());
    // println!("{:?}", g);
    // assert_eq!(c.to_tensor4(), g.to_tensor4());

    match g.to_circuit() {
        Ok(c1) => {
            println!("extracted ok");
            if Tensor4::scalar_compare(&c, &c1) {
                println!("Tensors match!");
                println!("After:\n{}", c1.stats());
            } else {
                println!("Tensors don't match. \n{}\n\n{}", c, c1);
            }
        }
        Err(ExtractError(msg, _c, _g)) => {
            println!("extract failed: {}", msg);
            println!("{}\n\n{}\n\n{}", msg, _c, _g.to_dot());
        }
    }
    Ok(())
}
