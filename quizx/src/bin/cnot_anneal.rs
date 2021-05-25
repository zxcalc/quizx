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
use quizx::hash_graph::*;
use quizx::simplify::*;
use quizx::extract::*;
use quizx::annealer::*;

use std::time::Instant;
use std::{thread, time};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let f = "../circuits/small/tof_10.qasm";
    let time = Instant::now();
    println!("{}", f);
    let c = Circuit::from_file(f)
        .expect(&format!("circuit failed to parse: {}", f))
        .to_basic_gates();
    println!("...done reading in {:.2?}", time.elapsed());

    println!("Before: {}", c.stats());

    println!("Simplifying circuit...");
    let time = Instant::now();
    let mut g: Graph = c.to_graph();
    flow_simp(&mut g);
    println!("Done simplifying in {:.2?}", time.elapsed());

    println!("Annealing...");
    let time = Instant::now();
    let mut annealer = Annealer::new(g);
    annealer
        .seed(1337)
        .temp(25.0)
        .anneal();
    g = annealer.g;
    println!("Done annealing in {:.2?}", time.elapsed());

    println!("Extracting circuit...");
    let time = Instant::now();
    let result = g.extractor()
        .gflow()
        // .up_to_perm()
        .extract();


    match result {
        Ok(_c1) => {
            println!("Done in {:.2?}", time.elapsed());
            println!("extracted ok");
            println!("After: {}", _c1.stats());
        },
        Err(ExtractError(msg, _c, _g)) => {
            println!("extract failed: {}", msg);
            // println!("\n\n{}", _g.to_dot());
        },
    }


    thread::sleep(time::Duration::from_millis(100));
    Ok(())
}
