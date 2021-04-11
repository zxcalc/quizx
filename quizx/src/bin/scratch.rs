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
use quizx::scalar::*;
use num::Rational;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1/(2+2j) = (2-2j)/8
    let mut global = Scalar4::Exact(-3, [2, 0, -2, 0]);
    // times -(7 + 5*rt2)
    global *= Scalar4::Exact(0, [-7, -5, 0, -5]);

    let b60 = Scalar4::Exact(-3, [-16, 12, 0, 12]) * &global;
    let b66 = Scalar4::Exact(-3, [-96, 68, 0, 68]) * &global;
    let o6 = Scalar4::Exact(1, [0, 10, -14, 10]) * &global;
    let mut k6 = Scalar4::Exact(0, [7, -5, 0, -5]) * &global;
    k6.mul_sqrt2_pow(5);
    let mut phi = Scalar4::Exact(0, [10, -7, 0, -7]) * &global;
    phi.mul_sqrt2_pow(9);
    phi.mul_phase(Rational::new(3,2));

    println!("let scalar = ScalarN::{:?}; // b60", b60);
    println!("let scalar = ScalarN::{:?}; // b66", b66);
    println!("let scalar = ScalarN::{:?}; // o6", o6);
    println!("let scalar = ScalarN::{:?}; // k6", k6);
    println!("let scalar = ScalarN::{:?}; // phi", phi);

    Ok(())
}
