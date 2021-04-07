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
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for e in fs::read_dir("../../circuits")? {
        if let Some(f) = e?.path().to_str() {
            let time = Instant::now();
            println!("{}", f);
            Circuit::from_file(f).expect(&format!("circuit failed to parse: {}", f));
            println!("...done in {:.2?}", time.elapsed());
        }
    }

    Ok(())
}
