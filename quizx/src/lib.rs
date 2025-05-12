//! [PyZX](https://github.com/zxlang/pyzx) is a Python library for quantum circuit optimisation and compiling using the [ZX-calculus](https://zxcalculus.com). It's great for hacking, learning, and trying things out in [Jupyter](https://jupyter.org/) notebooks. However, it's written to maximise clarity and fun, not performance.
//!
//! This is a port of some of the core functionality of PyZX to the [Rust](https://www.rust-lang.org/) programming language. This is a modern systems programming language, which enables writing software that is very fast and memory efficient.
//!
//! Check the [Rust Changelog](https://github.com/zxcalc/quizx/blob/master/quizx/CHANGELOG.md) for the latest updates.
//!
//! ## A bit about performance
//!
//! As a very anecdotal example of the performance difference, the program `spider_chain` builds a chain of 1 million green spiders and fuses them all. In PyZX, you can fuse all the spiders in a ZX-diagram as follows:
//!
//! ```python
//! from pyzx.basicrules import *
//!
//! success = True
//! while success:
//!     success = any(fuse(g, g.edge_s(e), g.edge_t(e)) for e in g.edges())
//! ```
//!
//! In QuiZX, the Rust code is slightly more verbose, but similar in spirit:
//! ```rust
//! # use quizx::graph::*;
//! # use quizx::vec_graph::Graph;
//! # use quizx::basic_rules::check_spider_fusion;
//! use quizx::basic_rules::*;
//!
//! let mut g = Graph::new();
//! let v0 = g.add_vertex(VType::Z);
//! let v1 = g.add_vertex(VType::Z);
//! let v2 = g.add_vertex(VType::X);
//! g.add_edge(v0, v1);
//! g.add_edge(v1, v2);
//!
//! loop {
//!     match g.find_edge(|v0,v1,_| check_spider_fusion(&g, v0, v1)) {
//!         Some((v0,v1,_)) => spider_fusion_unchecked(&mut g, v0, v1),
//!         None => break,
//!     };
//! }
//! ```
//!
//! On my laptop, the PyZX code takes about 98 seconds to fuse 1 million spiders, whereas the QuiZX code takes 17 milliseconds.

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

pub mod annealer;
pub mod basic_rules;
pub mod circuit;
pub mod decompose;
pub mod equality;
pub mod extract;
pub mod fscalar;
pub mod gate;
pub mod generate;
pub mod graph;
pub mod hash_graph;
pub mod json;
pub mod linalg;
pub mod optimize_circuit;
pub mod params;
pub mod phase;
pub mod random_graph;
// pub mod scalar;
pub mod scalar_traits;
pub mod simplify;
pub mod tensor;
pub mod util;
pub mod vec_graph;
