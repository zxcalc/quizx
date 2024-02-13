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

//! This module defines the causal structure digraph of open graphs.
//!
//! Rewriting regions of the graph that are convex in this structure is
//! guaranteed to preserve the existence of a causal flow, as long as the
//! replacement respects the causal structure.
//!
//! Given an open graph `G = (V, E, I, O)` where `I, O \subseteq V` and a causal
//! flow `(f, ≼)` for G, the causal structure is a directed graph `C = (V, E_c)`
//! where an edge `(u, v)` with `u≠v` is in `E_c \subseteq V x V` if and only if
//! any of the following conditions hold:
//!
//! - `u ≼ f⁻¹(v)`.
//! - `u ∈ I`, and `n ≼ f⁻¹(v)` for some neighbor `n` of `u`.
//! - `v ∈ O`, and `u ≼ v`.
//! - `u ∈ I`, `v ∈ O`, and `n ≼ v` for some neighbor `n` of `u`.

use crate::graph::GraphLike;

use super::causal::CausalFlow;

/// A causal structure digraph.
///
/// Uniquely defined for a given open graph with causal flow.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CausalStructure {
    /// For each vertex `u`, a list of vertices `v` such that `u -> v` is an edge.
    edges: Vec<Vec<usize>>,
}

impl CausalStructure {
    /// Constructs the causal structure of a graph from its causal flow.
    pub fn from_flow(flow: &CausalFlow, graph: &impl GraphLike) -> Self {
        let _ = (flow, graph);
        unimplemented!("CausalStructure::from_flow");
    }
}
