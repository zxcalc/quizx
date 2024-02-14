// QuiZX - Rust library for quantum circuit rewriting and optimization
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

//! Utilities for computing the cost of a graph.

use quizx::graph::GraphLike;

pub type Cost = usize;
pub type CostDelta = isize;

/// A trait for computing the cost of a graph.
pub trait CostMetric {
    /// Create a new instance of the cost metric.
    fn new() -> Self;

    /// Compute the cost of a graph.
    fn cost(&self, graph: &impl GraphLike) -> Cost;
}

/// A cost metric that counts the number of edges in a graph.
pub struct EdgeCount;

impl CostMetric for EdgeCount {
    fn new() -> Self {
        Self
    }

    fn cost(&self, graph: &impl GraphLike) -> Cost {
        graph.num_edges()
    }
}

/// A cost metric that counts the number of vertices in a graph.
pub struct VertexCount;

impl CostMetric for VertexCount {
    fn new() -> Self {
        Self
    }

    fn cost(&self, graph: &impl GraphLike) -> Cost {
        graph.num_vertices()
    }
}

/// For a graph with causal flow, the 2-qubit gate count corresponds to the number of
/// edges in the graph minus the number of vertices, plus the number of inputs.
pub struct TwoQubitGateCount;

impl CostMetric for TwoQubitGateCount {
    fn new() -> Self {
        Self
    }

    fn cost(&self, graph: &impl GraphLike) -> Cost {
        graph.inputs().len() + graph.num_edges() - graph.num_vertices()
    }
}
