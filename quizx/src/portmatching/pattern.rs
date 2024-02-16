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

//! Patterns for matching.

use std::{
    collections::{HashSet, VecDeque},
    iter::FromIterator,
};

use itertools::Itertools;
use portmatching::patterns::{Edge, LinePattern};
use serde::Deserialize;

use crate::{
    flow::causal::CausalFlow,
    graph::{GraphLike, V},
    json::{JsonGraph, VertexName},
};

use super::{PEdge, PNode};

/// A causal pattern in a graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalPattern<G> {
    graph: G,
    flow: CausalFlow,
    boundary: Vec<V>,
    inputs: HashSet<V>,
    outputs: HashSet<V>,
}

impl<G: GraphLike> CausalPattern<G> {
    /// Construct a pattern from a causal graph.
    ///
    /// Alongside a graph and a causal flow, a pattern must have a designated
    /// (ordered) set of boundary vertices, and IO vertices contained in the boundary.
    ///
    /// The inputs and outputs of a pattern (given by the flow) must be in its
    /// boundary.
    pub fn new(mut graph: G, boundary: Vec<V>, inputs: HashSet<V>, outputs: HashSet<V>) -> Self {
        // Don't use quizx input/outputs. They mean something else, it's confusing
        graph.set_inputs(vec![]);
        graph.set_outputs(vec![]);

        // Compute flow
        let flow = CausalFlow::from_graph_io(&graph, &inputs, &outputs).expect("invalid flow");

        // Check that inputs and outputs are in boundary
        for v in inputs.iter().chain(&outputs) {
            assert!(boundary.contains(v));
        }

        CausalPattern {
            graph,
            flow,
            boundary,
            inputs,
            outputs,
        }
    }

    pub fn with_graph_io(graph: G, boundary: Vec<V>) -> Self {
        let inputs = graph.inputs().iter().copied().collect();
        let outputs = graph.outputs().iter().copied().collect();
        Self::new(graph, boundary, inputs, outputs)
    }

    /// Return all edges in the pattern in a valid order.
    ///
    /// A valid order is one where the source of each edge is already known
    /// (the source of the first edge is considered the root).
    pub(super) fn edges(&self) -> Vec<Edge<V, PNode, PEdge>> {
        self.to_lines()
            .iter()
            .flat_map(|(_, vec)| vec)
            .map(|&(src, dst, etype)| Edge {
                source: Some(src),
                target: Some(dst),
                edge_prop: etype,
                source_prop: Some(self.graph.vertex_type(src)),
                target_prop: Some(self.graph.vertex_type(dst)),
            })
            .collect()
    }

    pub fn root(&self) -> V {
        self.flow.lines()[self.root_line()][0]
    }

    fn root_line(&self) -> usize {
        self.flow
            .lines()
            .iter()
            .position(|line| !line.is_empty())
            .expect("empty pattern")
    }

    /// The vertices are partitioned into boundary and internal vertices.
    ///
    /// Returns an iterator over the boundary vertices.
    pub(super) fn boundary(&self) -> impl Iterator<Item = V> + ExactSizeIterator + '_ {
        self.boundary.iter().copied()
    }

    /// The vertices are partitioned into boundary and internal vertices.
    ///
    /// Returns an iterator over the internal vertices.
    pub(super) fn internal(&self) -> impl Iterator<Item = V> + '_ {
        let boundary: HashSet<_> = self.boundary().collect();
        self.graph
            .vertices()
            .filter(move |&v| !boundary.contains(&v))
    }

    /// Express a causal flow pattern as a line pattern.
    ///
    /// Line patterns are used by `portmatching` to build matcher automata.
    pub(super) fn to_line_pattern(&self) -> LinePattern<V, PNode, PEdge> {
        let mut lp = LinePattern::new();
        // Add all vertices
        for v in self.graph.vertices() {
            lp.require(v, self.graph.vertex_type(v));
        }

        // Add all edges as lines
        for (root, line) in self.to_lines() {
            lp.add_line(root, line);
        }

        lp
    }

    fn to_lines(&self) -> Vec<(V, Vec<(V, V, PEdge)>)> {
        // Lines to return
        let mut lines = Vec::new();

        // Keep track of paths to visit / visited
        let mut seen_causal_paths = HashSet::new();
        // Pick first vertex of first non-empty line as root
        let mut causal_paths_queue = VecDeque::from_iter([(self.root_line(), self.root())]);

        // Every causal path in the flow is a line in the pattern.
        // Non-causal edges each form their own line.
        // Add one causal flow at a time, choosing roots as we go so that the
        // resulting line pattern is connected.
        while let Some((path_idx, line_root)) = causal_paths_queue.pop_front() {
            // Add causal flow path as line
            seen_causal_paths.insert(path_idx);
            lines.push(self.causal_path_to_line(path_idx, line_root));

            // Discover new causal flow paths
            // Add non causal flow edges as their own lines
            let causal_path = &self.flow.lines()[path_idx];
            for &u in causal_path {
                for (v, etype) in self.graph.incident_edges(u) {
                    if self.flow.is_non_causal_edge(u, v) {
                        let new_path_idx = self.flow.line_idx(v);
                        if !seen_causal_paths.contains(&new_path_idx) {
                            // We have found a new edge to a new causal path
                            causal_paths_queue.push_back((new_path_idx, v));
                            // (u, v) is a line in itself
                            let line = vec![(u, v, PEdge::non_causal(etype))];
                            lines.push((u, line));
                        }
                    }
                }
            }
        }
        if seen_causal_paths.len() != self.flow.lines().len() {
            panic!("The pattern is not connected");
        }
        lines
    }

    fn causal_path_to_line(&self, path_idx: usize, root: V) -> (V, Vec<(V, V, PEdge)>) {
        let causal_path = &self.flow.lines()[path_idx];
        let root_idx = causal_path
            .iter()
            .position(|&v| v == root)
            .expect("root not in path");
        let left_half = causal_path[0..=root_idx]
            .iter()
            .rev()
            .tuple_windows()
            .map(|(&u, &v)| (u, v, PEdge::rev_causal(self.graph.edge_type(u, v))));
        let right_half = causal_path[root_idx..]
            .iter()
            .tuple_windows()
            .map(|(&u, &v)| (u, v, PEdge::causal(self.graph.edge_type(u, v))));
        (root, left_half.chain(right_half).collect())
    }
}

impl<G: GraphLike> From<&CausalPattern<G>> for LinePattern<V, PNode, PEdge> {
    fn from(value: &CausalPattern<G>) -> Self {
        value.to_line_pattern()
    }
}

#[derive(serde::Serialize, serde::Deserialize)]
struct JSONCausalPattern {
    graph: JsonGraph,
    boundary: Vec<VertexName>,
    inputs: HashSet<VertexName>,
    outputs: HashSet<VertexName>,
}

impl<'de, G: GraphLike + Deserialize<'de>> serde::Deserialize<'de> for CausalPattern<G> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let JSONCausalPattern {
            graph,
            boundary,
            inputs,
            outputs,
        } = JSONCausalPattern::deserialize(deserializer)?;
        let (graph, map) = graph.to_graph(true);
        let boundary = boundary.into_iter().map(|v| map[&v]).collect();
        let inputs = inputs.into_iter().map(|v| map[&v]).collect();
        let outputs = outputs.into_iter().map(|v| map[&v]).collect();
        Ok(CausalPattern::new(graph, boundary, inputs, outputs))
    }
}

impl<G: GraphLike> serde::Serialize for CausalPattern<G> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let (graph, map) = JsonGraph::from_graph_with_map(&self.graph, true);
        let boundary = self.boundary.iter().map(|v| map[v].clone()).collect();
        let inputs = self.inputs.iter().map(|v| map[v].clone()).collect();
        let outputs = self.outputs.iter().map(|v| map[v].clone()).collect();
        JSONCausalPattern {
            graph,
            boundary,
            inputs,
            outputs,
        }
        .serialize(serializer)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use crate::{
        graph::{EType, GraphLike, VType},
        vec_graph::Graph,
    };

    use super::CausalPattern;

    #[test]
    fn serialize_pattern() {
        let mut g = Graph::new();
        let i = g.add_vertex(VType::B);
        let a = g.add_vertex(VType::Z);
        let b = g.add_vertex(VType::Z);
        let c = g.add_vertex(VType::Z);
        g.add_edge_with_type(i, a, EType::N);
        g.add_edge_with_type(a, b, EType::H);
        g.add_edge_with_type(b, c, EType::H);

        g.remove_vertex(i);

        let boundary = vec![a, c];
        let inputs = HashSet::from_iter([a]);
        let outputs = HashSet::from_iter([c]);
        let pattern = CausalPattern::new(g, boundary, inputs, outputs);
        let s = serde_json::to_string(&pattern).unwrap();
        let pattern2: CausalPattern<Graph> = serde_json::from_str(&s).unwrap();
        assert_eq!(pattern.graph.num_vertices(), pattern2.graph.num_vertices());
    }
}
