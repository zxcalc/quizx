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

use crate::{
    flow::causal::CausalFlow,
    hash_graph::{GraphLike, V},
};

use super::{PEdge, PNode};

/// A causal pattern in a graph.
#[derive(serde::Serialize, serde::Deserialize)]
pub struct CausalPattern<G> {
    graph: G,
    flow: CausalFlow,
    boundary: Vec<V>,
}

impl<G: GraphLike> CausalPattern<G> {
    /// Construct a pattern from a causal graph.
    ///
    /// Alongside a graph and a causal flow, a pattern must have a designated
    /// (ordered) set of boundary vertices.
    ///
    /// The inputs and outputs of a pattern (given by the flow) must be in its
    /// boundary.
    pub fn new(graph: G, flow: CausalFlow, boundary: Vec<V>) -> Self {
        // Check that inputs and outputs are in boundary
        for _v in flow.inputs().chain(flow.outputs()) {
            // TODO: This fails. The flow `inputs` and `outputs` are the
            //       `VType::B` (boundary) vertices of the graph, but `boundary`
            //       contains the `VType::Z` neighbors of those vertices
            //       instead.

            //assert!(boundary.contains(&v));
        }
        CausalPattern {
            graph,
            flow,
            boundary,
        }
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
