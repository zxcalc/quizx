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

//! Methods for computing the causal flow of a graph.
//!
//! Given an open graph `G = (V, E, I, O)` where `I, O \subseteq V` are
//! respectively the input and output node sets, a causal flow is a function
//! `f: (V / O) -> (V / I)` and a partial order `≼` such that:
//!
//! 1. `f` is injective
//! 2. `u ~ f(u)`
//! 3. `u ≼ f(u)`
//! 4. If `v ~ f(u)`, then `u ≼ v`

use std::collections::HashSet;
use std::ops::Range;

use itertools::Itertools;

use crate::graph::{GraphLike, V};

/// A causal flow of a graph.
///
/// Note that we only keep track of the flow function `f` and not the partial
/// order `≼`.
//
// TODO: Store the order too? Perhaps as an `Option`, and compute it on demand otherwise?
#[derive(serde::Serialize, serde::Deserialize, Debug, Clone, PartialEq, Eq)]
pub struct CausalFlow {
    /// The flow lines in the graph.
    ///
    /// The flow function `f` maps each line element to the next element in the vector.
    /// The last element of each line is an output node.
    lines: Vec<Vec<V>>,
    /// For each vertex, the line it is part of and its position in the line.
    positions: Vec<FlowPosition>,
}

/// The position of a vertex in a set of flow lines.
#[derive(
    serde::Serialize, serde::Deserialize, Debug, Copy, Clone, PartialEq, Eq, derive_more::From,
)]
pub struct FlowPosition {
    pub line: usize,
    pub pos: usize,
}

impl CausalFlow {
    pub fn from_graph(g: &impl GraphLike) -> Result<Self, CausalFlowError> {
        let inputs = g.inputs().iter().copied().collect();
        let outputs = g.outputs().iter().copied().collect();
        Self::from_graph_io(g, &inputs, &outputs)
    }
    /// Computes the causal flow of a graph.
    pub fn from_graph_io(
        g: &impl GraphLike,
        inputs: &HashSet<V>,
        outputs: &HashSet<V>,
    ) -> Result<Self, CausalFlowError> {
        // Sweep the graph from the inputs, consuming a candidate node with a
        // single non-candidate neighbour.

        // The candidates.
        let mut candidates: HashSet<V> = inputs.iter().copied().collect();

        let mut lines: Vec<Vec<V>> = Vec::with_capacity(outputs.len());
        let max_ind = g.vertices().max().unwrap_or(0);
        let mut positions: Vec<FlowPosition> = vec![(usize::MAX, usize::MAX).into(); max_ind + 1];
        let is_visited =
            |v: V, positions: &Vec<FlowPosition>| positions[v] != (usize::MAX, usize::MAX).into();

        let mut visited = 0;

        while !candidates.is_empty() {
            // Find a candidate node with a single non-candidate neighbor.
            let Some((candidate, neigh)) = candidates
                .iter()
                .filter_map(|&v| {
                    let single_neighbour = g
                        .neighbors(v)
                        .filter(|&n| !is_visited(n, &positions) && !candidates.contains(&n))
                        .exactly_one()
                        .ok()?;
                    Some((v, single_neighbour))
                })
                .next()
            else {
                return Err(CausalFlowError::NonCausal);
            };

            // Remove the candidate from the set, and add the neighbor (unless
            // it is an output).
            visited += 1;

            // If the candidate does not have a line yet, create a new line (it is an input).
            let (line, pos) = match positions[candidate] {
                FlowPosition {
                    line: usize::MAX,
                    pos: usize::MAX,
                } => {
                    let line = lines.len();
                    lines.push(vec![candidate]);
                    positions[candidate] = (line, 0).into();
                    (line, 1)
                }
                p => (p.line, p.pos + 1),
            };
            debug_assert_eq!(lines[line].len(), pos);
            lines[line].push(neigh);
            positions[neigh] = (line, pos).into();

            candidates.remove(&candidate);
            if !outputs.contains(&neigh) {
                candidates.insert(neigh);
            } else {
                visited += 1
            }
        }

        if visited != g.num_vertices() {
            return Err(CausalFlowError::NonCausal);
        }

        Ok(Self { lines, positions })
    }

    /// Returns the next node in the causal flow.
    pub fn next(&self, v: V) -> Option<V> {
        let p = self.positions[v];
        self.lines[p.line].get(p.pos + 1).copied()
    }

    /// Returns the previous node in the causal flow.
    pub fn pred(&self, v: V) -> Option<V> {
        let p = self.positions[v];
        p.pos.checked_sub(1).map(|pos| self.lines[p.line][pos])
    }

    pub fn inputs(&self) -> impl Iterator<Item = V> + '_ {
        self.lines.iter().map(|line| line[0])
    }

    pub fn outputs(&self) -> impl Iterator<Item = V> + '_ {
        self.lines.iter().map(|line| *line.last().unwrap())
    }

    /// The line index of a vertex.
    pub fn line_idx(&self, v: V) -> usize {
        self.positions[v].line
    }

    /// Returns the causal flow lines of the graph.
    pub fn lines(&self) -> &[Vec<V>] {
        &self.lines
    }

    /// Whether u -> v is a causal edge.
    pub fn is_causal_edge_dir(&self, u: V, v: V) -> bool {
        self.next(u) == Some(v)
    }

    /// Whether (u, v) is a causal edge.
    ///
    /// Either u -> v or v -> u is a causal edge.
    pub fn is_causal_edge(&self, u: V, v: V) -> bool {
        self.is_causal_edge_dir(u, v) || self.is_causal_edge_dir(v, u)
    }

    /// Whether (u, v) is a non-causal edge.
    ///
    /// This is true when both u -> v and v -> u are not causal.
    pub fn is_non_causal_edge(&self, u: V, v: V) -> bool {
        !self.is_causal_edge(u, v)
    }
}

/// The convex hull of a graph region, with respect to its causal flow.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConvexHull {
    /// The original vertices of the region.
    pub region: HashSet<V>,
    /// The additional vertices that are part of the convex hull.
    pub hull_vertices: HashSet<V>,
}

impl ConvexHull {
    /// Constructs the convex hull of a region given the graph's causal flow.
    pub fn from_region<'a>(region: impl IntoIterator<Item = &'a V>, flow: &CausalFlow) -> Self {
        let flow_lines = flow.lines();

        // For each flow line, the first and last vertices that are part of the region.
        let mut line_min_max: Vec<Option<Range<usize>>> = vec![None; flow_lines.len()];

        let region: HashSet<V> = region.into_iter().copied().collect();

        for &v in &region {
            let p = flow.positions[v];
            let min_max = line_min_max[p.line].get_or_insert_with(|| p.pos..p.pos);
            min_max.start = min_max.start.min(p.pos);
            min_max.end = min_max.end.max(p.pos);
        }

        let mut hull_vertices = HashSet::new();
        for (line, min_max) in line_min_max.iter().enumerate() {
            if let Some(range) = min_max {
                let line = &flow_lines[line];
                for v in line[range.clone()].iter().copied() {
                    if !region.contains(&v) {
                        hull_vertices.insert(v);
                    }
                }
            }
        }

        Self {
            region,
            hull_vertices,
        }
    }

    /// Returns all the vertices that are part of the region.
    pub fn vertices(&self) -> impl Iterator<Item = V> + '_ {
        self.region
            .iter()
            .copied()
            .chain(self.hull_vertices.iter().copied())
    }
}

#[derive(Debug, Clone, thiserror::Error)]
pub enum CausalFlowError {
    /// The graph does not have a causal flow.
    #[error("The graph does not have a causal flow")]
    NonCausal,
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::graph::VType;
    use crate::vec_graph::Graph;

    use rstest::{fixture, rstest};

    /// Makes a simple graph with a causal flow.
    ///
    /// The graph is:
    /// ```text
    /// 0 - 2 - 4 - 6
    ///     |   |
    /// 1 - 3 - 5 - 7
    /// ```
    ///
    /// With `0` and `1` as inputs, and `6` and `7` as outputs.
    #[fixture]
    fn simple_graph() -> (Graph, Vec<V>) {
        let mut g = Graph::new();
        let vs = vec![
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::X),
            g.add_vertex(VType::X),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];

        g.set_inputs(vec![vs[0], vs[1]]);
        g.set_outputs(vec![vs[6], vs[7]]);

        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[3]);
        g.add_edge(vs[2], vs[4]);
        g.add_edge(vs[2], vs[3]);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[4], vs[6]);
        g.add_edge(vs[5], vs[7]);
        (g, vs)
    }

    #[rstest]
    fn causal_flow(simple_graph: (Graph, Vec<V>)) {
        let (g, vs) = simple_graph;
        let flow = CausalFlow::from_graph(&g).unwrap();
        assert_eq!(flow.next(vs[0]), Some(vs[2]));
        assert_eq!(flow.next(vs[2]), Some(vs[4]));
        assert_eq!(flow.next(vs[4]), Some(vs[6]));
        assert_eq!(flow.next(vs[6]), None);
        assert_eq!(flow.next(vs[1]), Some(vs[3]));
        assert_eq!(flow.next(vs[3]), Some(vs[5]));
        assert_eq!(flow.next(vs[5]), Some(vs[7]));
        assert_eq!(flow.next(vs[7]), None);
    }

    #[rstest]
    fn convex_hull(simple_graph: (Graph, Vec<V>)) {
        let (g, vs) = simple_graph;
        let flow = CausalFlow::from_graph(&g).unwrap();

        let region = vec![vs[0], vs[3], vs[6], vs[4]];
        let expected_other = [vs[2]];

        let hull = ConvexHull::from_region(&region, &flow);

        assert_eq!(hull.region, region.iter().copied().collect());
        assert_eq!(hull.hull_vertices, expected_other.iter().copied().collect());
    }
}
