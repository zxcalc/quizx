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

//! Pattern Matching in Causal ZX Graphs
//!
//! This uses the `portmatching` crate

pub mod matcher;
pub mod pattern;

pub use matcher::CausalMatcher;
pub use pattern::CausalPattern;
pub use portmatching::PatternID;

use crate::graph::{EType, VType};
use portmatching::EdgeProperty;

// TODO: choose correct types
type PNode = VType;

/// An edge to be matched in a causal graph
#[derive(
    serde::Serialize, serde::Deserialize, Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash,
)]
struct PEdge {
    src: CausalPort,
    dst: CausalPort,
    etype: EType,
}

/// The types of port in causal ZX graphs.
#[derive(
    serde::Serialize, serde::Deserialize, Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash,
)]
enum CausalPort {
    /// The causal edge to the vertex predecessor.
    CausalInput,
    /// The causal edge to the vertex successor.
    CausalOutput,
    /// Any other edge
    Other,
}

impl EdgeProperty for PEdge {
    type OffsetID = CausalPort;

    fn reverse(&self) -> Option<Self> {
        Some(PEdge {
            src: self.dst,
            dst: self.src,
            etype: self.etype,
        })
    }

    fn offset_id(&self) -> Self::OffsetID {
        self.src
    }
}

impl PEdge {
    fn causal(etype: EType) -> Self {
        PEdge {
            src: CausalPort::CausalOutput,
            dst: CausalPort::CausalInput,
            etype,
        }
    }

    fn rev_causal(etype: EType) -> Self {
        PEdge {
            src: CausalPort::CausalInput,
            dst: CausalPort::CausalOutput,
            etype,
        }
    }

    fn non_causal(etype: EType) -> Self {
        PEdge {
            src: CausalPort::Other,
            dst: CausalPort::Other,
            etype,
        }
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use super::*;
    use crate::flow::causal::CausalFlow;
    use crate::graph::{GraphLike, VType, V};
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
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
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

    /// Makes a simple graph with non-unique edges.
    ///
    /// The graph is:
    /// ```text
    /// 0 --- 2 - 5
    ///      / \
    /// 1 - 3 - 4 - 6
    /// ```
    ///
    /// With `0` and `1` as inputs, and `5` and `6` as outputs.
    #[fixture]
    fn simple_multi_graph() -> (Graph, Vec<V>) {
        let mut g = Graph::new();
        let vs = vec![
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];

        g.set_inputs(vec![vs[0], vs[1]]);
        g.set_outputs(vec![vs[5], vs[6]]);

        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[2], vs[5]);
        g.add_edge(vs[1], vs[3]);
        g.add_edge(vs[3], vs[4]);
        g.add_edge(vs[4], vs[6]);
        g.add_edge(vs[2], vs[3]);
        g.add_edge(vs[2], vs[4]);
        (g, vs)
    }

    /// Makes a simple pattern that is a subgraph of `simple_multi_graph`.
    ///
    /// The graph is:
    /// ```text
    /// 0 --- 2
    ///      /
    ///     3 - 4
    /// ```
    #[fixture]
    fn p1() -> CausalPattern<Graph> {
        let mut p = Graph::new();
        let vs = vec![
            p.add_vertex(VType::B),
            p.add_vertex(VType::Z),
            p.add_vertex(VType::Z),
            p.add_vertex(VType::Z),
        ];

        p.set_inputs(vec![vs[0], vs[2]]);
        p.set_outputs(vec![vs[1], vs[3]]);

        p.add_edge(vs[0], vs[1]);
        p.add_edge(vs[1], vs[2]);
        p.add_edge(vs[2], vs[3]);
        CausalPattern::with_graph_io(p, vec![vs[0], vs[1], vs[2], vs[3]])
    }

    /// Makes a simple pattern that is a subgraph of `simple_multi_graph`.
    ///
    /// The graph is:
    /// ```text
    /// 0 --- 2
    ///        \
    ///         4 - 6
    /// ```
    #[fixture]
    fn p2() -> CausalPattern<Graph> {
        let mut p = Graph::new();
        let vs = vec![
            p.add_vertex(VType::B),
            p.add_vertex(VType::Z),
            p.add_vertex(VType::Z),
            p.add_vertex(VType::B),
        ];

        p.set_inputs(vec![vs[0], vs[2]]);
        p.set_outputs(vec![vs[1], vs[3]]);

        p.add_edge(vs[0], vs[1]);
        p.add_edge(vs[1], vs[2]);
        p.add_edge(vs[2], vs[3]);
        CausalPattern::with_graph_io(p, vec![vs[0], vs[1], vs[2], vs[3]])
    }

    #[rstest]
    fn pattern_match(simple_graph: (Graph, Vec<V>)) {
        let (g, vs) = simple_graph;
        let flow = CausalFlow::from_graph(&g).unwrap();
        let boundary = g
            .inputs()
            .iter()
            .copied()
            .chain(g.outputs().iter().copied())
            .collect();
        let p = CausalPattern::with_graph_io(g.clone(), boundary);

        let matcher = CausalMatcher::from_patterns(vec![p]);

        let res: HashSet<_> = matcher
            .find_matches(&g, &flow)
            .map(|m| {
                assert_eq!(
                    m.internal,
                    vec![vs[2], vs[3], vs[4], vs[5]]
                        .into_iter()
                        .collect::<HashSet<_>>()
                );
                m.boundary
            })
            .collect();
        assert_eq!(res.len(), 2);
        assert!(res.contains(&vec![vs[0], vs[1], vs[6], vs[7]]));
        assert!(res.contains(&vec![vs[1], vs[0], vs[7], vs[6]]));
    }

    #[rstest]
    fn pattern_match2(
        simple_multi_graph: (Graph, Vec<V>),
        p1: CausalPattern<Graph>,
        p2: CausalPattern<Graph>,
    ) {
        let (g, _) = simple_multi_graph;
        let flow = CausalFlow::from_graph(&g).unwrap();

        let matcher = CausalMatcher::from_patterns(vec![p1, p2]);

        let res: HashSet<_> = matcher
            .find_matches(&g, &flow)
            .map(|m| m.boundary)
            .collect();
        assert_eq!(res.len(), 3);
        assert!(res.contains(&vec![0, 2, 3, 4]));
        assert!(res.contains(&vec![0, 2, 4, 6]));
        assert!(res.contains(&vec![1, 3, 2, 5]));
    }
}
