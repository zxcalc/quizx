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

use crate::hash_graph::{EType, VType};
use portmatching::EdgeProperty;

// TODO: choose correct types
type PNode = VType;

/// An edge to be matched in a causal graph
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct PEdge {
    src: CausalPort,
    dst: CausalPort,
    etype: EType,
}

/// The types of port in causal ZX graphs.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

    fn property_id(&self) -> Option<Self::OffsetID> {
        match self.src {
            p @ (CausalPort::CausalInput | CausalPort::CausalOutput) => Some(p),
            CausalPort::Other => None,
        }
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
    use std::fs;

    use super::*;
    use crate::flow::causal::CausalFlow;
    use crate::hash_graph::V;
    use crate::vec_graph::Graph;
    use crate::{graph::VType, hash_graph::GraphLike};

    use itertools::Itertools;
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
        let p = CausalPattern::new(g.clone(), flow.clone(), boundary);

        let matcher = CausalMatcher::from_patterns(vec![p]);

        let res = matcher.find_matches(&g, &flow).collect_vec();
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].boundary, vec![vs[0], vs[1], vs[6], vs[7]]);
        assert_eq!(res[1].boundary, vec![vs[1], vs[0], vs[7], vs[6]]);
        assert_eq!(
            res[0].internal,
            vec![vs[2], vs[3], vs[4], vs[5]]
                .into_iter()
                .collect::<HashSet<_>>()
        );
    }
}
