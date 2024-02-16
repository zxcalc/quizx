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

//! Pattern matching is done using pre-built matcher objects.

use std::collections::HashSet;

use itertools::Itertools;
use portmatching::{
    automaton::{LineBuilder, ScopeAutomaton},
    PatternID, SinglePatternMatcher,
};

use crate::{
    flow::causal::CausalFlow,
    graph::{GraphLike, V},
};

use super::{CausalPattern, CausalPort, PEdge, PNode};

/// A pre-built pattern matcher for a set of patterns.
#[derive(serde::Serialize, serde::Deserialize)]
pub struct CausalMatcher<G: GraphLike> {
    automaton: ScopeAutomaton<PNode, PEdge, CausalPort>,
    patterns: Vec<CausalPattern<G>>,
}

impl<G: GraphLike> CausalMatcher<G> {
    /// Find all matches of the patterns in a causal graph.
    pub fn find_matches<'s, 'g: 's>(
        &'s self,
        graph: &'g impl GraphLike,
        flow: &'g CausalFlow,
    ) -> impl Iterator<Item = PatternMatch> + 's {
        graph
            .vertices()
            .flat_map(move |root| self.find_rooted_matches(root, graph, flow))
    }

    pub fn find_rooted_matches<'s, 'g: 's>(
        &'s self,
        root: V,
        graph: &'g impl GraphLike,
        flow: &'g CausalFlow,
    ) -> impl Iterator<Item = PatternMatch> + 's {
        self.run_automaton(root, graph, flow)
            .unique()
            .flat_map(move |pattern_id| self.get_pattern_matches(pattern_id, root, graph, flow))
    }

    pub(super) fn run_automaton<'s, 'g: 's>(
        &'s self,
        root: V,
        graph: &'g impl GraphLike,
        flow: &'g CausalFlow,
    ) -> impl Iterator<Item = PatternID> + 's {
        self.automaton
            .run(root, vertex_predicate(graph), edge_predicate(graph, flow))
    }

    /// Build a matcher from a set of patterns.
    pub fn from_patterns(patterns: Vec<CausalPattern<G>>) -> Self {
        let lps = patterns.iter().map_into().collect();
        let builder = LineBuilder::from_patterns(lps);
        Self {
            automaton: builder.build(),
            patterns,
        }
    }

    /// Get a pattern by its ID.
    pub fn get_pattern(&self, pattern_id: PatternID) -> &CausalPattern<G> {
        &self.patterns[pattern_id.0]
    }

    pub fn dot_string(&self) -> String {
        self.automaton.dot_string()
    }

    pub(super) fn get_pattern_matches(
        &self,
        pattern_id: PatternID,
        root: V,
        graph: &impl GraphLike,
        flow: &CausalFlow,
    ) -> Vec<PatternMatch> {
        let pattern = self.get_pattern(pattern_id);
        let matches = SinglePatternMatcher::new(pattern, pattern.edges(), pattern.root())
            .get_match_map(root, vertex_predicate(graph), edge_predicate(graph, flow));
        matches
            .into_iter()
            .map(|match_map| {
                let boundary = pattern
                    .boundary()
                    .map(|v| *match_map.get_by_left(&v).unwrap())
                    .collect();
                let internal = pattern
                    .internal()
                    .map(|v| *match_map.get_by_left(&v).unwrap())
                    .collect();
                PatternMatch {
                    pattern_id,
                    boundary,
                    internal,
                }
            })
            .collect()
    }
}

#[derive(Clone, Debug)]
pub struct PatternMatch {
    /// The ID of the pattern
    pub pattern_id: PatternID,
    /// The list of vertices on the boundary
    pub boundary: Vec<V>,
    /// The set of internal vertices
    pub internal: HashSet<V>,
}

fn edge_predicate<'g>(
    graph: &'g impl GraphLike,
    flow: &'g CausalFlow,
) -> impl for<'a> Fn(V, &'a PEdge) -> Vec<Option<V>> + 'g {
    move |v, &PEdge { src, dst, etype }| {
        let is_correct_port = |v, neigh, port| match port {
            CausalPort::CausalInput => flow.is_causal_edge_dir(neigh, v),
            CausalPort::CausalOutput => flow.is_causal_edge_dir(v, neigh),
            CausalPort::Other => flow.is_non_causal_edge(v, neigh),
        };
        graph
            .incident_edges(v)
            .map(|(neigh, t)| {
                (is_correct_port(v, neigh, src) && is_correct_port(neigh, v, dst) && t == etype)
                    .then_some(neigh)
            })
            .collect()
    }
}

// TODO: currently not checking anything for vertices
fn vertex_predicate(graph: &impl GraphLike) -> impl for<'a> Fn(V, &'a PNode) -> bool + '_ {
    move |v, vtype| graph.vertex_type(v) == *vtype
}
