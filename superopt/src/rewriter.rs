//! Rewriter for the SuperOptimizer.

use std::collections::{HashMap, HashSet};

use quizx::{
    flow::causal::CausalFlow,
    graph::GraphLike,
    portmatching::{CausalMatcher, CausalPattern, PatternID},
    vec_graph::V,
};

use crate::{
    cost::{CostDelta, CostMetric},
    rewrite_sets::{RewriteRhs, RewriteSet},
};

pub trait Rewriter {
    type Rewrite;

    /// Get the rewrites that can be applied to the graph.
    fn get_rewrites(&self, graph: &impl GraphLike) -> Vec<Self::Rewrite>;
}

pub trait Strategy<R: Rewriter> {
    type CostMetric: CostMetric;

    /// Apply the rewrites to the graph.
    fn apply_rewrites<G: GraphLike>(
        &self,
        rewrites: Vec<R::Rewrite>,
        graph: &G,
    ) -> impl Iterator<Item = RewriteResult<G>>;
}

pub struct RewriteResult<G> {
    pub graph: G,
    pub cost_delta: CostDelta,
}

#[derive(serde::Serialize, serde::Deserialize, Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct RhsIdx(usize);

/// A rewriter that applies causal flow preserving rewrites.
///
/// The set of possible rewrites are given as a list of `RewriteSet`s.
#[derive(serde::Serialize, serde::Deserialize)]
pub struct CausalRewriter<G: GraphLike> {
    matcher: CausalMatcher<G>,
    lhs_to_rhs: HashMap<PatternID, RhsIdx>,
    all_rhs: Vec<Vec<RewriteRhs<G>>>,
}

pub struct Rewrite<G> {
    /// lhs_boundary.len() == rhs_boundary.len()
    lhs_boundary: Vec<V>,
    rhs_boundary: Vec<V>,
    lhs_internal: HashSet<V>,
    rhs: G,
}

impl<G: GraphLike> Rewriter for CausalRewriter<G> {
    type Rewrite = Rewrite<G>;

    fn get_rewrites(&self, graph: &impl GraphLike) -> Vec<Self::Rewrite> {
        let flow = CausalFlow::from_graph(graph).expect("no causal flow");
        self.matcher
            .find_matches(graph, &flow)
            .flat_map(|m| {
                self.get_rhs(m.pattern_id).iter().map(move |rhs| {
                    let lhs_boundary = m.boundary.clone();
                    let lhs_internal = m.internal.clone();
                    let rhs_boundary = rhs.boundary();
                    let rhs = rhs.g.clone();
                    Rewrite {
                        lhs_boundary,
                        rhs_boundary,
                        lhs_internal,
                        rhs,
                    }
                })
            })
            .collect()
    }
}

impl<G: GraphLike + Clone> CausalRewriter<G> {
    fn get_rhs(&self, lhs_idx: PatternID) -> &[RewriteRhs<G>] {
        let idx = &self.lhs_to_rhs[&lhs_idx];
        &self.all_rhs[idx.0]
    }

    pub fn from_rewrite_rules(rules: impl IntoIterator<Item = RewriteSet<G>>) -> Self {
        let mut patterns = Vec::new();
        let mut map_to_rhs = HashMap::new();
        let mut all_rhs = Vec::new();
        for RewriteSet { lhs, lhs_ios, rhss } in rules {
            let rhs_idx = RhsIdx(all_rhs.len());
            all_rhs.push(rhss);
            let boundary = lhs.boundary();
            let lhs = lhs.g;
            for lhs_io in lhs_ios {
                let mut p = lhs.clone();
                p.set_inputs(lhs_io.inputs().to_owned());
                p.set_outputs(lhs_io.outputs().to_owned());
                let flow = CausalFlow::from_graph(&p).expect("invalid causal flow in pattern");
                patterns.push(CausalPattern::new(p, flow, boundary.clone()));
                map_to_rhs.insert(PatternID(patterns.len() - 1), rhs_idx);
            }
        }
        CausalRewriter {
            matcher: CausalMatcher::from_patterns(patterns),
            lhs_to_rhs: map_to_rhs,
            all_rhs,
        }
    }
}
