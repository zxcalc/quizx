//! Rewriter for the SuperOptimizer.

use quizx::graph::GraphLike;

use crate::cost::{CostDelta, CostMetric};

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
