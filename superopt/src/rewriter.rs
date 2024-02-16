//! Rewriter for the SuperOptimizer.

use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use quizx::vec_graph::{EType, VType};
use quizx::{
    flow::causal::{CausalFlow, ConvexHull},
    graph::GraphLike,
    portmatching::{CausalMatcher, CausalPattern, PatternID},
    vec_graph::V,
};

use crate::rewrite_sets::RuleSide;
use crate::{
    cost::CostDelta,
    rewrite_sets::{RewriteRhs, RewriteSet},
};

pub trait Rewriter {
    type Rewrite;

    /// Get the rewrites that can be applied to the graph.
    fn get_rewrites(&self, graph: &impl GraphLike) -> Vec<Self::Rewrite>;

    /// Apply the rewrites to the graph.
    fn apply_rewrite<G: GraphLike>(&self, rewrite: Self::Rewrite, graph: &G) -> RewriteResult<G>;
}

pub struct RewriteResult<G> {
    pub graph: G,
    pub cost_delta: CostDelta,
}

#[derive(serde::Serialize, serde::Deserialize, Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct RhsIdx(usize);

/// A rewriter that applies causal flow preserving rewrites.
///
/// The set of possible rewrite rule are given as a list of `RewriteSet`s.
#[derive(serde::Serialize, serde::Deserialize)]
pub struct CausalRewriter<G: GraphLike> {
    matcher: CausalMatcher<G>,
    lhs_to_rhs: HashMap<PatternID, RhsIdx>,
    all_rhs: Vec<Vec<RewriteRhs<G>>>,
}

#[derive(Clone, Debug)]
pub struct Rewrite<G> {
    /// The nodes matching the LHS boundary in the matched graph.
    lhs_boundary: Vec<V>,
    /// The nodes matching the RHS boundary in `rhs`.
    rhs_boundary: Vec<V>,
    /// The internal nodes of the LHS in the matched graph.
    lhs_internal: HashSet<V>,
    /// The replacement graph.
    rhs: G,
    /// The cost delta of the rewrite.
    ///
    /// Negative delta is an improvement (cost decrease).
    cost_delta: CostDelta,
}

impl<G: GraphLike> Rewrite<G> {
    fn lhs_vertices(&self) -> impl Iterator<Item = &V> + '_ {
        self.lhs_boundary.iter().chain(&self.lhs_internal)
    }

    /// Whether the rewrite is flow preserving when applied on `graph`.
    ///
    /// TODO: This can be done faster by pre-computing a "causal structure"
    fn is_flow_preserving(&self, graph: &impl GraphLike, flow: &CausalFlow) -> bool {
        let hull = ConvexHull::from_region(self.lhs_vertices(), flow);
        let mut subgraph = graph.induced_subgraph(hull.vertices());
        subgraph.set_inputs(hull.inputs().to_owned());
        subgraph.set_outputs(hull.outputs().to_owned());
        self.apply(&mut subgraph);
        CausalFlow::from_graph(&subgraph).is_ok()
    }

    fn apply<H: GraphLike>(&self, g: &mut H) -> CostDelta {
        let mut new_r_names: HashMap<V, V> = HashMap::new();

        // Remove the internal nodes of the LHS.
        for &v in &self.lhs_internal {
            g.remove_vertex(v);
        }

        // Replace the LHS boundary nodes with the RHS's.
        for (&l, &r) in self.lhs_boundary.iter().zip(self.rhs_boundary.iter()) {
            new_r_names.insert(r, l);
            g.set_phase(l, self.rhs.phase(r));
            g.set_vertex_type(l, self.rhs.vertex_type(r));
        }

        // Insert the internal nodes of the RHS.
        for r in self.rhs.vertices() {
            if new_r_names.contains_key(&r) {
                // It was already added as a boundary node.
                continue;
            }

            let vtype = self.rhs.vertex_type(r);
            if vtype == VType::B {
                continue;
            }

            let l = g.add_vertex_with_phase(vtype, self.rhs.phase(r));
            new_r_names.insert(r, l);
        }

        // Reconnect the edges.
        for (u, v, ty) in self.rhs.edges() {
            let (Some(&u), Some(&v)) = (new_r_names.get(&u), new_r_names.get(&v)) else {
                // Ignore the boundary edges.
                continue;
            };
            assert_eq!(ty, EType::H);
            g.add_edge_smart(u, v, ty);
        }

        self.cost_delta
    }
}

impl<G: GraphLike> Rewriter for CausalRewriter<G> {
    type Rewrite = Rewrite<G>;

    fn get_rewrites(&self, graph: &impl GraphLike) -> Vec<Self::Rewrite> {
        let flow = CausalFlow::from_graph(graph).expect("no causal flow");
        let mut rewrites = self
            .matcher
            .find_matches(graph, &flow)
            .flat_map(|m| {
                self.get_rhs(m.pattern_id).iter().map(move |rhs| {
                    let lhs_boundary = m.boundary.clone();
                    let lhs_internal = m.internal.clone();
                    let rhs_boundary = rhs.boundary().collect_vec();
                    let cost_delta = -rhs.reduction;
                    let rhs = rhs.graph();
                    assert_eq!(lhs_boundary.len(), rhs_boundary.len());
                    Rewrite {
                        lhs_boundary,
                        rhs_boundary,
                        lhs_internal,
                        rhs,
                        cost_delta,
                    }
                })
            })
            .filter(|rw| rw.is_flow_preserving(graph, &flow))
            .collect_vec();
        rewrites.sort_by_key(|rw| rw.cost_delta);
        rewrites
    }

    fn apply_rewrite<H: GraphLike>(&self, rewrite: Self::Rewrite, graph: &H) -> RewriteResult<H> {
        let mut graph = graph.clone();
        let cost_delta = rewrite.apply(&mut graph);
        RewriteResult { graph, cost_delta }
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
        for rw_set in rules {
            let rhs_idx = RhsIdx(all_rhs.len());
            all_rhs.push(rw_set.rhss().to_owned());
            let boundary = rw_set.lhs().boundary().collect_vec();
            for (inputs, outputs) in rw_set.lhs().ios() {
                let p = rw_set.lhs().graph();
                let inputs = HashSet::from_iter(inputs);
                let outputs = HashSet::from_iter(outputs);
                patterns.push(CausalPattern::new(p, boundary.clone(), inputs, outputs));
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

#[cfg(test)]
pub(crate) mod test {
    use crate::cost::{CostMetric, TwoQubitGateCount};
    use crate::rewrite_sets::test::rewrite_set_2qb_lc;

    use super::*;
    use quizx::json::decode_graph;
    use quizx::vec_graph::Graph;
    use rstest::{fixture, rstest};

    /// Makes a simple graph.
    ///
    /// The graph is:
    /// ```text
    /// 0 - 8 - 2 - 3 - 6
    ///        / \
    /// 1 --- 4 - 5 - 7
    /// ```
    ///
    /// with inputs 0, 1 and outputs 6, 7.
    #[fixture]
    pub(crate) fn small_graph() -> Graph {
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
            g.add_vertex(VType::Z),
        ];

        g.set_inputs(vec![vs[0], vs[1]]);
        g.set_outputs(vec![vs[6], vs[7]]);

        g.add_edge_with_type(vs[0], vs[8], EType::N);
        g.add_edge_with_type(vs[1], vs[4], EType::N);

        g.add_edge_with_type(vs[8], vs[2], EType::H);
        g.add_edge_with_type(vs[2], vs[3], EType::H);
        g.add_edge_with_type(vs[2], vs[4], EType::H);
        g.add_edge_with_type(vs[2], vs[5], EType::H);
        g.add_edge_with_type(vs[4], vs[5], EType::H);

        g.add_edge_with_type(vs[3], vs[6], EType::N);
        g.add_edge_with_type(vs[5], vs[7], EType::N);

        g
    }

    #[fixture]
    pub(crate) fn json_simple_graph() -> Graph {
        const SIMPLE_GRAPH_JSON: &str = include_str!("../../test_files/simple-graph.json");
        decode_graph(SIMPLE_GRAPH_JSON).unwrap()
    }

    #[fixture]
    pub(crate) fn compiled_rewriter() -> CausalRewriter<Graph> {
        let rw_set = rewrite_set_2qb_lc();
        CausalRewriter::from_rewrite_rules(rw_set)
    }

    #[fixture]
    pub(crate) fn pre_compiled_rewriter() -> CausalRewriter<Graph> {
        const REWRITE_2QB_LC: &[u8] = include_bytes!("../../test_files/rewrites-2qb-lc.rwr");
        rmp_serde::from_slice(REWRITE_2QB_LC).unwrap()
    }

    #[rstest]
    #[case::small_compiled(small_graph(), compiled_rewriter())]
    #[case::json_compiled(json_simple_graph(), compiled_rewriter())]
    #[case::small_precompiled(small_graph(), pre_compiled_rewriter())]
    #[case::json_precompiled(json_simple_graph(), pre_compiled_rewriter())]
    fn test_match_apply(
        #[case] graph: Graph,
        #[case] rewriter: CausalRewriter<Graph>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let cost_metric = TwoQubitGateCount::new();
        let graph_cost = cost_metric.cost(&graph);

        let rewrites = rewriter.get_rewrites(&graph);

        println!("Orig cost {graph_cost}");
        for rw in rewrites {
            let r = rewriter.apply_rewrite(rw, &graph);
            let new_cost = cost_metric.cost(&r.graph);

            println!("New cost {new_cost}");
            assert_eq!(graph_cost.saturating_add_signed(r.cost_delta), new_cost);
        }

        Ok(())
    }
}
