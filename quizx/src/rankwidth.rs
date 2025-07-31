//! Uses simulated annealing to find good rank decompositions of a graph
pub mod annealer;
pub mod decomp_tree;

use crate::graph::GraphLike;
use annealer::RankwidthAnnealer;
use decomp_tree::DecompTree;

pub fn rank_decomp<G: GraphLike>(graph: &G) -> DecompTree<'_, G> {
    let mut rng = rand::rng();
    let mut annealer = RankwidthAnnealer::from_graph(graph, &mut rng);
    annealer.run()
}
