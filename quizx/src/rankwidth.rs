//! This module uses simulated annealing to find good rank decompositions of a graph
//!
//! This module provides an implementation of a simulated annealing algorithm
//! to find good rank decompositions of a graph. It is based on the master's thesis
//! of Florian Nouwt:
//!
//! "A simulated annealing method for computing rank-width". University of Utrecht, 2022.
//! https://doi.org/20.500.12932/41566
//!
//! Some specifics are ported from Florian's C# implementation, which is available here:
//!
//! https://github.com/Gericom/RankWidthApproximate/
//!
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
