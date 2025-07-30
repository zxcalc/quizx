//! Uses simulated annealing to find good rank decompositions of a graph
pub mod decomp_tree;

use crate::graph::{GraphLike, V};
use bitgauss::BitMatrix;
use decomp_tree::DecompTree;
use rand::Rng;
use rustc_hash::FxHashMap;

pub struct RankWidthDecomposer<T: GraphLike> {
    graph: T,
    tree: DecompTree,
    ranks: FxHashMap<(usize, usize), usize>,
}

impl<T: GraphLike> RankWidthDecomposer<T> {
    pub fn new(graph: T, tree: DecompTree) -> Self {
        Self {
            graph,
            tree,
            ranks: FxHashMap::default(),
        }
    }

    pub fn compute_ranks(&mut self) {
        for e in self.tree.edges() {
            if !self.ranks.contains_key(&e) {
                let (p0, p1) = self.tree.partition(e);
                let rank = BitMatrix::build(p0.len(), p1.len(), |i, j| {
                    self.graph.connected(p0[i], p1[j])
                })
                .rank();
                self.ranks.insert(e, rank);
            }
        }
    }

    fn add_random_partition(
        parent: usize,
        tree: &mut DecompTree,
        mut vs: Vec<V>,
        rng: &mut impl Rng,
    ) -> usize {
        if vs.len() == 1 {
            tree.add_leaf(parent, vs[0])
        } else if vs.len() >= 2 {
            // split 'vs' into two random, non-empty subsets
            let mut left = vec![vs.remove(rng.random_range(0..vs.len()))];
            let mut right = vec![vs.remove(rng.random_range(0..vs.len()))];
            while !vs.is_empty() {
                if rng.random_bool(0.5) {
                    left.push(vs.remove(rng.random_range(0..vs.len())));
                } else {
                    right.push(vs.remove(rng.random_range(0..vs.len())));
                }
            }

            let i = tree.add_interior([parent, 0, 0]);
            let l = Self::add_random_partition(i, tree, left, rng);
            let r = Self::add_random_partition(i, tree, right, rng);
            tree.nodes_mut()[i].nhd_mut()[1] = l;
            tree.nodes_mut()[i].nhd_mut()[2] = r;

            i
        } else {
            panic!("Attempted to decompose an empty list of vertices");
        }
    }

    pub fn random_decomp(graph: T, rng: &mut impl Rng) -> Self {
        let mut tree = DecompTree::new();

        let mut vs: Vec<V> = graph.vertices().collect();
        if vs.len() >= 2 {
            let v = vs.pop().unwrap();
            tree.add_leaf(0, v);
            let i = Self::add_random_partition(0, &mut tree, vs, rng);
            tree.nodes_mut()[0].nhd_mut()[0] = i;
        }

        RankWidthDecomposer::new(graph, tree)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::VType;
    use crate::vec_graph::Graph;

    use rand::{rngs::SmallRng, SeedableRng};
    use rustc_hash::FxHashSet;

    #[test]
    fn test_random_decomp() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut graph = Graph::new();
        for _ in 0..10 {
            graph.add_vertex(VType::Z);
        }

        let vertex_set: FxHashSet<V> = graph.vertices().collect();

        let decomposer = RankWidthDecomposer::random_decomp(graph, &mut rng);

        for e in decomposer.tree.edges() {
            let (p1, p2) = decomposer.tree.partition(e);
            let p1_set: FxHashSet<V> = p1.into_iter().collect();
            let p2_set: FxHashSet<V> = p2.into_iter().collect();
            assert!(
                p1_set.is_disjoint(&p2_set),
                "Partitions {} and {} overlap",
                p1_set.len(),
                p2_set.len()
            );
            assert!(
                p1_set.union(&p2_set).eq(&vertex_set),
                "{:?} union {:?} != {:?}",
                p1_set,
                p2_set,
                vertex_set
            );
        }
    }
}
