//! Uses simulated annealing to find good rank decompositions of a graph
pub mod annealer;
pub mod decomp_tree;

use crate::graph::{GraphLike, V};
use bitgauss::BitMatrix;
use decomp_tree::DecompTree;
use rand::Rng;
use rustc_hash::FxHashMap;

pub struct RankWidthDecomposer<'a, T: GraphLike> {
    graph: &'a T,
    tree: DecompTree,
    ranks: FxHashMap<(usize, usize), usize>,
}

impl<T: GraphLike> Clone for RankWidthDecomposer<'_, T> {
    fn clone(&self) -> Self {
        RankWidthDecomposer {
            graph: self.graph,
            tree: self.tree.clone(),
            ranks: self.ranks.clone(),
        }
    }
}

impl<'a, T: GraphLike> RankWidthDecomposer<'a, T> {
    pub fn new(graph: &'a T, tree: DecompTree) -> Self {
        Self {
            graph,
            tree,
            ranks: FxHashMap::default(),
        }
    }

    pub fn compute_ranks(&mut self) {
        if self.tree.num_edges() == self.ranks.len() {
            return;
        }

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

    pub fn rankwidth(&mut self) -> usize {
        self.compute_ranks();
        self.ranks.values().max().copied().unwrap_or(0)
    }

    pub fn rankwidth_score(&mut self) -> usize {
        self.compute_ranks();
        self.ranks.values().map(|r| r * r).sum()
    }

    pub fn set_rank(&mut self, e: (usize, usize), rank: usize) {
        let e = if e.0 < e.1 { e } else { (e.1, e.0) };
        self.ranks.insert(e, rank);
    }

    pub fn clear_rank(&mut self, e: (usize, usize)) {
        let e = if e.0 < e.1 { e } else { (e.1, e.0) };
        self.ranks.remove(&e);
    }

    pub fn rank(&mut self, e: (usize, usize)) -> Option<usize> {
        let e = if e.0 < e.1 { e } else { (e.1, e.0) };
        self.ranks.get(&e).copied()
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
            tree.nodes[i].nhd_mut()[1] = l;
            tree.nodes[i].nhd_mut()[2] = r;

            i
        } else {
            panic!("Attempted to decompose an empty list of vertices");
        }
    }

    pub fn swap_random_leaves(&mut self, rng: &mut impl rand::Rng) {
        if self.tree.leaves.len() < 2 {
            return; // Not enough leaves to swap
        }

        let i1 = rng.random_range(0..self.tree.leaves.len());
        let mut i2 = rng.random_range(0..self.tree.leaves.len() - 1);
        if i2 >= i1 {
            i2 += 1;
        }

        let l1 = self.tree.leaves[i1];
        let p1 = self.tree.nodes[l1].parent();
        let l2 = self.tree.leaves[i2];
        let p2 = self.tree.nodes[l2].parent();
        self.tree.swap_subtrees((p1, l1), (p2, l2));

        let path = self.tree.path(l1, l2);
        let mut x = path[0];
        for &y in &path[1..] {
            self.clear_rank((x, y));
            x = y;
        }
    }

    pub fn move_random_subtree(&mut self, rng: &mut impl rand::Rng) {
        // Need to have at least 2 nodes with no common neighbors, which can only
        // happen for well-formed trees with at least 6 nodes
        if self.tree.nodes.len() < 6 {
            return;
        }

        let mut path;
        loop {
            let a = rng.random_range(0..self.tree.nodes.len());
            let b = rng.random_range(0..self.tree.nodes.len());
            path = self.tree.path(a, b);

            if path.len() >= 4 {
                break;
            }
        }

        self.tree.move_subtree(&path);

        let mut x = path[0];
        for &y in &path[1..] {
            self.clear_rank((x, y));
            x = y;
        }
    }

    pub fn random_local_swap(&mut self, rng: &mut impl rand::Rng) {
        // Need to have at least 2 adjacent interior nodes, which can only
        // happen for well-formed trees with at least 6 nodes
        if self.tree.nodes.len() < 6 {
            return;
        }

        let c = self.tree.interior[rng.random_range(0..self.tree.interior.len())];
        let n1 = rng.random_range(0..3);
        let mut n2 = rng.random_range(0..2);
        if n2 >= n1 {
            n2 += 1;
        }
        let mut a = self.tree.nodes[c].nhd()[n1];
        let mut b = self.tree.nodes[c].nhd()[n2];

        // ensure b is always an interior node
        if !self.tree.nodes[b].is_interior() {
            if self.tree.nodes[a].is_interior() {
                (b, a) = (a, b);
            } else {
                b = self.tree.nodes[c].other_neighbor(&[a, b]);
            }
        }

        assert!(
            self.tree.nodes[b].is_interior(),
            "Node b should be interior (malformed tree)"
        );

        let d = self.tree.nodes[b].other_neighbor(&[c]);
        self.tree.swap_subtrees((c, a), (b, d));
        self.clear_rank((b, c));
    }

    pub fn random_decomp(graph: &'a T, rng: &mut impl Rng) -> Self {
        let mut tree = DecompTree::new();

        let mut vs: Vec<V> = graph.vertices().collect();
        if vs.len() >= 2 {
            let v = vs.pop().unwrap();
            tree.add_leaf(0, v);
            let i = Self::add_random_partition(0, &mut tree, vs, rng);
            tree.nodes[0].nhd_mut()[0] = i;
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

    #[test]
    fn test_random_decomp() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut graph = Graph::new();
        for _ in 0..10 {
            graph.add_vertex(VType::Z);
        }

        let mut vs: Vec<V> = graph.vertices().collect();
        vs.sort();

        let decomposer = RankWidthDecomposer::random_decomp(&graph, &mut rng);

        for e in decomposer.tree.edges() {
            let (mut p1, p2) = decomposer.tree.partition(e);
            assert!(p1.len() > 0 && p2.len() > 0, "Partition must be non-empty");
            p1.extend_from_slice(&p2);
            p1.sort();
            assert_eq!(p1, vs, "Partition does not match original vertex set");
        }
    }
}
