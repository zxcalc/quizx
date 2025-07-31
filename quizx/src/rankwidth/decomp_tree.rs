use bitgauss::BitMatrix;
use rand::Rng;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::graph::{GraphLike, V};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecompNode {
    Leaf([usize; 1], V),
    Interior([usize; 3]),
}

impl DecompNode {
    pub fn is_leaf(&self) -> bool {
        matches!(self, DecompNode::Leaf(_, _))
    }

    pub fn is_interior(&self) -> bool {
        matches!(self, DecompNode::Interior(_))
    }

    pub fn parent(self) -> usize {
        if let DecompNode::Leaf([n], _) = self {
            n
        } else {
            panic!("Called parent on an Interior node");
        }
    }

    pub fn nhd(&self) -> &[usize] {
        match self {
            DecompNode::Leaf(nhd, _) => nhd,
            DecompNode::Interior(nhd) => nhd,
        }
    }

    pub fn nhd_mut(&mut self) -> &mut [usize] {
        match self {
            DecompNode::Leaf(nhd, _) => nhd,
            DecompNode::Interior(nhd) => nhd,
        }
    }

    pub fn other_neighbor(&self, ns: &[usize]) -> usize {
        self.nhd()
            .iter()
            .find(|&&n| !ns.contains(&n))
            .copied()
            .expect("No other neighbor found")
    }

    pub fn replace_neighbor(&mut self, old: usize, new: usize) {
        for n in self.nhd_mut().iter_mut() {
            if *n == old {
                *n = new;
                return;
            }
        }
        panic!("Old neighbor {} not found in node {:?}", old, self);
    }
}

#[derive(Debug, Clone)]
pub struct DecompTree {
    pub nodes: Vec<DecompNode>,
    pub leaves: Vec<usize>,
    pub interior: Vec<usize>,
    ranks: FxHashMap<(usize, usize), usize>,
}

impl DecompTree {
    pub fn new() -> Self {
        DecompTree {
            nodes: Vec::new(),
            leaves: Vec::new(),
            interior: Vec::new(),
            ranks: FxHashMap::default(),
        }
    }

    pub fn add_leaf(&mut self, parent: usize, v: V) -> usize {
        let node = DecompNode::Leaf([parent], v);
        let index = self.nodes.len();
        self.nodes.push(node);
        self.leaves.push(index);
        index
    }

    pub fn add_interior(&mut self, nhd: [usize; 3]) -> usize {
        let node = DecompNode::Interior(nhd);
        let index = self.nodes.len();
        self.nodes.push(node);
        self.interior.push(index);
        index
    }

    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        for (i, node) in self.nodes.iter().enumerate() {
            for &j in node.nhd() {
                if i <= j {
                    edges.push((i, j));
                }
            }
        }
        edges
    }

    pub fn num_edges(&self) -> usize {
        (3 * self.interior.len() + self.leaves.len()) / 2
    }

    /// Find a path from vertex `v1` to vertex `v2` in the decomposition tree using a depth-first search.
    fn dfs(
        &self,
        start: usize,
        avoid: &[usize],
        mut visit: impl FnMut(usize) -> bool,
    ) -> Vec<usize> {
        let mut path = Vec::new();
        let mut seen: FxHashSet<usize> = avoid.iter().copied().collect();

        let mut current = start;
        path.push(current);
        seen.insert(current);

        'dfs: while !visit(current) {
            for &n in self.nodes[current].nhd() {
                if !seen.contains(&n) {
                    current = n;
                    path.push(current);
                    seen.insert(current);
                    continue 'dfs;
                }
            }

            path.pop();
            if !path.is_empty() {
                current = *path.last().unwrap();
            } else {
                break;
            }
        }

        path
    }

    /// Find a path from vertex `n1` to vertex `n2` in the decomposition tree using a depth-first search.
    pub fn path(&self, n1: usize, n2: usize) -> Vec<usize> {
        self.dfs(n1, &[], |n| n == n2)
    }

    /// Return the vertex partition defined by the given edge
    pub fn partition(&self, edge: (usize, usize)) -> (Vec<V>, Vec<V>) {
        let mut p1 = Vec::new();
        self.dfs(edge.0, &[edge.1], |n| {
            if let DecompNode::Leaf(_, v) = self.nodes[n] {
                p1.push(v);
            }
            false
        });

        let mut p2 = Vec::new();
        self.dfs(edge.1, &[edge.0], |n| {
            if let DecompNode::Leaf(_, v) = self.nodes[n] {
                p2.push(v);
            }
            false
        });

        (p1, p2)
    }

    pub fn sort_nhds(&mut self) {
        for node in &mut self.nodes {
            node.nhd_mut().sort();
        }
    }

    pub fn swap_subtrees(&mut self, t1: (usize, usize), t2: (usize, usize)) {
        let (p1, c1) = t1;
        let (p2, c2) = t2;

        // swap subtrees by updating connections from parent to child in both directions
        self.nodes[p1].replace_neighbor(c1, c2);
        self.nodes[p2].replace_neighbor(c2, c1);
        self.nodes[c1].replace_neighbor(p1, p2);
        self.nodes[c2].replace_neighbor(p2, p1);
    }

    pub fn move_subtree(&mut self, path: &[usize]) {
        let a = path[0];
        let a1 = path[1];
        let a2 = path[2];

        let b = path[path.len() - 1];
        let b1 = path[path.len() - 2];

        // the node adjacent to a1 that is not in the path
        let ao = self.nodes[a1].other_neighbor(&[a, a2]);

        // connect ao to a2
        self.nodes[a2].replace_neighbor(a1, ao);
        self.nodes[ao].replace_neighbor(a1, a2);

        // connect b to a1
        self.nodes[b].replace_neighbor(b1, a1);
        self.nodes[a1].replace_neighbor(a2, b);

        // connect a1 to b1
        self.nodes[a1].replace_neighbor(ao, b1);
        self.nodes[b1].replace_neighbor(b, a1);
    }

    pub fn compute_ranks(&mut self, graph: &impl GraphLike) {
        if self.num_edges() == self.ranks.len() {
            return;
        }

        for (i, n) in self.nodes.iter().enumerate() {
            for &j in n.nhd() {
                if i <= j && !self.ranks.contains_key(&(i, j)) {
                    let (p0, p1) = self.partition((i, j));
                    let rank =
                        BitMatrix::build(p0.len(), p1.len(), |i, j| graph.connected(p0[i], p1[j]))
                            .rank();
                    self.ranks.insert((i, j), rank);
                }
            }
        }
    }

    pub fn rankwidth(&mut self, graph: &impl GraphLike) -> usize {
        self.compute_ranks(graph);
        self.ranks.values().max().copied().unwrap_or(0)
    }

    pub fn rankwidth_score(&mut self, graph: &impl GraphLike) -> usize {
        self.compute_ranks(graph);
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

    pub fn swap_random_leaves(&mut self, rng: &mut impl rand::Rng) {
        if self.leaves.len() < 2 {
            return; // Not enough leaves to swap
        }

        let i1 = rng.random_range(0..self.leaves.len());
        let mut i2 = rng.random_range(0..self.leaves.len() - 1);
        if i2 >= i1 {
            i2 += 1;
        }

        let l1 = self.leaves[i1];
        let p1 = self.nodes[l1].parent();
        let l2 = self.leaves[i2];
        let p2 = self.nodes[l2].parent();
        self.swap_subtrees((p1, l1), (p2, l2));

        let path = self.path(l1, l2);
        let mut x = path[0];
        for &y in &path[1..] {
            self.clear_rank((x, y));
            x = y;
        }
    }

    pub fn move_random_subtree(&mut self, rng: &mut impl rand::Rng) {
        // Need to have at least 2 nodes with no common neighbors, which can only
        // happen for well-formed trees with at least 6 nodes
        if self.nodes.len() < 6 {
            return;
        }

        let mut path;
        loop {
            let a = rng.random_range(0..self.nodes.len());
            let b = rng.random_range(0..self.nodes.len());
            path = self.path(a, b);

            if path.len() >= 4 {
                break;
            }
        }

        self.move_subtree(&path);

        let mut x = path[0];
        for &y in &path[1..] {
            self.clear_rank((x, y));
            x = y;
        }
    }

    pub fn random_local_swap(&mut self, rng: &mut impl rand::Rng) {
        // Need to have at least 2 adjacent interior nodes, which can only
        // happen for well-formed trees with at least 6 nodes
        if self.nodes.len() < 6 {
            return;
        }

        let c = self.interior[rng.random_range(0..self.interior.len())];
        let n1 = rng.random_range(0..3);
        let mut n2 = rng.random_range(0..2);
        if n2 >= n1 {
            n2 += 1;
        }
        let mut a = self.nodes[c].nhd()[n1];
        let mut b = self.nodes[c].nhd()[n2];

        // ensure b is always an interior node
        if !self.nodes[b].is_interior() {
            if self.nodes[a].is_interior() {
                (b, a) = (a, b);
            } else {
                b = self.nodes[c].other_neighbor(&[a, b]);
            }
        }

        assert!(
            self.nodes[b].is_interior(),
            "Node b should be interior (malformed tree)"
        );

        let d = self.nodes[b].other_neighbor(&[c]);
        self.swap_subtrees((c, a), (b, d));
        self.clear_rank((b, c));
    }

    /// Check the decomposition tree gives a valid partition of the graph
    pub fn is_valid_for_graph(&self, graph: &impl GraphLike) -> bool {
        let mut vs: Vec<V> = graph.vertices().collect();
        vs.sort();

        for e in self.edges() {
            let (mut p1, p2) = self.partition(e);
            if p1.is_empty() || p2.is_empty() {
                return false;
            }
            p1.extend_from_slice(&p2);
            p1.sort();
            if p1 != vs {
                return false;
            }
        }
        true
    }

    fn add_random_partition(&mut self, parent: usize, mut vs: Vec<V>, rng: &mut impl Rng) -> usize {
        if vs.len() == 1 {
            self.add_leaf(parent, vs[0])
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

            let i = self.add_interior([parent, 0, 0]);
            let l = self.add_random_partition(i, left, rng);
            let r = self.add_random_partition(i, right, rng);
            self.nodes[i].nhd_mut()[1] = l;
            self.nodes[i].nhd_mut()[2] = r;

            i
        } else {
            panic!("Attempted to decompose an empty list of vertices");
        }
    }

    pub fn random_decomp(graph: &impl GraphLike, rng: &mut impl Rng) -> Self {
        let mut tree = Self::new();

        let mut vs: Vec<V> = graph.vertices().collect();
        if vs.len() >= 2 {
            let v = vs.pop().unwrap();
            tree.add_leaf(0, v);
            let i = tree.add_random_partition(0, vs, rng);
            tree.nodes[0].nhd_mut()[0] = i;
        }

        tree
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::VType;
    use crate::vec_graph::Graph;

    use rand::{rngs::SmallRng, SeedableRng};

    #[test]
    fn path_simple_tree() {
        let mut tree = DecompTree::new();

        // Create a simple tree structure:
        //      0
        //    / | \
        //   1  2  3
        tree.add_interior([1, 2, 3]); // node 0
        tree.add_leaf(0, 10); // node 1, vertex 10, connects back to 0
        tree.add_leaf(0, 20); // node 2, vertex 20, connects back to 0
        tree.add_leaf(0, 30); // node 3, vertex 30, connects back to 0

        // Test path from node 1 to node 2
        let path = tree.path(1, 2);
        assert_eq!(path, vec![1, 0, 2]);

        // Test path from node 1 to node 3
        let path = tree.path(1, 3);
        assert_eq!(path, vec![1, 0, 3]);

        // Test path from node 2 to node 3
        let path = tree.path(2, 3);
        assert_eq!(path, vec![2, 0, 3]);
    }

    #[test]
    fn path_deeper_tree() {
        let mut tree = DecompTree::new();

        // Create a deeper tree structure:
        //        0
        //      / | \
        //     1  2  3
        //    /\    /\
        //   4 5   6 7
        tree.add_interior([1, 2, 3]); // node 0
        tree.add_interior([0, 4, 5]); // node 1
        tree.add_leaf(0, 20); // node 2
        tree.add_interior([0, 6, 7]); // node 3
        tree.add_leaf(1, 40); // node 4
        tree.add_leaf(1, 50); // node 5
        tree.add_leaf(3, 60); // node 6
        tree.add_leaf(3, 70); // node 7

        // Test path from deep leaf to another deep leaf
        let path = tree.path(4, 6);
        assert_eq!(path, vec![4, 1, 0, 3, 6]);

        // Test path from leaf to intermediate node
        let path = tree.path(4, 3);
        assert_eq!(path, vec![4, 1, 0, 3]);

        // Test path in the opposite direction
        let path = tree.path(3, 4);
        assert_eq!(path, vec![3, 0, 1, 4]);
    }

    #[test]
    fn path_same_node() {
        let mut tree = DecompTree::new();
        tree.add_leaf(0, 10);

        // Path from a node to itself should just be the node
        let path = tree.path(0, 0);
        assert_eq!(path, vec![0]);
    }

    #[test]
    fn partition() {
        // Create this tree:
        //        0
        //      / | \
        //     1  2  3
        //    /\    /\
        //   4 5   6 7
        let mut tree = DecompTree::new();
        tree.add_interior([1, 2, 3]); // node 0
        tree.add_interior([0, 4, 5]); // node 1
        tree.add_leaf(0, 20); // node 2
        tree.add_interior([0, 6, 7]); // node 3
        tree.add_leaf(1, 40); // node 4
        tree.add_leaf(1, 50); // node 5
        tree.add_leaf(3, 60); // node 6
        tree.add_leaf(3, 70); // node 7
        let edge = (0, 3);
        let (mut p1, mut p2) = tree.partition(edge);
        p1.sort();
        p2.sort();
        assert_eq!(p1, vec![20, 40, 50]);
        assert_eq!(p2, vec![60, 70]);
    }

    #[test]
    fn move_subtree() {
        // example from Florian Nouwt's thesis "A simulated annealing method for computing rank-width", p.30
        let mut tree = DecompTree::new();
        tree.add_leaf(11, 0); // node 0 = v1
        tree.add_leaf(12, 0); // node 1 = v2
        tree.add_leaf(13, 0); // node 2 = v3
        tree.add_leaf(15, 0); // node 3 = v4
        tree.add_leaf(15, 0); // node 4 = v5
        tree.add_leaf(9, 0); // node 5 = v6
        tree.add_leaf(10, 0); // node 6 = v7
        tree.add_leaf(9, 0); // node 7 = v8
        tree.add_leaf(13, 0); // node 8 = v9
        tree.add_interior([5, 7, 10]); // node 9 = a
        tree.add_interior([6, 9, 14]); // node 10 = a'
        tree.add_interior([0, 12, 13]); // node 11 = b
        tree.add_interior([1, 11, 14]); // node 12 = b'
        tree.add_interior([2, 8, 11]); // node 13 = x
        tree.add_interior([10, 12, 15]); // node 14 = y
        tree.add_interior([3, 4, 14]); // node 15 = z

        let mut tree1 = DecompTree::new();
        tree1.add_leaf(11, 0); // node 0 = v1
        tree1.add_leaf(12, 0); // node 1 = v2
        tree1.add_leaf(13, 0); // node 2 = v3
        tree1.add_leaf(15, 0); // node 3 = v4
        tree1.add_leaf(15, 0); // node 4 = v5
        tree1.add_leaf(9, 0); // node 5 = v6
        tree1.add_leaf(14, 0); // node 6 = v7
        tree1.add_leaf(9, 0); // node 7 = v8
        tree1.add_leaf(13, 0); // node 8 = v9
        tree1.add_interior([5, 7, 10]); // node 9 = a
        tree1.add_interior([9, 11, 12]); // node 10 = a'
        tree1.add_interior([0, 10, 13]); // node 11 = b
        tree1.add_interior([1, 10, 14]); // node 12 = b'
        tree1.add_interior([2, 8, 11]); // node 13 = x
        tree1.add_interior([6, 12, 15]); // node 14 = y
        tree1.add_interior([3, 4, 14]); // node 15 = z

        let p = tree.path(9, 11);
        tree.move_subtree(&p);
        tree.sort_nhds();
        for i in 0..tree.nodes.len() {
            assert_eq!(
                tree.nodes[i], tree1.nodes[i],
                "Mismatch at node {}:\n{:?} != {:?}",
                i, tree.nodes[i], tree1.nodes[i]
            );
        }
    }

    #[test]
    fn edge_vec() {
        // Create this tree:
        //        0
        //      / | \
        //     1  2  3
        //    /\    /\
        //   4 5   6 7
        let mut tree = DecompTree::new();
        tree.add_interior([1, 2, 3]); // node 0
        tree.add_interior([0, 4, 5]); // node 1
        tree.add_leaf(0, 20); // node 2
        tree.add_interior([0, 6, 7]); // node 3
        tree.add_leaf(1, 40); // node 4
        tree.add_leaf(1, 50); // node 5
        tree.add_leaf(3, 60); // node 6
        tree.add_leaf(3, 70); // node 7

        let edges: Vec<_> = tree.edges();
        let expected_edges = vec![(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (3, 6), (3, 7)];
        assert_eq!(edges, expected_edges);
    }

    #[test]
    fn test_random_decomp() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut graph = Graph::new();
        for _ in 0..10 {
            graph.add_vertex(VType::Z);
        }

        let mut vs: Vec<V> = graph.vertices().collect();
        vs.sort();

        let tree = DecompTree::random_decomp(&graph, &mut rng);

        for e in tree.edges() {
            let (mut p1, p2) = tree.partition(e);
            assert!(
                !p1.is_empty() && !p2.is_empty(),
                "Partition must be non-empty"
            );
            p1.extend_from_slice(&p2);
            p1.sort();
            assert_eq!(p1, vs, "Partition does not match original vertex set");
        }
    }
}
