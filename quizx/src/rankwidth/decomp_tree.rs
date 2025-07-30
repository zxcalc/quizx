use rustc_hash::FxHashSet;

use crate::graph::V;

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
}

pub struct EdgeIter<'a> {
    tree: &'a DecompTree,
    i: usize,
    j: usize,
}

impl<'a> Iterator for EdgeIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while self.i < self.tree.nodes.len() {
            let n = self.tree.nodes[self.i];
            while self.j < n.nhd().len() {
                let j = n.nhd()[self.j];
                self.j += 1;

                // only return one copy of each edge
                if j > self.i {
                    return Some((self.i, j));
                }
            }
            self.i += 1;
            self.j = 0;
        }
        None
    }
}

impl DecompTree {
    pub fn new() -> Self {
        DecompTree {
            nodes: Vec::new(),
            leaves: Vec::new(),
            interior: Vec::new(),
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

    pub fn edges(&self) -> EdgeIter {
        EdgeIter {
            tree: self,
            i: 0,
            j: 0,
        }
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
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn edge_iterator() {
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

        let edges: Vec<_> = tree.edges().collect();
        let expected_edges = vec![(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (3, 6), (3, 7)];
        assert_eq!(edges, expected_edges);
    }
}
