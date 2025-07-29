//! Uses simulated annealing to find good rank decompositions of a graph

use rustc_hash::FxHashSet;

use crate::graph::V;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecompNode {
    Leaf([usize; 1], V),
    Interior([usize; 3]),
}

impl DecompNode {
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

pub struct DecompTree {
    pub nodes: Vec<DecompNode>,
    pub leaves: Vec<usize>,
}

impl DecompTree {
    pub fn new() -> Self {
        DecompTree {
            nodes: Vec::new(),
            leaves: Vec::new(),
        }
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

    pub fn swap_random_leaves(&mut self, rng: &mut impl rand::Rng) {
        if self.leaves.len() < 2 {
            return; // Not enough leaves to swap
        }

        let i1 = rng.gen_range(0..self.leaves.len());
        let mut i2 = rng.gen_range(0..self.leaves.len() - 1);
        if i2 >= i1 {
            i2 += 1;
        }

        let l1 = self.leaves[i1];
        let l2 = self.leaves[i2];
        let p1 = self.nodes[l1].nhd()[0];
        let p2 = self.nodes[l2].nhd()[0];

        // swap leaves by updating connections in both directions
        self.nodes[l1].replace_neighbor(p1, p2);
        self.nodes[l2].replace_neighbor(p2, p1);
        self.nodes[p1].replace_neighbor(l1, l2);
        self.nodes[p2].replace_neighbor(l2, l1);
    }

    pub fn move_random_subtree(&mut self, rng: &mut impl rand::Rng) {
        // Need to have at least 2 nodes with no common neighbors, which can only
        // happen for well-formed trees with at least 6 nodes
        if self.nodes.len() < 6 {
            return;
        }

        let mut a;
        let mut b;
        let mut path;

        loop {
            a = rng.gen_range(0..self.nodes.len());
            b = rng.gen_range(0..self.nodes.len());
            path = self.path(a, b);

            if path.len() >= 4 {
                break;
            }
        }

        let a1 = path[1];
        let b1 = path[path.len() - 2];

        // Get the remaining nodes that are not in the path
        let a2 = self.nodes[a1]
            .nhd()
            .iter()
            .find(|&&n| n != a && n != path[2])
            .unwrap();
        let b2 = self.nodes[b1]
            .nhd()
            .iter()
            .find(|&&n| n != b && n != path[path.len() - 3])
            .unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_simple_tree() {
        let mut tree = DecompTree::new();

        // Create a simple tree structure:
        //      0
        //    / | \
        //   1  2  3
        tree.nodes.push(DecompNode::Interior([1, 2, 3])); // node 0
        tree.nodes.push(DecompNode::Leaf([0], 10)); // node 1, vertex 10, connects back to 0
        tree.nodes.push(DecompNode::Leaf([0], 20)); // node 2, vertex 20, connects back to 0
        tree.nodes.push(DecompNode::Leaf([0], 30)); // node 3, vertex 30, connects back to 0
        tree.leaves = vec![1, 2, 3];

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
    fn test_path_deeper_tree() {
        let mut tree = DecompTree::new();

        // Create a deeper tree structure:
        //        0
        //      / | \
        //     1  2  3
        //    /\    /\
        //   4 5   6 7
        tree.nodes.push(DecompNode::Interior([1, 2, 3])); // node 0
        tree.nodes.push(DecompNode::Interior([0, 4, 5])); // node 1
        tree.nodes.push(DecompNode::Leaf([0], 20)); // node 2
        tree.nodes.push(DecompNode::Interior([0, 6, 7])); // node 3
        tree.nodes.push(DecompNode::Leaf([1], 40)); // node 4
        tree.nodes.push(DecompNode::Leaf([1], 50)); // node 5
        tree.nodes.push(DecompNode::Leaf([3], 60)); // node 6
        tree.nodes.push(DecompNode::Leaf([3], 70)); // node 7
        tree.leaves = vec![2, 4, 5, 6, 7];

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
    fn test_path_same_node() {
        let mut tree = DecompTree::new();
        tree.nodes.push(DecompNode::Leaf([0], 10));
        tree.leaves = vec![10];

        // Path from a node to itself should just be the node
        let path = tree.path(0, 0);
        assert_eq!(path, vec![0]);
    }
}
