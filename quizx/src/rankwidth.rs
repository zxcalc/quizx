//! Uses simulated annealing to find good rank decompositions of a graph

use rustc_hash::FxHashSet;

use crate::graph::V;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecompNode {
    Leaf(usize, V),
    Interior([usize; 3]),
}

impl DecompNode {
    pub fn parent(self) -> usize {
        if let DecompNode::Leaf(n, _) = self {
            n
        } else {
            panic!("Called parent on an Interior node");
        }
    }

    pub fn children(self) -> [usize; 3] {
        if let DecompNode::Interior(nhd) = self {
            nhd
        } else {
            panic!("Called children on a Leaf node");
        }
    }

    pub fn replace_child(&mut self, old: usize, new: usize) {
        if let DecompNode::Interior(ref mut nhd) = self {
            for n in nhd.iter_mut() {
                if *n == old {
                    *n = new;
                    return;
                }
            }
            panic!("Old child not found in Interior node");
        } else {
            panic!("Called replace_child on a Leaf node");
        }
    }

    pub fn set_parent(&mut self, new_parent: usize) {
        if let DecompNode::Leaf(_, v) = self {
            *self = DecompNode::Leaf(new_parent, *v);
        } else {
            panic!("Called set_parent on an Interior node");
        }
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
    pub fn path(&self, v1: usize, v2: usize) -> Vec<usize> {
        let mut path = Vec::new();
        let mut seen = FxHashSet::default();

        let mut current = v1;
        path.push(current);
        seen.insert(current);

        if let DecompNode::Leaf(n, _) = self.nodes[current] {
            current = n;
            path.push(current);
            seen.insert(current);
        }

        'dfs: while current != v2 {
            if let DecompNode::Interior(nhd) = self.nodes[current] {
                for n in nhd {
                    if !seen.contains(&n) {
                        current = n;
                        path.push(current);
                        seen.insert(current);
                        continue 'dfs;
                    }
                }
            }

            path.pop();
            if !path.is_empty() {
                current = *path.last().unwrap();
            } else {
                panic!("No path found from {} to {}", v1, v2);
            }
        }

        path
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
        let p1 = self.nodes[l1].parent();
        let p2 = self.nodes[l2].parent();

        // swap leaves by updating connections in both directions
        self.nodes[l1].set_parent(p2);
        self.nodes[l2].set_parent(p1);
        self.nodes[p1].replace_child(l1, l2);
        self.nodes[p2].replace_child(l2, l1);
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
            .children()
            .iter()
            .find(|&&n| n != a && n != path[2])
            .unwrap();
        let b2 = self.nodes[b1]
            .children()
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
        tree.nodes.push(DecompNode::Leaf(0, 10)); // node 1, vertex 10, connects back to 0
        tree.nodes.push(DecompNode::Leaf(0, 20)); // node 2, vertex 20, connects back to 0
        tree.nodes.push(DecompNode::Leaf(0, 30)); // node 3, vertex 30, connects back to 0
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
        tree.nodes.push(DecompNode::Leaf(0, 20)); // node 2
        tree.nodes.push(DecompNode::Interior([0, 6, 7])); // node 3
        tree.nodes.push(DecompNode::Leaf(1, 40)); // node 4
        tree.nodes.push(DecompNode::Leaf(1, 50)); // node 5
        tree.nodes.push(DecompNode::Leaf(3, 60)); // node 6
        tree.nodes.push(DecompNode::Leaf(3, 70)); // node 7
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
        tree.nodes.push(DecompNode::Leaf(0, 10));
        tree.leaves = vec![10];

        // Path from a node to itself should just be the node
        let path = tree.path(0, 0);
        assert_eq!(path, vec![0]);
    }
}
