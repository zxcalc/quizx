// ACKNOWLEDGEMENT: This entire module is based on a similar proposal in Borghan's master thesis
// for graph states and was suggested to me by Rodatz

use crate::graph::VType;
use crate::hash_graph::{Graph, GraphLike};
use crate::linalg::Mat2;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Pauli {
    X,
    Y,
    Z,
}

/// Represents a Pauli web in a ZX diagram
#[derive(Debug, Default, Clone)]
pub struct PauliWeb {
    /// Maps edge (from, to) to Pauli operator
    /// Note: from < to to ensure consistent ordering
    pub edge_operators: HashMap<(usize, usize), Pauli>,
}

impl PauliWeb {
    /// Create a new empty PauliWeb
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the Pauli operator for an edge between two nodes
    pub fn set_edge(&mut self, from: usize, to: usize, pauli: Pauli) {
        self.edge_operators
            .insert((from.min(to), from.max(to)), pauli);
    }

    /// Get the Pauli operator for an edge between two nodes
    pub fn get_edge(&self, from: usize, to: usize) -> Option<Pauli> {
        self.edge_operators
            .get(&(from.min(to), from.max(to)))
            .copied()
    }

    /// Get the color to use when drawing an edge
    pub fn get_edge_color(&self, from: usize, to: usize) -> Option<&'static str> {
        self.get_edge(from, to).map(|pauli| match pauli {
            Pauli::X => "green", // Green for X operators
            Pauli::Y => "blue",  // Blue for Y operators
            Pauli::Z => "red",   // Red for Z operators
        })
    }
}

fn ordered_nodes(g: &Graph) -> (Vec<usize>, HashMap<usize, usize>) {
    // Get all vertices and sort them for consistent ordering
    let mut original: Vec<usize> = g.vertices().collect();
    original.sort();

    // First put outputs (nodes that are neither inputs nor outputs in the original graph)
    let outputs: Vec<usize> = original
        .iter()
        .filter(|&&v| !g.inputs().contains(&v) && !g.outputs().contains(&v))
        .cloned()
        .collect();

    // Then add the rest (inputs and outputs) that have type != 0 (B type is 0 in Python)
    let mut vertices = outputs.clone();
    vertices.extend(
        original
            .iter()
            .filter(|&&v| {
                let vtype = g.vertex_type(v);
                vtype != VType::B && !outputs.contains(&v)
            })
            .cloned(),
    );

    // Create index map (matrix index -> original node index)
    let index_map: HashMap<usize, usize> =
        vertices.iter().enumerate().map(|(i, &v)| (i, v)).collect();

    // log::debug!("Ordered vertices: {:?}", vertices);
    // log::debug!("Index map: {:?}", index_map);

    (vertices, index_map)
}

pub fn get_pw(index_map: &HashMap<usize, usize>, v: &Mat2, g: &Graph) -> PauliWeb {
    let n_outs = g.inputs().len() + g.outputs().len();
    let mut red_edges = BTreeSet::new();
    let mut green_edges = BTreeSet::new();
    let mut pw = PauliWeb::new();

    // Process each non-zero entry in the matrix row
    // Assuming v is a row vector (1 x n matrix)
    for col in 0..v.num_cols() {
        if v[(0, col)] == 1 {
            let node = *index_map
                .get(&(col - n_outs))
                .expect("Node index not found in index map.");
            let node_color = g.vertex_type(node);

            // Find all edges connected to this node
            for edge in g.edges() {
                if node == edge.0 || node == edge.1 {
                    if node_color == VType::Z {
                        green_edges.insert(edge);
                    } else if node_color == VType::X {
                        red_edges.insert(edge);
                    } else {
                        unreachable!("Unexpected Node color: {:?}", node_color);
                    }
                }
            }
        }
    }

    // Add edges to PauliWeb
    for e in red_edges {
        pw.set_edge(e.0, e.1, Pauli::Z);
    }
    for e in green_edges {
        pw.set_edge(e.0, e.1, Pauli::X);
    }

    pw
}

fn draw_mat(_name: &str, _mat: &Mat2) {
    // log::debug!("Matrix {} ({}x{}):\n{}", _name, _mat.rows(), _mat.cols(), _mat);
}

/// Returns all detection webs of a quizx graph
/// Will inplace convert the graph to rg form
///
/// TODO: perhaps handle the input/output stuff, currently we break it and just assume thats not a set
/// property
pub fn get_detection_webs(g: &mut Graph) -> Vec<PauliWeb> {
    // First convert to RG form
    g.make_rg();

    // Lets make the whole outputs thing native:
    let mut outputs = Vec::new();
    for v in g.vertices() {
        if g.vertex_type(v) == VType::B {
            outputs.push(v);
        }
    }
    g.set_outputs(outputs);

    // Get number of inputs + outputs
    let outs = g.inputs().len() + g.outputs().len();

    // Get ordered nodes and index map
    let (nodelist, index_map) = ordered_nodes(g);
    // log::debug!("Ordered nodes: {:?}", nodelist);
    // log::debug!("outs: {}", outs);

    // Get adjacency matrix in the specified node order
    let big_n = g.adjacency_matrix(Some(&nodelist));
    draw_mat("N (adjacency)", &big_n);

    // Create I_n (identity matrix of size outs x outs)
    let i_n = Mat2::id(outs);
    draw_mat("I_n", &i_n);

    // Create zero block of size (n - outs) x outs
    let zeroblock = Mat2::zeros(big_n.num_rows() - outs, outs);
    draw_mat("zeroblock", &zeroblock);

    // Stack I_n on top of zeroblock vertically
    let mdl = i_n.vstack(&zeroblock);
    draw_mat("mdl", &mdl);

    // Horizontally concatenate mdl and big_n
    let md = mdl.hstack(&big_n);
    draw_mat("md", &md);

    // Create the no_output matrix that will be stacked below md
    // This is [I_{2*outs} | 0] where I is identity and 0 is zero matrix
    let eye_part = Mat2::id(2 * outs);
    let zero_part = Mat2::zeros(2 * outs, md.num_cols() - 2 * outs);
    let no_output = eye_part.hstack(&zero_part);

    // Vertically stack md and no_output
    let md_no_output = md.vstack(&no_output);
    draw_mat("md_no_output", &md_no_output);

    // Compute nullspace
    let mdnons = md_no_output.nullspace();
    // log::debug!("Number of basis vectors in nullspace: {}", mdnons.len());

    // Convert each basis vector to a PauliWeb
    let mut pws = Vec::with_capacity(mdnons.len());
    for basis in mdnons.into_iter() {
        // log::debug!("Basis vector: {}", basis);
        // Create and store the PauliWeb
        let pw = get_pw(&index_map, &basis, g);
        pws.push(pw);
    }

    pws
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_detection_webs() {
        let mut g = Graph::new();
        let webs = get_detection_webs(&mut g);
        assert_eq!(webs.len(), 0);
    }
    #[test]
    fn test_new_pauliweb() {
        let pw = PauliWeb::new();
        assert!(pw.edge_operators.is_empty());
    }

    #[test]
    fn test_set_and_get_edge() {
        let mut pw = PauliWeb::new();

        // Test setting and getting an edge
        pw.set_edge(1, 2, Pauli::X);
        assert_eq!(pw.get_edge(1, 2), Some(Pauli::X));
        assert_eq!(pw.get_edge(2, 1), Some(Pauli::X)); // Should work in both directions

        // Test updating an edge
        pw.set_edge(1, 2, Pauli::Z);
        assert_eq!(pw.get_edge(1, 2), Some(Pauli::Z));

        // Test non-existent edge
        assert_eq!(pw.get_edge(1, 3), None);
    }

    #[test]
    fn test_get_edge_color() {
        let mut pw = PauliWeb::new();

        // Test colors for each Pauli operator
        pw.set_edge(1, 2, Pauli::X);
        pw.set_edge(2, 3, Pauli::Y);
        pw.set_edge(3, 4, Pauli::Z);

        assert_eq!(pw.get_edge_color(1, 2), Some("green"));
        assert_eq!(pw.get_edge_color(2, 3), Some("blue"));
        assert_eq!(pw.get_edge_color(3, 4), Some("red"));
        assert_eq!(pw.get_edge_color(4, 5), None); // Non-existent edge
    }

    #[test]
    fn test_edge_ordering() {
        let mut pw = PauliWeb::new();

        // Test that edge ordering doesn't matter for get/set
        pw.set_edge(2, 1, Pauli::X);
        assert_eq!(pw.get_edge(1, 2), Some(Pauli::X));

        // Test that updating with different order works
        pw.set_edge(1, 2, Pauli::Z);
        assert_eq!(pw.get_edge(2, 1), Some(Pauli::Z));
    }
}
