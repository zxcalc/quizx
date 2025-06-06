use quizx::graph::*;
use quizx::decompose::*;

fn main() {
    // Create a graph with a single T gate
    let mut g = Graph::new();
    let v = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
    
    // Add some vertices and edges to make it more interesting
    for i in 1..9 {
        g.add_vertex(VType::Z);
    }
    
    // Add some edges
    g.add_edge_with_type(1, 5, EType::H);
    g.add_edge_with_type(2, 7, EType::H);
    g.add_edge_with_type(3, 8, EType::H);
    g.add_edge_with_type(4, 6, EType::H);
    
    // Create a SingleDecomp for the T gate
    let decomp = SingleDecomp(vec![v]);
    
    // Call eff_alpha
    let alpha = eff_alpha(&g, &decomp);
    println!("Final alpha: {}", alpha);
} 