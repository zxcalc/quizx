use quizx::graph::*;
use quizx::vec_graph::Graph;
// use quizx::hash_graph::Graph;
use quizx::basic_rules::*;
use std::time::Instant;

fn main() {
    let sz = 1_000_000;
    println!("Building Z-spider chain of size: {}...", sz);
    let time = Instant::now();
    let mut g = Graph::new();
    // let mut g = SparseGraph::new();
    g.add_vertex(VType::B);
    g.add_vertex(VType::B);
    g.add_vertex(VType::B);
    g.add_vertex(VType::Z);
    g.add_edge(0, 3);
    g.add_edge(1, 3);
    g.add_edge(2, 3);

    for i in 4..(sz+4) {
        g.add_vertex(VType::Z);
        g.add_edge(i-1, i);
    }

    println!("Done in {:.2?}", time.elapsed());

    println!("Fusing all spiders...");
    let time = Instant::now();

    for i in 4..(sz+4) {
        // let success = spider_fusion(&mut g, 3, i);
        // assert!(success, "Spider fusion failed for v[3] -> v[{}]!", i);
        spider_fusion_unsafe(&mut g, 3, i);
    }

    println!("Done in {:.2?}", time.elapsed());
}
