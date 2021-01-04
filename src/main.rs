use quizx::graph::*;
use quizx::vec_graph::Graph;
// use quizx::hash_graph::Graph;
use quizx::basic_rules::*;
use std::time::Instant;

fn main() {
    let sz = 50_000;
    println!("Building Z-spider chain of size: {}...", sz);
    let time = Instant::now();
    let mut g = Graph::new();
    // let mut g = SparseGraph::new();
    g.add_vertex(VType::Z);

    for i in 1..sz {
        g.add_vertex(VType::Z);
        g.add_edge(i-1, i);
    }

    println!("Done in {:.2?}", time.elapsed());

    println!("Fusing all spiders...");
    let time = Instant::now();

    for i in 1..sz {
        // let success = spider_fusion(&mut g, 3, i);
        // assert!(success, "Spider fusion failed for v[3] -> v[{}]!", i);
        spider_fusion(&mut g, 0, i);
    }

    println!("Done in {:.2?}", time.elapsed());
}
