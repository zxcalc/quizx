use quizx::graph::*;
// use quizx::vec_graph::Graph;
use quizx::hash_graph::Graph;
use quizx::basic_rules::*;
use std::time::Instant;

fn main() {
    let sz = 100_000;
    println!("Building Z-spider chain of size: {}...", sz);
    let time = Instant::now();
    let mut g = Graph::new();
    g.add_vertex(VType::Z);

    for i in 1..sz {
        g.add_vertex(VType::Z);
        g.add_edge(i-1, i);
    }

    println!("Done in {:.2?}", time.elapsed());

    assert_eq!(g.num_vertices(), sz);

    println!("Fusing all spiders...");
    let time = Instant::now();

    loop {
        match g.find_edge(|v0,v1,_| check_spider_fusion(&g, v0, v1)) {
            Some((v0,v1,_)) => spider_fusion_unsafe(&mut g, v0, v1),
            None => break,
        };
    }

    println!("Done in {:.2?}", time.elapsed());
    assert_eq!(g.num_vertices(), 1);
}
