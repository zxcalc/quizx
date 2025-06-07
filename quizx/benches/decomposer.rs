use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use quizx::circuit::Circuit;
use quizx::decompose::Decomposer;
use quizx::hash_graph::Graph as HashGraph;
use quizx::vec_graph::Graph as VecGraph;

fn get_test_files() -> Vec<String> {
    vec!["../circuits/small/barenco_tof_3.qasm".to_string()]
}

fn benchmark_graph_scalar(c: &mut Criterion) {
    for file in get_test_files() {
        let file_name = file.split('/').next_back().unwrap_or("unknown_file");
        let qasm = std::fs::read_to_string(&file)
            .unwrap_or_else(|_| panic!("Failed to read QASM file: {}", file))
            .replace("\r\n", "\n");
        // Path to the graph file
        let circ = Circuit::from_qasm(&qasm).expect("Failed to create circuit from QASM");
        let vec_graph: VecGraph = circ.to_graph();
        let hash_graph: HashGraph = circ.to_graph();

        // Benchmark the scalar evaluation
        c.bench_function(&format!("evaluate_vec_graph_scalar_{}", file_name), |b| {
            b.iter_batched_ref(
                || vec_graph.clone(), // Clone the graph before timing
                |g| {
                    let mut decomposer = Decomposer::new(g);
                    decomposer.decompose();
                    let scalar = decomposer.scalar();
                    std::hint::black_box(scalar); // Prevent optimization
                },
                BatchSize::SmallInput,
            );
        });

        // Benchmark the scalar evaluation
        c.bench_function(&format!("evaluate_hash_graph_scalar_{}", file_name), |b| {
            b.iter_batched_ref(
                || hash_graph.clone(), // Clone the graph before timing
                |g| {
                    let mut decomposer = Decomposer::new(g);
                    decomposer.decompose();
                    let scalar = decomposer.scalar();
                    std::hint::black_box(scalar); // Prevent optimization
                },
                BatchSize::SmallInput,
            );
        });
    }
}

criterion_group!(benches, benchmark_graph_scalar);
criterion_main!(benches);
