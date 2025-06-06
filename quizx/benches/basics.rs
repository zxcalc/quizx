use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use quizx::circuit::Circuit;
use quizx::graph::GraphLike;
use quizx::vec_graph::Graph;

fn get_test_files() -> Vec<String> {
    vec![
        // "../circuits/small/tof_3.qasm".to_string(),
        // "../circuits/small/tof_5.qasm".to_string(),
        "../circuits/small/tof_10.qasm".to_string(),
    ]
}

fn benchmark_loading_saving_cloning(c: &mut Criterion) {
    for file in get_test_files() {
        let file_name = file.split('/').next_back().unwrap_or("unknown_file");
        let qasm = std::fs::read_to_string(&file)
            .unwrap_or_else(|_| panic!("Failed to read QASM file: {}", file))
            .replace("\r\n", "\n");

        c.bench_function(&format!("loading_saving_circuit_{}", file_name), |b| {
            b.iter(|| {
                let circuit = Circuit::from_qasm(&qasm).unwrap();
                std::hint::black_box(circuit.to_qasm());
            });
        });

        c.bench_function(&format!("converting_circuit_to_graph_{}", file_name), |b| {
            b.iter_batched_ref(
                || Circuit::from_qasm(&qasm).unwrap(),
                |circuit| {
                    let graph: Graph = circuit.to_graph();
                    assert!(graph.num_vertices() > 0);
                },
                BatchSize::SmallInput,
            );
        });

        c.bench_function(&format!("cloning_circuit_{}", file_name), |b| {
            b.iter_batched_ref(
                || Circuit::from_qasm(&qasm).unwrap(),
                |circuit| {
                    let _clone = circuit.clone();
                },
                BatchSize::SmallInput,
            );
        });
    }
}

criterion_group!(benches, benchmark_loading_saving_cloning);
criterion_main!(benches);
