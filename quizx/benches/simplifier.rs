use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use quizx::circuit::Circuit;
use quizx::simplify::{clifford_simp, flow_simp, full_simp, interior_clifford_simp};
use quizx::vec_graph::*;

fn simp_surface_code(c: &mut Criterion) {
    // initial setup
    let distance = 30;
    let rounds = 30;
    let circuit = Circuit::surface_code()
        .distance(distance)
        .rounds(rounds)
        .build();
    let g: Graph = circuit.to_graph();

    // benchmarking code
    let mut group = c.benchmark_group("surface_code");
    group.sample_size(10); // 10 is the minimum, 100 is default

    group.bench_function("surface_code_interior_clifford_simp", |b| {
        b.iter_batched_ref(
            || g.clone(), // clone the graph before timing
            |g1| {
                // timed application of the simplifier
                interior_clifford_simp(g1);
            },
            // use SmallInput for smaller graphs for less overhead
            BatchSize::LargeInput,
        )
    });

    group.bench_function("surface_code_clifford_simp", |b| {
        b.iter_batched_ref(
            || g.clone(), // clone the graph before timing
            |g1| {
                // timed application of the simplifier
                clifford_simp(g1);
            },
            // use SmallInput for smaller graphs for less overhead
            BatchSize::LargeInput,
        )
    });

    group.bench_function("surface_code_full_simp", |b| {
        b.iter_batched_ref(
            || g.clone(), // clone the graph before timing
            |g1| {
                // timed application of the simplifier
                full_simp(g1);
            },
            // use SmallInput for smaller graphs for less overhead
            BatchSize::LargeInput,
        )
    });

    group.bench_function("surface_code_flow_simp", |b| {
        b.iter_batched_ref(
            || g.clone(), // clone the graph before timing
            |g1| {
                // timed application of the simplifier
                flow_simp(g1);
            },
            // use SmallInput for smaller graphs for less overhead
            BatchSize::LargeInput,
        )
    });
}

criterion_group!(benches, simp_surface_code);
criterion_main!(benches);
