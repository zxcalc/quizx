use num::Zero;
use quizx::params::Expr;
use rustc_hash::FxHashSet;
use std::time::Instant;

use quizx::circuit::*;
use quizx::simplify::interior_clifford_simp;
use quizx::{graph::*, vec_graph::Graph};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let distance = 50;
    let rounds = 50;
    let time = Instant::now();
    println!(
        "Building surface code circuit with distance {} and {} rounds...",
        distance, rounds
    );
    let c = Circuit::surface_code()
        .distance(distance)
        .rounds(rounds)
        .build();
    println!("Done in {:?}", time.elapsed());

    let time = Instant::now();
    println!("Converting to graph...");
    let mut g: Graph = c.to_graph();
    println!("Done in {:?}", time.elapsed());

    println!("Initial ZX diagram has {} vertices", g.num_vertices());

    let time = Instant::now();
    println!("Performing clifford simp...");
    interior_clifford_simp(&mut g);
    println!("Done in {:?}", time.elapsed());

    println!(
        "Result has {} vertices, {} inputs, {} outputs",
        g.num_vertices(),
        g.inputs().len(),
        g.outputs().len()
    );

    let factors: Vec<Expr> = g
        .scalar_factors()
        .filter_map(|(e, s)| if s.is_zero() { Some(e.clone()) } else { None })
        .collect();

    let mut vars = FxHashSet::default();

    let mut size = 0.0;
    for f in factors.iter() {
        size += f[0].len() as f64;
        for v in f[0].iter() {
            vars.insert(v);
        }
    }

    size /= factors.len() as f64;

    assert!(factors.iter().all(|e| e.is_linear()));
    println!(
        "Result has {} linear constraints on {} vars (avg size: {:.2} vars).",
        factors.len(),
        vars.len(),
        size,
    );

    Ok(())
}
