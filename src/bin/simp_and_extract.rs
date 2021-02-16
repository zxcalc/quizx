// use std::time::Instant;
use quizx::circuit::*;
use quizx::vec_graph::Graph;
use quizx::simplify::clifford_simp;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Circuit::from_file("../../circuits/adder_8.qasm")?;
    let mut g: Graph = c.to_graph();
    // clifford_simp(&mut g);

    Ok(())
}
