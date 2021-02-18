// use std::time::Instant;
use quizx::circuit::*;
use quizx::vec_graph::*;
use quizx::simplify::*;
use quizx::extract::*;
// use quizx::tensor::ToTensor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Circuit::from_file("../../circuits/adder_8.qasm")?;
    let mut g: Graph = c.to_graph();
    interior_clifford_simp(&mut g);

    match g.to_circuit() {
        Ok(_) => {},
        Err(ExtractError(msg, _, g1)) => {
            println!("{}\n\n{}", msg, g1.to_dot());
        },
    }

    //assert_eq!(c.to_tensor4(), c1.to_tensor4());

    Ok(())
}
