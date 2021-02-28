// use std::time::Instant;
use quizx::circuit::*;
use quizx::vec_graph::*;
use quizx::simplify::*;
use quizx::extract::*;
use quizx::tensor::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Circuit::from_file("../../circuits/mod5_4.qasm")?;
    let c1 = c.to_basic_gates();
    if !Tensor4::scalar_compare(&c, &c1) {
        panic!("Tensors don't match: c, c1");
    }

    // let c = Circuit::random()
    //     .seed(1337)
    //     .qubits(5)
    //     .depth(30)
    //     .p_t(0.2)
    //     .with_cliffords()
    //     .build();
    let mut g: Graph = c1.to_graph();
    interior_clifford_simp(&mut g);
    println!("{}", g.to_dot());
    println!("{}", c1);

    // if !Tensor4::scalar_compare(&g, &c1) {
    //     panic!("Tensors don't match: g, c2");
    // }

    // return Ok(());

    // assert_eq!(c.to_tensor4(), g.to_tensor4());

    // println!("g={}", g.to_dot());

    match g.to_circuit() {
        Ok(c1) => {
            println!("extracted ok");
            if Tensor4::scalar_compare(&c, &c1) {
                println!("Tensors match!");
            } else {
                println!("Tensors don't match!");
                // println!("Tensors don't match. \nc={}\n\nc1={}", c, c1);
                // println!("g scalar: {}", g.scalar());
            }
        },
        Err(ExtractError(msg, _, _)) => {
            println!("extract failed: {}", msg);
            //println!("{}\n\n{}", msg, g1.to_dot());
        },
    }

    //assert_eq!(c.to_tensor4(), c1.to_tensor4());

    Ok(())
}
