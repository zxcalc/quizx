// use std::time::Instant;
use quizx::circuit::*;
use quizx::vec_graph::*;
use quizx::simplify::*;
use quizx::extract::*;
use quizx::tensor::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Circuit::from_file("../../circuits/mod5_4.qasm")?;
    // let c = c.to_basic_gates();
    // let c1 = c.to_basic_gates();
    // if !Tensor4::scalar_compare(&c, &c1) {
    //     panic!("Tensors don't match: c, c1");
    // }

    // let c = Circuit::random()
    //     .seed(1337)
    //     .qubits(5)
    //     .depth(30)
    //     .p_t(0.2)
    //     .with_cliffords()
    //     .build();
    let mut g: Graph = c.to_graph();
    clifford_simp(&mut g);
    // println!("{}", g.to_dot());
    // println!("{:?}", g);
    // assert_eq!(c.to_tensor4(), g.to_tensor4());

    match g.to_circuit() {
        Ok(c1) => {
            println!("extracted ok");
            if Tensor4::scalar_compare(&c, &c1) {
                println!("Tensors match!");
            } else {
                println!("Tensors don't match. \n{}\n\n{}", c, c1);
            }
        },
        Err(ExtractError(msg, _c, _g)) => {
            println!("extract failed: {}", msg);
            println!("{}\n\n{}\n\n{}", msg, _c, _g.to_dot());
        },
    }
    Ok(())
}
