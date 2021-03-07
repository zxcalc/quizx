// use std::time::Instant;
use quizx::circuit::*;
use quizx::vec_graph::*;
use quizx::simplify::*;
use quizx::extract::*;
use quizx::tensor::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let c = Circuit::from_file("../../circuits/mod5_4.qasm")?;
    let c = Circuit::from_qasm(r#"
qreg qubits[5];
ccx qubits[0],qubits[3],qubits[4];
cx qubits[3],qubits[4];
ccx qubits[0],qubits[1],qubits[4];
    "#)?;
    let c = c.to_basic_gates();
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
    println!("{}", g.to_dot());
    // return Ok(());
    // println!("{}", c1);

    if !Tensor4::scalar_compare(&g, &c) {
         panic!("Tensors don't match: g, c");
    }

    // return Ok(());

    // assert_eq!(c.to_tensor4(), g.to_tensor4());

    // println!("g={}", g.to_dot());

    match g.to_circuit() {
        Ok(c1) => {
            println!("extracted ok");
            if Tensor4::scalar_compare(&c, &c1) {
                println!("Tensors match!");
            } else {
                // println!("Tensors don't match!");
                println!("Tensors don't match. \nc={}\n\nc1={}", c, c1);
                // println!("g scalar: {}", g.scalar());
            }
        },
        Err(ExtractError(msg, _c, _g)) => {
            println!("extract failed: {}", msg);
            println!("{}\n\n{}\n\n{}", msg, _c, _g.to_dot());
        },
    }

    //assert_eq!(c.to_tensor4(), c1.to_tensor4());

    Ok(())
}
