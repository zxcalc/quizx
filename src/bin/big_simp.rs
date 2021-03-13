// use std::time::Instant;
use quizx::circuit::*;
use quizx::gate::{CNOT,CZ};
use quizx::vec_graph::*;
use quizx::simplify::*;
use quizx::extract::*;
use std::time::Instant;
// use quizx::tensor::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Circuit::from_file("../../circuits/hwb8.qasm")?;
    let mut cx = 0;
    for g in c.gates.iter() {
        if g.t == CNOT || g.t == CZ {
            cx += 1;
        }
    }
    println!("had {} CNOT+CZ", cx);
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

    println!("simplifying...");
    let time = Instant::now();
    interior_clifford_simp(&mut g);
    println!("Done in {:.2?}", time.elapsed());

    let time = Instant::now();
    println!("extracting...");
    // println!("{}", g.to_dot());
    // println!("{:?}", g);
    // assert_eq!(c.to_tensor4(), g.to_tensor4());

    match g.to_circuit() {
        Ok(c1) => {
            println!("Done in {:.2?}", time.elapsed());
            println!("extracted ok");
            let mut cx = 0;
            for g in c1.gates {
                if g.t == CNOT || g.t == CZ {
                    cx += 1;
                }
            }
            println!("got {} CNOT+CZ", cx);
        },
        Err(ExtractError(msg, _c, _g)) => {
            println!("extract failed: {}", msg);
            // println!("{}\n\n{}\n\n{}", msg, _c, _g.to_dot());
        },
    }
    Ok(())
}
