// use std::time::Instant;
use quizx::circuit::*;
use std::fs::File;
use std::io::Read;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let mut c = Circuit::new(4);
    // c.add_gate("cz", vec![0,1]);
    // c.add_gate_with_phase("rz", vec![2], Rational::new(1,3));
    // println!("{}", c);

    let mut f = File::open("../../circuits/mod5_4.qasm")?;
    let mut source = String::new();
    f.read_to_string(&mut source)?;
    let c = Circuit::from_qasm(&source)?;

    println!("{:?}", c);
    println!("{}", c);

    Ok(())
}
