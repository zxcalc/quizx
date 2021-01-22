// use std::time::Instant;
use quizx::circuit::*;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let c = Circuit::from_file("../../circuits/mod5_4.qasm")?;
    // println!("{:?}", c);
    // println!("{}", c);

    for e in fs::read_dir("../../circuits")? {
        if let Some(f) = e?.path().to_str() {
            println!("{}", f);
            Circuit::from_file(f).expect(&format!("circuit failed to parse: {}", f));
        }
    }

    Ok(())
}
