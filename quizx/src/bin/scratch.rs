use std::time::Instant;
use quizx::circuit::*;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let c = Circuit::from_file("../../circuits/mod5_4.qasm")?;
    // let time = Instant::now();
    // println!("loading circuit");
    // let _c = Circuit::from_file("../../circuits/hwb12.qasm")?;
    // println!("done in {:.2?}", time.elapsed());
    // println!("{:?}", c);
    // println!("{}", c);

    for e in fs::read_dir("../../circuits")? {
        if let Some(f) = e?.path().to_str() {
            let time = Instant::now();
            println!("{}", f);
            Circuit::from_file(f).expect(&format!("circuit failed to parse: {}", f));
            println!("...done in {:.2?}", time.elapsed());
        }
    }

    Ok(())
}
