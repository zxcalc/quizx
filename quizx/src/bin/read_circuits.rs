use std::time::Instant;
use quizx::circuit::*;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
