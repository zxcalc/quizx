// use std::time::Instant;
use quizx::circuit::*;

fn main() {
   let mut c = Circuit::new(4);
   c.add_gate("cz", vec![0,1]);
   println!("{}", c);
}
