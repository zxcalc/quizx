// use std::time::Instant;
use quizx::circuit::*;
use num::Rational;

fn main() {
   let mut c = Circuit::new(4);
   c.add_gate("cz", vec![0,1]);
   c.add_gate_with_phase("rz", vec![2], Rational::new(1,3));
   println!("{}", c);
}
