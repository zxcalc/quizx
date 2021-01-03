// use quizx::graph::*;
use quizx::scalar::*;
// use approx::{assert_relative_eq, assert_abs_diff_eq};

fn main() {
    let s = Scalar::one();
    let t = Scalar::zero();
    let mut u = &s * t;
    u *= &s;
    u *= s;
    println!("eps: {}", f64::EPSILON);
}
