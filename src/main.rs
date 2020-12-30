// use quizx::graph::*;
use quizx::scalar::*;
use num::rational::Rational;

fn main() {
    let s = Scalar::one_plus_phase(Rational::new(1,4));
    let t = Scalar::one_plus_phase(Rational::new(5,4));
    let st = s.mult_with(&t);
    let sf = s.as_float();
    let tf = t.as_float();
    println!("s: {}", s);
    println!("t: {}", t);
    println!("s*t: {}", &s * &t);
    println!("float(s * t): {}", st.as_float());
    println!("float(s) * float(t): {}", sf * tf);
}
