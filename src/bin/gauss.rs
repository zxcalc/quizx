use quizx::linalg::*;

fn main() {
    println!("ok");
    let z = Mat2::zero(3, 4);
    println!("{}", z);

    let id = Mat2::id(5);
    println!("{}", id);

    println!("{}", id[0][1]);
}
