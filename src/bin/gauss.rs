use quizx::linalg::*;

fn main() {
    println!("ok");
    let z = MatF2::zero(3, 4);
    println!("{}", z);

    let id = MatF2::id(5);
    println!("{}", id);

    println!("{}", id[0][1]);
}
