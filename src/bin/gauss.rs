use quizx::linalg::*;
use rustc_hash::FxHashMap;

fn main() {
    let z = Mat2::zeros(3, 4);
    println!("{}", z);

    let id = Mat2::id(5);
    println!("{}", id);

    let u = Mat2::unit_vector(3, 1);
    println!("{}", u);

    let mut chunks: FxHashMap<&[u32],usize> = FxHashMap::default();

    let v = [1, 0, 1, 1, 0, 1];
    chunks.insert(&v[0..3], 1);

    let vc: Vec<u32> = v[1..3].to_vec();

    println!("{:?}", vc);
}
