
// use quizx::graph::*;
// use quizx::vec_graph::Graph;
// use quizx::tensor::ToTensor;
// use quizx::basic_rules::*;
// use std::time::Instant;
use ndarray::prelude::*;
use ndarray::{Axis,stack};

fn main() {
    // let a1 = array![[1, 0],[0, 1]];
    let a = array![[1, 2],[3, 4]];
    let b = stack![Axis(0),a,a];
    // let sh = Vec::from(b.shape());
    println!("{}", b);
    println!("{:?}", b.shape());
}
