// use crate::scalar::*;
use crate::graph::*;
use crate::scalar::*;
use num::{Complex,Rational};
use ndarray::prelude::*;
use ndarray::{Array, Axis, stack, ScalarOperand};
use std::collections::VecDeque;
use rustc_hash::FxHashMap;

pub type Tensor<A> = Array<A,IxDyn>;
pub type Matrix<A> = Array<A,Ix2>;

impl Sqrt2 for Complex<f64> {
    fn sqrt2() -> Complex<f64> { Complex::new(f64::sqrt(2.0), 0.0) }
    fn one_over_sqrt2() -> Complex<f64> { Complex::new(1.0 / f64::sqrt(2.0), 0.0) }
}

impl FromPhase for Complex<f64> {
    fn from_phase(p: Rational) -> Complex<f64> {
        let exp = (*p.numer() as f64) / (*p.denom() as f64);
        Complex::new(-1.0, 0.0).powf(exp)
    }
}

fn compute_tensor<G,A>(graph: &G) -> Tensor<A>
    where G: IsGraph + Clone,
          A: Copy + Zero + One + Sqrt2 + FromPhase + ScalarOperand
{
    let mut g = graph.clone();
    g.x_to_z();
    // H-boxes are not implemented yet
    for v in g.vertices() {
        let t = g.vertex_type(v);
        if t != VType::B && t != VType::Z {
            panic!("Vertex type currently unsupported: {:?}", t);
        }
    }

    let mut a = array![A::one()].into_dyn();
    let inp = g.inputs().iter().copied();
    let mid = g.vertices().filter(|&v| g.vertex_type(v) != VType::B);
    let outp = g.outputs().iter().copied();
    let mut vs: Vec<V> = inp.chain(mid.chain(outp)).collect();

    if vs.len() < g.num_vertices() {
        panic!("All boundary vertices must be an input or an output");
    }

    vs.reverse();
    // TODO: pick a good sort order for mid

    let mut indexv: VecDeque<V> = VecDeque::new();
    let mut seenv: FxHashMap<V,usize> = FxHashMap::default();

    let (delta, had): (Matrix<A>,Matrix<A>) = {
        let o = A::one();
        let z = A::zero();
        let minus_one = A::from_phase(Rational::new(1,1));
        let one_over_rt2 = A::one_over_sqrt2();
        (array![[o,z],[z,o]],
         array![[o,o],[o,minus_one]] * one_over_rt2)
    };

    let mut fst = true;

    for v in vs {
        let p = g.phase(v);
        if p == Rational::new(0,1) {
            if fst {
                a = array![A::one(), A::one()].into_dyn();
                fst = false;
            } else {
                a = stack![Axis(0), a, a];
            }
        } else {
            let f = A::from_phase(p);
            if fst {
                a = array![A::one(), f].into_dyn();
                fst = false;
            } else {
                a = stack![Axis(0), a, &a * f];
            }
        }


        indexv.push_front(v);
        let mut deg_v = 0;

        for (w, et) in g.incident_edges(v) {
            match seenv.get_mut(&w) {
                Some (deg_w) => {
                    deg_v += 1;
                    *deg_w += 1;

                    let vi = indexv.iter().position(|x| *x == v).unwrap();
                    let mut wi = indexv.iter().position(|x| *x == w).unwrap();

                    // treat delta or hadamard as a K-tensor, with trivial dimensions except
                    // at the indexes of v and w, and multiply it in
                    let mut shape: Vec<usize> = vec![1; indexv.len()];
                    shape[vi] = 2;
                    shape[wi] = 2;
                    println!("{:?}", shape);

                    let m = if et == EType::N { &delta } else { &had }
                        .clone()
                        .into_shape(shape)
                        .expect("Bad tensor indices");

                    a = &m * &a;

                    // if v and w now have all their edges in the tensor, contract away the
                    // index

                    if g.vertex_type(v) != VType::B && g.degree(v) == deg_v {
                        a.sum_axis(Axis(vi));
                        indexv.remove(vi);
                        if wi > vi { wi -= 1; }
                    }

                    if g.vertex_type(w) != VType::B && g.degree(w) == *deg_w {
                        a.sum_axis(Axis(wi));
                        indexv.remove(wi);
                    }
                },
                None => {}
            }
        }
        seenv.insert(v, deg_v);
    }

    a
}

impl crate::vec_graph::Graph {
    pub fn to_tensor<A>(&self) -> Tensor<A>
    where A: Copy + Zero + One + Sqrt2 + FromPhase + ScalarOperand
    { compute_tensor(self) }
}

impl crate::hash_graph::Graph {
    pub fn to_tensor<A>(&self) -> Tensor<A>
    where A: Copy + Zero + One + Sqrt2 + FromPhase + ScalarOperand
    { compute_tensor(self) }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use crate::graph::*;
    use crate::vec_graph::Graph;
    use crate::scalar4::*;

    #[test]
    fn tensor_1() {
        let mut g = Graph::new();
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_edge(0,1);
        let t: Tensor<Scalar4> = g.to_tensor();
        println!("{}", t);
    }
}
