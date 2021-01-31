// QuiZX - Rust library for quantum circuit rewriting and optimisation
//         using the ZX-calculus
// Copyright (C) 2021 - Aleks Kissinger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// use crate::scalar::*;
use crate::graph::*;
use crate::scalar::*;
use crate::circuit::*;
use num::{Complex,Rational};
use ndarray::prelude::*;
use ndarray::parallel::prelude::*;
use ndarray::{Array, Axis, stack, ScalarOperand, IxDyn, SliceInfo, SliceOrIndex};
use std::collections::VecDeque;
use rustc_hash::FxHashMap;

pub type Tensor<A> = Array<A,IxDyn>;
pub type Matrix<A> = Array<A,Ix2>;

impl Sqrt2 for Complex<f64> {
    fn sqrt2_pow(p: i32) -> Complex<f64> {
        let rt2 = Complex::new(f64::sqrt(2.0), 0.0);
        if p == 1 { rt2 } else { rt2.powi(p) }
    }
}

impl FromPhase for Complex<f64> {
    fn from_phase(p: Rational) -> Complex<f64> {
        let exp = (*p.numer() as f64) / (*p.denom() as f64);
        Complex::new(-1.0, 0.0).powf(exp)
    }
}

/// Wraps all the traits we need to compute tensors from ZX-diagrams.
pub trait TensorElem: Copy + Send + Sync +
    Zero + One + Sqrt2 + FromPhase + FromScalar<ScalarN> +
    ScalarOperand + std::ops::MulAssign + std::fmt::Debug {}
impl<T> TensorElem for T
where T: Copy + Send + Sync +
    Zero + One + Sqrt2 + FromPhase + FromScalar<ScalarN> +
    ScalarOperand + std::ops::MulAssign + std::fmt::Debug {}

/// Trait that implements conversion of graphs to tensors
///
/// This implements a generic method [ToTensor::to_tensor] for any number type that
/// implements [TensorElem], as well as two convenience methods [ToTensor::to_tensor4]
/// and [ToTensor::to_tensorf] for [Scalar4] and floating-point [Complex] numbers,
/// respectively.
pub trait ToTensor {
    fn to_tensor<A: TensorElem>(&self) -> Tensor<A>;

    /// Shorthand for `to_tensor::<Scalar4>()`
    fn to_tensor4(&self) -> Tensor<Scalar4> { self.to_tensor() }

    /// Shorthand for `to_tensor::<Complex<f64>>()`
    fn to_tensorf(&self) -> Tensor<Complex<f64>> { self.to_tensor() }
}

pub trait QubitOps {
    fn ident(q: usize) -> Self;
    fn delta(q: usize) -> Self;
    fn cphase(p: Rational, q: usize) -> Self;
    fn hadamard() -> Self;
    fn delta_at(&mut self, qs: &[usize]);
    fn cphase_at(&mut self, p: Rational, qs: &[usize]);
    fn hadamard_at(&mut self, i: usize);
}

impl<A: TensorElem> QubitOps for Tensor<A> {
    fn ident(q: usize) -> Tensor<A> {
        Tensor::from_shape_fn(vec![2;q*2], |ix| {
            if (0..q).all(|i| ix[i] == ix[q+i]) { A::one() } else { A::zero() }
        })
    }

    fn delta(q: usize) -> Tensor<A> {
        Tensor::from_shape_fn(vec![2;q], |ix| {
            if (0..q).all(|i| ix[i] == 0) || (0..q).all(|i| ix[i] == 1) { A::one() }
            else { A::zero() }
        })
    }

    fn cphase(p: Rational, q: usize) -> Tensor<A> {
        Tensor::from_shape_fn(vec![2;q], |ix| {
            if (0..q).all(|i| ix[i] == 1) { A::from_phase(p) } else { A::one() }
        })
    }

    fn hadamard() -> Tensor<A> {
        let n = A::one_over_sqrt2();
        let minus = A::from_phase(Rational::one());
        array![[n, n], [n, minus * n]].into_dyn()
    }

    fn delta_at(&mut self, qs: &[usize]) {
        let mut shape: Vec<usize> = vec![1; self.ndim()];
        for &q in qs { shape[q] = 2; }
        let del: Tensor<A> = Tensor::delta(qs.len())
            .into_shape(shape).expect("Bad indices for delta_at");
        *self *= &del;
    }

    fn cphase_at(&mut self, p: Rational, qs: &[usize]) {
        let mut shape: Vec<usize> = vec![1; self.ndim()];
        for &q in qs { shape[q] = 2; }
        let cp: Tensor<A> = Tensor::cphase(p, qs.len())
            .into_shape(shape).expect("Bad indices for cphase_at");
        *self *= &cp;
    }

    fn hadamard_at(&mut self, q: usize) {
        let n = A::one_over_sqrt2();
        let minus = A::from_phase(Rational::one()); // -1 = e^(i pi)

        let mut slice0 = vec![SliceOrIndex::from(..); self.ndim()];
        let mut slice1 = slice0.clone();
        slice0[q] = SliceOrIndex::from(0);
        slice1[q] = SliceOrIndex::from(1);

        // split into two non-overlapping pieces, where index q=0 and q=1
        let (mut ma, mut mb) = self.multi_slice_mut((
                SliceInfo::<_,IxDyn>::new(slice0).unwrap().as_ref(),
                SliceInfo::<_,IxDyn>::new(slice1).unwrap().as_ref(),
                ));

        // iterate over the pieces together and apply a hadamard to each of the
        // pairs of elements
        par_azip!((a in &mut ma, b in &mut mb) {
            let a1 = *a;
            *a = n * (*a + *b);
            *b = n * (a1 + minus * *b);
        });
    }
}

impl<G: GraphLike + Clone> ToTensor for G {
    fn to_tensor<A: TensorElem>(&self) -> Tensor<A> {
        let mut g = self.clone();
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

        let mut fst = true;
        let mut num_had = 0;

        for v in vs {
            let p = g.phase(v);
            if fst {
                if p == Rational::new(0,1) {
                    a = array![A::one(), A::one()].into_dyn();
                } else {
                    a = array![A::one(), A::from_phase(p)].into_dyn();
                }
                fst = false;
            } else {
                if p == Rational::new(0,1) {
                    a = stack![Axis(0), a, a];
                } else {
                    let f = A::from_phase(p);
                    a = stack![Axis(0), a, &a * f];
                }
            }


            indexv.push_front(v);
            let mut deg_v = 0;

            for (w, et) in g.incident_edges(v) {
                if let Some(deg_w) = seenv.get_mut(&w) {
                    deg_v += 1;
                    *deg_w += 1;

                    let vi = indexv.iter()
                        .position(|x| *x == v)
                        .expect("v should be in indexv");
                    let mut wi = indexv.iter()
                        .position(|x| *x == w)
                        .expect("w should be in indexv");

                    if et == EType::N {
                        a.delta_at(&[vi, wi]);
                    } else {
                        a.cphase_at(Rational::one(), &[vi, wi]);
                        num_had += 1;
                    }

                    // if v and w now have all their edges in the tensor, contract away the
                    // index

                    if g.vertex_type(v) != VType::B && g.degree(v) == deg_v {
                        // println!("contracting v={}, deg_v={}", v, deg_v);
                        a = a.sum_axis(Axis(vi));
                        indexv.remove(vi);
                        if wi > vi { wi -= 1; }
                    }

                    if g.vertex_type(w) != VType::B && g.degree(w) == *deg_w {
                        // println!("contracting w={}, deg_w={}", w, *deg_w);
                        a = a.sum_axis(Axis(wi));
                        indexv.remove(wi);
                    }
                }
            }
            seenv.insert(v, deg_v);
        }

        let s = A::from_scalar(g.scalar()) * A::sqrt2_pow(-num_had);
        a * s
    }
}

impl ToTensor for Circuit {
    fn to_tensor<A: TensorElem>(&self) -> Tensor<A> {
        let q = self.num_qubits();
        let mut a = Tensor::ident(q);
        a
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    // use crate::graph::*;
    use crate::vec_graph::Graph;

    #[test]
    fn tensor_1() {
        let mut g = Graph::new();
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_edge(0,1);
        let t: Tensor<Scalar4> = g.to_tensor();
        println!("{}", t);
    }

    #[test]
    fn had_at() {
        let mut arr: Tensor<Scalar4> = Tensor::ident(1);
        arr.hadamard_at(0);
        assert_eq!(arr, Tensor::hadamard());
        let mut arr: Tensor<Scalar4> = Tensor::ident(2);
        arr.hadamard_at(0);
        arr.hadamard_at(1);
        arr.hadamard_at(0);
        arr.hadamard_at(1);
        assert_eq!(arr, Tensor::ident(2));
    }
}
