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
use crate::circuit::*;
use crate::graph::*;
use crate::phase::Phase;
use crate::scalar::*;
use ndarray::parallel::prelude::*;
use ndarray::prelude::*;
use ndarray::*;
use num::{Complex, Rational64};
use rustc_hash::FxHashMap;
use std::collections::VecDeque;
use std::convert::TryFrom;
use std::iter::FromIterator;

/// Generic tensor type used by quizx
pub type Tensor<A> = Array<A, IxDyn>;

/// Shorthand for tensors over [`Scalar4`]
pub type Tensor4 = Tensor<Scalar4>;

/// Shorthand for tensors over floating point complex numbers
pub type TensorF = Tensor<Complex<f64>>;

impl Sqrt2 for Complex<f64> {
    fn sqrt2_pow(p: i32) -> Complex<f64> {
        let rt2 = Complex::new(f64::sqrt(2.0), 0.0);
        if p == 1 {
            rt2
        } else {
            rt2.powi(p)
        }
    }
}

impl FromPhase for Complex<f64> {
    fn from_phase(p: impl Into<Phase>) -> Complex<f64> {
        let p = p.into().to_rational();
        let exp = (*p.numer() as f64) / (*p.denom() as f64);
        Complex::new(-1.0, 0.0).powf(exp)
    }

    fn minus_one() -> Complex<f64> {
        Self::from_phase(Phase::one())
    }
}

/// Wraps all the traits we need to compute tensors from ZX-diagrams.
pub trait TensorElem:
    Copy
    + Send
    + Sync
    + PartialEq
    + Zero
    + One
    + Sqrt2
    + FromPhase
    + TryFrom<Scalar4, Error: std::fmt::Debug>
    + ScalarOperand
    + std::ops::MulAssign
    + std::fmt::Debug
{
}
impl<T> TensorElem for T where
    T: Copy
        + Send
        + Sync
        + PartialEq
        + Zero
        + One
        + Sqrt2
        + FromPhase
        + TryFrom<Scalar4, Error: std::fmt::Debug>
        + ScalarOperand
        + std::ops::MulAssign
        + std::fmt::Debug
{
}

/// Trait that implements conversion of graphs to tensors
///
/// This implements a generic method [ToTensor::to_tensor] for any number type that
/// implements [TensorElem], as well as two convenience methods [ToTensor::to_tensorf]
/// and [ToTensor::to_tensorf] for [`Scalar4`] and floating-point [`Complex`] numbers,
/// respectively.
pub trait ToTensor {
    fn to_tensor<A: TensorElem>(&self) -> Tensor<A>;

    /// Shorthand for `to_tensor::<Tensor4>()`
    fn to_tensor4(&self) -> Tensor4 {
        self.to_tensor()
    }

    /// Shorthand for `to_tensor::<TensorF>()`
    fn to_tensorf(&self) -> TensorF {
        self.to_tensor()
    }
}

pub trait QubitOps<A: TensorElem> {
    fn ident(q: usize) -> Self;
    fn delta(q: usize) -> Self;
    fn cphase(p: impl Into<Phase>, q: usize) -> Self;
    fn hadamard() -> Self;
    fn delta_at(&mut self, qs: &[usize]);
    fn cphase_at(&mut self, p: impl Into<Phase>, qs: &[usize]);
    fn hadamard_at(&mut self, i: usize);

    /// split into two non-overlapping pieces, where index q=0 and q=1
    fn slice_qubit_mut(&mut self, q: usize) -> (ArrayViewMut<A, IxDyn>, ArrayViewMut<A, IxDyn>);

    /// contract the last n qubit indices with the first n qubits of other
    ///
    /// panics if n is greater than the number of qubits of self or other.
    fn plug_n_qubits(self, n: usize, other: &Tensor<A>) -> Tensor<A>;
}

pub trait CompareTensors {
    fn scalar_eq(t0: &Self, t1: &Self) -> bool;
    fn compare(x0: &impl ToTensor, x1: &impl ToTensor) -> bool;
    fn scalar_compare(x0: &impl ToTensor, x1: &impl ToTensor) -> bool;
}

impl<A: TensorElem> CompareTensors for Tensor<A> {
    fn scalar_eq(t0: &Tensor<A>, t1: &Tensor<A>) -> bool {
        // if dimensions are different, tensors are different
        if t0.dim() != t1.dim() {
            return false;
        }

        // find the first non-zero element of each tensor
        let a0 = t0.iter().find(|s| !s.is_zero());
        let a1 = t1.iter().find(|s| !s.is_zero());
        match (a0, a1) {
            // if both tensors have a non-zero element, there are 2 cases
            (Some(b0), Some(b1)) => {
                // if the non-zero element is equal, tensors should be equal on the nose
                if b0 == b1 {
                    t0 == t1
                }
                // if they are different, we cross-multiply to check scalar equivalence
                else {
                    t0 * *b1 == t1 * *b0
                }
            }
            // all-zero tensors of the same dimension are equal
            (None, None) => true,
            // otherwise, one is a zero tensor and the other is non-zero
            _ => false,
        }
    }

    fn scalar_compare(x0: &impl ToTensor, x1: &impl ToTensor) -> bool {
        let t0 = x0.to_tensor::<A>();
        let t1 = x1.to_tensor::<A>();
        Tensor::scalar_eq(&t0, &t1)
    }

    fn compare(x0: &impl ToTensor, x1: &impl ToTensor) -> bool {
        x0.to_tensor::<A>() == x1.to_tensor::<A>()
    }
}

impl<A: TensorElem> QubitOps<A> for Tensor<A> {
    fn slice_qubit_mut(&mut self, q: usize) -> (ArrayViewMut<A, IxDyn>, ArrayViewMut<A, IxDyn>) {
        let slice0: SliceInfo<_, IxDyn, IxDyn> =
            SliceInfo::try_from(Vec::from_iter((0..self.ndim()).map(|i| {
                if i == q {
                    SliceInfoElem::from(0)
                } else {
                    SliceInfoElem::from(..)
                }
            })))
            .unwrap();

        let slice1: SliceInfo<_, IxDyn, IxDyn> =
            SliceInfo::try_from(Vec::from_iter((0..self.ndim()).map(|i| {
                if i == q {
                    SliceInfoElem::from(1)
                } else {
                    SliceInfoElem::from(..)
                }
            })))
            .unwrap();

        self.multi_slice_mut((slice0.as_ref(), slice1.as_ref()))
    }

    fn ident(q: usize) -> Tensor<A> {
        Tensor::from_shape_fn(vec![2; q * 2], |ix| {
            if (0..q).all(|i| ix[i] == ix[q + i]) {
                A::one()
            } else {
                A::zero()
            }
        })
    }

    fn delta(q: usize) -> Tensor<A> {
        Tensor::from_shape_fn(vec![2; q], |ix| {
            if (0..q).all(|i| ix[i] == 0) || (0..q).all(|i| ix[i] == 1) {
                A::one()
            } else {
                A::zero()
            }
        })
    }

    fn cphase(p: impl Into<Phase>, q: usize) -> Tensor<A> {
        let mut t = Tensor::ident(q);
        let qs: Vec<_> = (0..q).collect();
        t.cphase_at(p, &qs);
        t
    }

    fn hadamard() -> Tensor<A> {
        let n = A::one_over_sqrt2();
        let minus = A::from_phase(1);
        array![[n, n], [n, minus * n]].into_dyn()
    }

    fn delta_at(&mut self, qs: &[usize]) {
        let mut shape: Vec<usize> = vec![1; self.ndim()];
        for &q in qs {
            shape[q] = 2;
        }
        let del: Tensor<A> = Tensor::delta(qs.len())
            .into_shape_with_order(shape)
            .expect("Bad indices for delta_at");
        *self *= &del;
    }

    fn cphase_at(&mut self, p: impl Into<Phase>, qs: &[usize]) {
        let p = p.into();
        let mut shape: Vec<usize> = vec![1; self.ndim()];
        for &q in qs {
            shape[q] = 2;
        }
        let cp: Tensor<A> = Tensor::from_shape_fn(vec![2; qs.len()], |ix| {
            if (0..qs.len()).all(|i| ix[i] == 1) {
                A::from_phase(p)
            } else {
                A::one()
            }
        })
        .into_shape_with_order(shape)
        .expect("Bad indices for cphase_at");
        *self *= &cp;
    }

    fn hadamard_at(&mut self, q: usize) {
        let n = A::one_over_sqrt2();
        let minus = A::from_phase(1); // -1 = e^(i pi)

        // split into two non-overlapping pieces, where index q=0 and q=1
        let (mut ma, mut mb) = self.slice_qubit_mut(q);

        // iterate over the pieces together and apply a hadamard to each of the
        // pairs of elements
        par_azip!((a in &mut ma, b in &mut mb) {
            let a1 = *a;
            *a = n * (*a + *b);
            *b = n * (a1 + minus * *b);
        });
    }

    fn plug_n_qubits(self, n: usize, other: &Tensor<A>) -> Tensor<A> {
        let d1 = self.shape().len();
        let d2 = other.shape().len();
        let shape1: Vec<usize> = (0..(d1 + d2 - n))
            .map(|i| if i < d1 { 2 } else { 1 })
            .collect();
        let shape2: Vec<usize> = (0..(d1 + d2 - n))
            .map(|i| if i < d1 - n { 1 } else { 2 })
            .collect();

        let t1 = self
            .into_shared()
            .into_shape_with_order(shape1)
            .expect("Invalid tensor reshape");
        let t1p = t1.broadcast(vec![2; d1 + n]).unwrap();
        let t2 = other
            .clone()
            .into_shared()
            .into_shape_with_order(shape2)
            .expect("Invalid tensor reshape");
        let mut t3 = &t1p * &t2;
        for _ in 0..n {
            t3 = t3.sum_axis(Axis(d1 - n));
        }

        t3
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
                panic!("Vertex type currently unsupported: {t:?}");
            }
        }

        // initialise the trivial tensor
        let mut a = Tensor::from_shape_vec(vec![], vec![A::one()]).unwrap();
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
        let mut seenv: FxHashMap<V, usize> = FxHashMap::default();

        // let mut num_had = 0;
        // let mut i = 1;
        // let tot = vs.len();

        for v in vs {
            // println!("contracting {} ({}/{})", v, i, tot);
            // i += 1;
            let p = g.phase(v);

            // the stack! call computes the tensor product of a new spider
            // (1, e^(i pi p)) with the existing tensor 'a'
            if p.is_zero() {
                a = stack![Axis(0), a, a];
            } else {
                let f = A::from_phase(p);
                a = stack![Axis(0), a, &a * f];
            }

            indexv.push_front(v);
            let mut deg_v = 0;

            for (w, et) in g.incident_edges(v) {
                if let Some(deg_w) = seenv.get_mut(&w) {
                    deg_v += 1;
                    *deg_w += 1;

                    let wi = indexv
                        .iter()
                        .position(|x| *x == w)
                        .expect("w should be in indexv");

                    if et == EType::N {
                        a.delta_at(&[0, wi]);
                    } else {
                        a.cphase_at(1, &[0, wi]);
                        // TODO incorporate with cphase_at
                        a *= A::one_over_sqrt2();
                        // num_had += 1;
                    }

                    if g.vertex_type(w) != VType::B && g.degree(w) == *deg_w {
                        a = a.sum_axis(Axis(wi));
                        indexv.remove(wi);
                    }
                }
            }

            if g.vertex_type(v) != VType::B && g.degree(v) == deg_v {
                a = a.sum_axis(Axis(0));
                indexv.remove(0);
            }

            seenv.insert(v, deg_v);
        }

        let s = A::try_from(*g.scalar()).unwrap(); // * A::sqrt2_pow(-num_had);
        a * s
    }
}

impl ToTensor for Circuit {
    fn to_tensor<A: TensorElem>(&self) -> Tensor<A> {
        use crate::gate::GType::*;
        let q = self.num_qubits();

        // start with the identity matrix
        let mut a = Tensor::ident(q);

        // since we are applying the gates to the input indices, this actually
        // computes the transpose of the circuit, but all the gates are self-
        // transposed, so we can get the circuit itself if we just reverse the order.
        for g in self.gates.iter().rev() {
            match g.t {
                ZPhase => a.cphase_at(g.phase, &g.qs),
                Z | CZ | CCZ => a.cphase_at(1, &g.qs),
                S => a.cphase_at(Rational64::new(1, 2), &g.qs),
                T => a.cphase_at(Rational64::new(1, 4), &g.qs),
                Sdg => a.cphase_at(Rational64::new(-1, 2), &g.qs),
                Tdg => a.cphase_at(Rational64::new(-1, 4), &g.qs),
                HAD => a.hadamard_at(g.qs[0]),
                NOT => {
                    a.hadamard_at(g.qs[0]);
                    a.cphase_at(Rational64::one(), &g.qs);
                    a.hadamard_at(g.qs[0]);
                }
                XPhase => {
                    a.hadamard_at(g.qs[0]);
                    a.cphase_at(g.phase, &g.qs);
                    a.hadamard_at(g.qs[0]);
                }
                CNOT => {
                    a.hadamard_at(g.qs[1]);
                    a.cphase_at(Rational64::one(), &g.qs);
                    a.hadamard_at(g.qs[1]);
                }
                TOFF => {
                    a.hadamard_at(g.qs[2]);
                    a.cphase_at(Rational64::one(), &g.qs);
                    a.hadamard_at(g.qs[2]);
                }
                SWAP => a.swap_axes(g.qs[0], g.qs[1]),
                // n.b. these are pyzx-specific gates
                XCX => {
                    a.hadamard_at(g.qs[0]);
                    a.hadamard_at(g.qs[1]);
                    a.cphase_at(g.phase, &g.qs);
                    a.hadamard_at(g.qs[0]);
                    a.hadamard_at(g.qs[1]);
                }
                // TODO: these "gates" are not implemented yet
                ParityPhase => {
                    panic!("Unsupported gate: ParityPhase")
                }
                InitAncilla => {
                    panic!("Unsupported gate: InitAncilla")
                }
                PostSelect => {
                    panic!("Unsupported gate: PostSelect")
                }
                Measure => {
                    panic!("Unsupported gate: Measure")
                }
                MeasureReset => {
                    panic!("Unsupported gate: MeasureReset")
                }
                UnknownGate => {} // unknown gates are quietly ignored
            }
        }
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
        g.add_edge(0, 1);
        let t: Tensor<Scalar4> = g.to_tensor();
        println!("{t}");
    }

    #[test]
    fn tensor_id() {
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_edge(0, 1);
        g.set_inputs(vec![0]);
        g.set_outputs(vec![1]);
        let t: Tensor<Scalar4> = g.to_tensor();
        assert_eq!(t, Tensor::ident(1));

        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_edge(0, 2);
        g.add_edge(2, 1);
        g.set_inputs(vec![0]);
        g.set_outputs(vec![1]);
        let t: Tensor<Scalar4> = g.to_tensor();
        assert_eq!(t, Tensor::ident(1));
    }

    #[test]
    fn tensor_delta() {
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_edge(0, 4);
        g.add_edge(1, 5);
        g.add_edge_with_type(4, 5, EType::N);
        g.add_edge(2, 4);
        g.add_edge(3, 5);
        g.set_inputs(vec![0, 1]);
        g.set_outputs(vec![2, 3]);
        let t = g.to_tensor4();
        println!("{t}");
        assert_eq!(t, Tensor::delta(4));
    }

    #[test]
    fn tensor_cz() {
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_edge(0, 2);
        g.add_edge(1, 3);
        g.add_edge_with_type(2, 3, EType::H);
        g.add_edge(2, 4);
        g.add_edge(3, 5);
        g.set_inputs(vec![0, 1]);
        g.set_outputs(vec![4, 5]);
        g.scalar_mut().mul_sqrt2_pow(1);
        let t = g.to_tensor4();
        println!("{t}");
        assert_eq!(t, Tensor::cphase(Rational64::one(), 2));
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

    #[test]
    fn circuit_eqs() {
        let c1 = Circuit::from_qasm(
            r#"
        qreg q[2];
        cx q[0], q[1];
        cx q[1], q[0];
        cx q[0], q[1];
        "#,
        )
        .unwrap();

        let c2 = Circuit::from_qasm(
            r#"
        qreg q[2];
        swap q[0], q[1];
        "#,
        )
        .unwrap();

        println!("{}", c1.to_tensor4());
        println!("{}", c2.to_tensor4());
        assert_eq!(c1.to_tensor4(), c2.to_tensor4());
    }

    #[test]
    fn tensor_plug() {
        let c1 = Circuit::from_qasm(
            r#"
        qreg q[2];
        cz q[0], q[1];
        "#,
        )
        .unwrap();

        let c2 = Circuit::from_qasm(
            r#"
        qreg q[2];
        cx q[0], q[1];
        "#,
        )
        .unwrap();

        let c3 = &c1 + &c2;

        let t1 = c1.to_tensor4();
        let t2 = c2.to_tensor4();
        let t3 = t1.plug_n_qubits(2, &t2);

        println!("{t3}");
        println!("{}", c3.to_tensor4());

        assert_eq!(t3, c3.to_tensor4());
    }
}
