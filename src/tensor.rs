// use crate::scalar::*;
use crate::graph::*;
use num::{Complex,Rational};
use num::traits::identities::*;
use ndarray::prelude::*;
use ndarray::{Array, Axis, stack};
use std::collections::VecDeque;
use rustc_hash::FxHashMap;

type ComplexTensor = Array<Complex<f64>,IxDyn>;
type ComplexMatrix = Array<Complex<f64>,Ix2>;

pub trait ToTensor: IsGraph {
    fn to_tensor(&self) -> ComplexTensor {
        let g = self.clone();
        // TODO: g must only contain Z and Boundary vertices
        let mut a = array![Complex::one()].into_dyn();
        let vs: Vec<_> = g.vertices().collect();
        // TODO: pick a good sort order for vs

        let mut indexv: VecDeque<V> = VecDeque::new();
        let mut seenv: FxHashMap<V,usize> = FxHashMap::default();

        let (delta, had): (ComplexMatrix,ComplexMatrix) = {
            let o = Complex::one();
            let z = Complex::zero();
            let one_over_rt2 = Complex::new(1.0 / f64::sqrt(2.0), 0.0);
            (array![[o,z],[z,o]], one_over_rt2 * array![[o,o],[o,-o]])
        };

        for v in vs {
            let p = g.phase(v);
            if p == Rational::new(0,1) {
                a = stack![Axis(0), a, a];
            } else {
                let f = Complex::new(-1.0, 0.0).powf((*p.numer() as f64) / (*p.denom() as f64));
                a = stack![Axis(0), a, f * &a];
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

                        let m = if et == EType::N { &delta } else { &had }
                            .view()
                            .into_shape(shape)
                            .expect("Bad tensor indices");

                        a *= &m;

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
}
