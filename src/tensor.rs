use crate::scalar::*;
use crate::graph::*;
use num::Complex;
use ndarray::prelude::*;
use ndarray::{Array, Axis, stack};

type ComplexTensor = Array<Complex<f64>,IxDyn>;

pub trait ToTensor: IsGraph {
    fn to_tensor(&self) -> ComplexTensor {
        let a = array![Complex::new(1.0,0.0)].into_dyn();
        a
    }
}
