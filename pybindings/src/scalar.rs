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

//! Python-side scalar definition.

use derive_more::{Add, Mul};
use num::complex::Complex;
use num::rational::Rational64;
use num::{FromPrimitive, One, Zero};
use pyo3::prelude::*;
use quizx::scalar::{FromPhase, ScalarN, Sqrt2};

/// A type for exact and approximate representation of complex
/// numbers.
///
/// The exact representation of a scalar is given as an element of
/// D\[omega\], where D is the ring if dyadic rationals and omega is
/// the 2N-th root of unity, represented by its first N coefficients.
/// Addition for this type is O(N) and multiplication O(N^2). Ring
/// elements are stored as a global power of 2 and a list of integer
/// coefficients. This is effectively a floating point number, but
/// with a shared exponent and different behaviour w.r.t. limited
/// precision (namely it panics if big numbers are added to small
/// ones rather than approximating).
///
/// The float representation of a scalar is given as a 64-bit
/// floating point complex number.
#[pyclass]
#[derive(Debug, Clone, Add, Mul, PartialEq)]
pub struct Scalar {
    /// Rust representation of the scalar.
    ///
    /// Fixed to variable length coefficient lengths, as
    /// pyo3 cannot bind to generic types.
    s: ScalarN,
}

impl From<ScalarN> for Scalar {
    fn from(s: ScalarN) -> Self {
        Self { s }
    }
}

impl From<Scalar> for ScalarN {
    fn from(s: Scalar) -> Self {
        s.s
    }
}

#[pymethods]
impl Scalar {
    /// Create a new complex scalar from a pair of floats.
    #[staticmethod]
    pub fn complex(complex: Complex<f64>) -> Self {
        Self { s: complex.into() }
    }

    /// Create a new real scalar from a float number.
    #[staticmethod]
    pub fn real(real: f64) -> Self {
        Self {
            s: ScalarN::real(real),
        }
    }

    /// Returns a scalar value of 1^{i \pi phase}.
    //
    // TODO: Phase should be a `Rational64` instead of a float.
    //       That requires pyo3 to map python's `fractions.Fraction` to rationals.
    //       See https://github.com/PyO3/pyo3/pull/4148
    #[staticmethod]
    pub fn from_phase(phase: f64) -> Self {
        let phase =
            Rational64::from_f64(phase).unwrap_or_else(|| panic!("Invalid phase value {phase}"));
        Self {
            s: ScalarN::from_phase(phase),
        }
    }

    /// Create a scalar from a list of integer coefficients.
    #[staticmethod]
    pub fn from_int_coeffs(coeffs: Vec<isize>) -> Self {
        Self {
            s: ScalarN::from_int_coeffs(&coeffs),
        }
    }

    /// Returns a scalar value of 1 + 1^{i \pi phase}.
    //
    // TODO: Phase should be a `Rational64` instead of a float.
    // See `from_phase`
    #[staticmethod]
    pub fn one_plus_phase(phase: f64) -> Self {
        Self::one() + Self::from_phase(phase)
    }

    /// Returns a scalar for the p-th power of sqrt(2).
    #[staticmethod]
    pub fn sqrt2_pow(p: i32) -> Self {
        Self {
            s: ScalarN::sqrt2_pow(p),
        }
    }

    /// Returns the float representation of the scalar.
    pub fn complex_value(&self) -> Complex<f64> {
        self.s.complex_value()
    }

    /// Returns the scalar multiplied by the n-th power of sqrt(2).
    pub fn mul_sqrt2_pow(&self, n: i32) -> Self {
        let mut s = self.clone();
        s.s.mul_sqrt2_pow(n);
        s
    }

    /// Returns the scalar multiplied by a phase.
    //
    // TODO: Phase should be a `Rational64` instead of a float.
    // See `from_phase`
    pub fn mul_phase(&mut self, phase: f64) -> Self {
        let mut s = self.clone();
        let phase =
            Rational64::from_f64(phase).unwrap_or_else(|| panic!("Invalid phase value {phase}"));
        s.s.mul_phase(phase);
        s
    }

    /// Returns a zero scalar.
    #[staticmethod]
    pub fn zero() -> Self {
        Self { s: ScalarN::zero() }
    }

    /// Returns a one scalar.
    #[staticmethod]
    pub fn one() -> Self {
        Self { s: ScalarN::one() }
    }

    /// Returns `True` if the scalar is zero.
    pub fn is_zero(&self) -> bool {
        self.s.is_zero()
    }

    /// Returns `True` if the scalar is one.
    pub fn is_one(&self) -> bool {
        self.s.is_one()
    }

    /// Encode the scalar in pyzx-compatible JSON format.
    pub fn to_json(&self) -> String {
        let json_scalar = quizx::json::JsonScalar::from_scalar(&self.s);
        serde_json::to_string(&json_scalar).unwrap()
    }

    /// Decode the scalar from pyzx-compatible JSON format.
    #[staticmethod]
    pub fn from_json(json: &str) -> Self {
        let json_scalar: quizx::json::JsonScalar = serde_json::from_str(json).unwrap();
        Self {
            s: json_scalar.to_scalar().unwrap_or_else(|e| panic!("{}", e)),
        }
    }

    /// Returns the complex conjugate of the scalar.
    pub fn conjugate(&self) -> Self {
        Self { s: self.s.conj() }
    }

    pub fn __repr__(&self) -> String {
        format!("{:?}", self.s)
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.s)
    }

    pub fn __add__(&self, other: &Self) -> Self {
        Self {
            s: &self.s + &other.s,
        }
    }

    pub fn __radd__(&self, other: &Self) -> Self {
        Self {
            s: &self.s + &other.s,
        }
    }

    pub fn __iadd__(&mut self, other: &Self) {
        self.s = &self.s + &other.s;
    }

    pub fn __sub__(&self, other: &Self) -> Self {
        Self {
            s: &self.s + (&other.s * ScalarN::minus_one()),
        }
    }

    pub fn __rsub__(&self, other: &Self) -> Self {
        self.__sub__(other)
    }

    pub fn __isub__(&mut self, other: &Self) {
        self.s = self.__sub__(other).s;
    }

    pub fn __mul__(&self, other: &Self) -> Self {
        Self {
            s: &self.s * &other.s,
        }
    }

    pub fn __rmul__(&self, other: &Self) -> Self {
        Self {
            s: &self.s * &other.s,
        }
    }

    pub fn __imul__(&mut self, other: &Self) {
        self.s *= &other.s;
    }

    pub fn __pos__(&self) -> Self {
        self.clone()
    }

    pub fn __neg__(&self) -> Self {
        Self {
            s: &self.s * ScalarN::minus_one(),
        }
    }

    pub fn __float__(&self) -> f64 {
        self.complex_value().re
    }
}
