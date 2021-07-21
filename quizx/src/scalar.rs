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

use num::{integer,Integer};
use num::complex::Complex;
use num::rational::Rational;
pub use num::traits::identities::{Zero,One};
use std::fmt;
use std::cmp::min;
use approx::AbsDiffEq;

/// A type for exact and approximate representation of complex
/// numbers.
///
/// The [Exact] representation of a scalar is given as an element of
/// D\[omega\], where D is the ring if dyadic rationals and omega is
/// the 2N-th root of unity, represented by its first N coefficients.
/// Addition for this type is O(N) and multiplication O(N^2). Ring
/// elements are stored as a global power of 2 and a list of integer
/// coefficients. This is effectively a floating point number, but
/// with a shared exponent and different behaviour w.r.t. limited
/// precision (namely it panics if big numbers are added to small
/// ones rather than approximating).
///
/// The type of the coefficient list is given as a type parameter
/// implementing a trait [Coeffs].  This is to allow fixed N (with an
/// array) or variable N (with a [Vec]).  Only the former is allowed
/// to implement the [Copy] trait, needed for tensor/matrix elements.
///
/// The [Float] representation of a scalar is given as a 64-bit
/// floating point [Complex] number.
#[derive(Debug,Clone)]
pub enum Scalar<T: Coeffs> {
    Exact(i32, T),
    Float(Complex<f64>),
}

/// Adds the ability to take non-integer types modulo 2. The output
/// should be normalised to be in the range (-1,1].
pub trait Mod2 {
    fn mod2(&self) -> Self;
}

impl Mod2 for Rational {
    fn mod2(&self) -> Rational {
        let mut num = self.numer().rem_euclid(2 * *self.denom());
        if num > *self.denom() { num -= 2 * *self.denom(); }
        Rational::new(num, *self.denom())
    }
}

/// Produce a number from rational root of -1.
pub trait FromPhase {
    fn from_phase(p: Rational) -> Self;
    fn minus_one() -> Self;
}

/// Contains the numbers sqrt(2) and 1/sqrt(2), often used for
/// renormalisation of qubit tensors and matrices.
pub trait Sqrt2: Sized {
    fn sqrt2() -> Self { Self::sqrt2_pow(1) }
    fn one_over_sqrt2() -> Self { Self::sqrt2_pow(-1) }
    fn sqrt2_pow(p: i32) -> Self;
}

/// A list of coefficients. We give this as a parameter to allow
/// either fixed-size lists (e.g. [i32;4]) or dynamic ones (e.g.
/// [Vec]\<i32\>). Only the former can be used in tensors and
/// matrices, because they have to implement Copy (the size must be
/// known at compile time).
pub trait Coeffs: Clone + std::ops::IndexMut<usize,Output=isize> {
    fn len(&self) -> usize;
    fn zero() -> Self;
    fn one() -> Self;
    fn new(sz: usize) -> Option<(Self,usize)>;
}

/// Implement Copy whenever our coefficient list allows us to.
impl<T: Coeffs + Copy> Copy for Scalar<T> {}

use Scalar::{Exact,Float};

/// Allows transformation from a scalar.
///
/// We do not use the standard library's [From] trait to avoid a clash
/// when converting Scalar\<S\> to Scalar\<T\>, which is already
/// implemented as a noop for [From] when S = T.
pub trait FromScalar<T> {
    fn from_scalar(s: &T) -> Self;
}

fn lcm_with_padding(n1: usize, n2: usize) -> (usize,usize,usize) {
    if n1 == n2 {
        (n1, 1, 1)
    } else {
        let lcm0 = integer::lcm(n1, n2);
        (lcm0, lcm0 / n1, lcm0 / n2)
    }
}

impl<T: Coeffs> Scalar<T> {
    pub fn complex(re: f64, im: f64) -> Scalar<T> {
        Float(Complex::new(re, im))
    }

    pub fn real(re: f64) -> Scalar<T> {
        Float(Complex::new(re, 0.0))
    }

    pub fn float_value(&self) -> Complex<f64> {
        match self {
            Exact(pow, coeffs) => {
                let omega = Complex::new(-1f64, 0f64).powf(1f64 / (coeffs.len() as f64));
                let pow2 = 2f64.powi(*pow);

                let mut num = Complex::new(0f64, 0f64);
                for i in 0..coeffs.len() {
                    num += pow2 * (coeffs[i] as f64) * omega.powu(i as u32);
                }
                num
            },
            Float(c) => *c
        }
    }

    pub fn mul_sqrt2_pow(&mut self, p: i32) {
        *self *= Scalar::sqrt2_pow(p);
    }

    pub fn mul_phase(&mut self, phase: Rational) {
        *self *= Scalar::from_phase(phase);
    }

    pub fn to_float(&self) -> Scalar<T> {
        Float(self.float_value())
    }

    pub fn one_plus_phase(p: Rational) -> Scalar<T> {
        Scalar::one() + Scalar::from_phase(p)
    }

    pub fn from_int_coeffs(coeffs: &[isize]) -> Scalar<T> {
        match T::new(coeffs.len()) {
            Some((mut coeffs1, pad)) => {
                for i in 0..coeffs.len() {
                    coeffs1[i*pad] = coeffs[i];
                }
                Exact(0, coeffs1).reduce()
            },
            None => panic!("Wrong number of coefficients for scalar type")
        }
    }

    /// Compute the reduced form of the scalar value
    ///
    /// For non-zero scalars, increment the power of 2 as long as the last bit in
    /// every coefficient is 0. For the zero scalar, set the power of 2 to 0.
    fn reduce(mut self) -> Scalar<T> {
        if let Exact(pow, coeffs) = &mut self {
            let mut all_zero = true;
            for i in 0..coeffs.len() {
                if coeffs[i] != 0 { all_zero = false; break; }
            }

            if all_zero {
                *pow = 0;
            } else {
                let one: isize = 1;
                'outer: loop {
                    for i in 0..coeffs.len() {
                        if one & coeffs[i] == one { break 'outer; }
                    }

                    for i in 0..coeffs.len() { coeffs[i] = coeffs[i] >> 1; }
                    *pow += 1;
                }
            }
        }

        self
    }


    /// Compute the complex conjugate of a scalar and return it
    pub fn conj(&self) -> Scalar<T> {
        match self {
            Exact(pow, coeffs) => {
                // create a new coeff list. n.b. this should always be a good size for T, so we unwrap()
                let mut new_coeffs = T::new(coeffs.len()).unwrap().0;

                // copy the real coeff
                new_coeffs[0] = coeffs[0];

                // for each complex coeff, invert the index mod N and add the negative coeff
                // to that position
                for i in 1..coeffs.len() {
                    new_coeffs[coeffs.len() - i] = -coeffs[i];
                }

                Exact(*pow, new_coeffs)
            },
            Float(c) => { Float(c.conj()) },
        }
    }

    // TODO: this is non-trivial (code below is wrong). Think about this some more.
    // /// Returns true if scalar is real
    // pub fn is_real(&self) -> bool {
    //     match self {
    //         Exact(_, coeffs) => {
    //             for i in 0..coeffs.len() { if coeffs[i] != 0 { return false; }}
    //             true
    //         },
    //         Float(c) => c.im == 0.0,
    //     }
    // }

    // /// Returns true if scalar is real and >= 0
    // pub fn is_non_negative(&self) -> bool {
    //     if !self.is_real() { return false; }
    //     match self {
    //         Exact(_, coeffs) => coeffs[0] >= 0,
    //         Float(c) => c.re >= 0.0,
    //     }
    // }
}

impl<T: Coeffs> Zero for Scalar<T> {
    fn zero() -> Scalar<T> {
        Exact(0, T::zero())
    }

    fn is_zero(&self) -> bool {
        *self == Scalar::zero()
    }
}

impl<T: Coeffs> One for Scalar<T> {
    fn one() -> Scalar<T> {
        Exact(0, T::one())
    }

    fn is_one(&self) -> bool {
        *self == Scalar::one()
    }
}

impl<T: Coeffs> Sqrt2 for Scalar<T> {
    fn sqrt2_pow(p: i32) -> Scalar<T> {
        match T::new(4) {
            Some((mut coeffs,pad)) => {
                // we use the fact that when omega = e^(i pi/4), omega - omega^3 = sqrt(2)

                if p % 2 == 0 {
                    // for even p, use: sqrt(2)^p = 2^(p/2)
                    coeffs[0] = 1;
                    Exact(p/2, coeffs)
                } else {
                    // for odd p, use:
                    // sqrt(2)^p = sqrt(2)^(p-1) * sqrt(2) = 2^((p-1)/2) * (omega - omega^3)
                    coeffs[pad] = 1;
                    coeffs[3*pad] = -1;
                    Exact((p-1)/2, coeffs)
                }
            }
            None => Float(Complex::new(2.0f64.powi(p), 0.0f64))
        }
    }
}

impl<T: Coeffs> FromPhase for Scalar<T> {
    fn from_phase(p: Rational) -> Scalar<T> {
        let mut rnumer = *p.numer();
        let mut rdenom = *p.denom();
        match T::new(rdenom as usize) {
            Some((mut coeffs,pad)) => {
                rnumer *= pad as isize;
                rdenom *= pad as isize;
                rnumer = rnumer.rem_euclid(2 * rdenom);
                let sgn = if rnumer >= rdenom {
                    rnumer = rnumer - rdenom;
                    -1
                } else {
                    1
                };
                coeffs[rnumer as usize] = sgn;
                Exact(0, coeffs)
            },
            None => {
                let f = (*p.numer() as f64) / (*p.denom() as f64);
                Float(Complex::new(-1.0,0.0).powf(f))
            }
        }
    }

    fn minus_one() -> Scalar<T> {
        Scalar::from_phase(Rational::one())
    }
}


impl<T: Coeffs> fmt::Display for Scalar<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Exact(pow, coeffs) => {
                // special output for real clifford+T
                if coeffs.len() == 4 && coeffs[1] == -coeffs[3] && coeffs[2] == 0 {
                    if *pow != 0 { write!(f, "2^{} * (", pow)?; }
                    write!(f, "{}", coeffs[0])?;
                    if coeffs[1] != 0 { write!(f, " + {} * sqrt2", coeffs[1])?; }
                    if *pow != 0 { write!(f, ")")?; }
                    return Ok(());
                }

                // otherwise normal exact output
                let mut fst = true;
                for i in 0..coeffs.len() {
                    if !coeffs[i].is_zero() {
                        if fst {
                            fst = false;
                            if *pow != 0 { write!(f, "2^{} * (", pow)?; }
                        } else {
                            write!(f, " + ")?;
                        }

                        write!(f, "{}", coeffs[i])?;
                        // if *pow != 0 { write!(f, " * 2^{}", pow)?; }
                        if i != 0 { write!(f, " * om^{}", i)?; }
                    }
                }

                if fst { write!(f, "0") }
                else {
                    if *pow != 0 { write!(f, ")")?; }
                    Ok(())
                }
            },
            Float(c) => write!(f, "{}", c),
        }
    }
}

// The main implementation of the Mul trait uses references, so
// we don't need to make a copy of the scalars to multiply them.
impl<'a, 'b, T: Coeffs> std::ops::Mul<&'b Scalar<T>> for &'a Scalar<T> {
    type Output = Scalar<T>;

    fn mul(self, rhs: &Scalar<T>) -> Self::Output {
        match (self,rhs) {
            (Float(c), x) => Float(c * x.float_value()),
            (x, Float(c)) => Float(x.float_value() * c),
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                let (lcm, pad0, pad1) = lcm_with_padding(coeffs0.len(), coeffs1.len());
                match T::new(lcm) {
                    Some((mut coeffs,pad)) => {
                        for i in 0..coeffs0.len() {
                            for j in 0..coeffs1.len() {
                                let pos = (i*pad*pad0 + j*pad*pad1).rem_euclid(2*lcm);
                                if pos < lcm {
                                    coeffs[pos] += coeffs0[i] * coeffs1[j];
                                } else {
                                    coeffs[pos - lcm] += -1 * coeffs0[i] * coeffs1[j];
                                }
                            }
                        }

                        Exact(pow0 + pow1, coeffs).reduce()
                    },
                    None => {
                        Float(self.float_value() * rhs.float_value())
                    }
                }
            },
        }
    }
}

// These 3 variations take ownership of one or both args
impl<T: Coeffs> std::ops::Mul<Scalar<T>> for Scalar<T> {
    type Output = Scalar<T>;
    fn mul(self, rhs: Scalar<T>) -> Self::Output { &self * &rhs } }
impl<'a, T: Coeffs> std::ops::Mul<Scalar<T>> for &'a Scalar<T> {
    type Output = Scalar<T>;
    fn mul(self, rhs: Scalar<T>) -> Self::Output { self * &rhs } }
impl<'a, T: Coeffs> std::ops::Mul<&'a Scalar<T>> for Scalar<T> {
    type Output = Scalar<T>;
    fn mul(self, rhs: &Scalar<T>) -> Self::Output { &self * rhs } }

/// Implements *=
impl<'a, T: Coeffs> std::ops::MulAssign<Scalar<T>> for Scalar<T> {
    fn mul_assign(&mut self, rhs: Scalar<T>) {
        *self = &*self * &rhs;
    }
}

// Variation takes ownership of rhs
impl<'a, T: Coeffs> std::ops::MulAssign<&'a Scalar<T>> for Scalar<T> {
    fn mul_assign(&mut self, rhs: &Scalar<T>) { *self = &*self * rhs; } }

// The main implementation of the Add trait uses references, so we
// don't need to make a copy of the scalars to add them.
impl<'a, 'b, T: Coeffs> std::ops::Add<&'b Scalar<T>> for &'a Scalar<T> {
    type Output = Scalar<T>;

    fn add(self, rhs: &Scalar<T>) -> Self::Output {
        // catch zeros early to prevent overflows for very large numbers
        if rhs.is_zero() { return self.clone(); }
        if self.is_zero() { return rhs.clone(); }
        match (self,rhs) {
            (Float(c), x) => Float(c + x.float_value()),
            (x, Float(c)) => Float(x.float_value() + c),
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                let (lcm, pad0, pad1) = lcm_with_padding(coeffs0.len(), coeffs1.len());

                // nb. this can overflow if powers are too different
                let minpow = min(*pow0, *pow1);
                let base0 = 2isize.pow((*pow0 - minpow) as u32);
                let base1 = 2isize.pow((*pow1 - minpow) as u32);

                match T::new(lcm) {
                    Some((mut coeffs, pad)) => {
                        for i in 0..coeffs0.len() {
                            coeffs[i*pad*pad0] += coeffs0[i] * base0;
                        }

                        for i in 0..coeffs1.len() {
                            coeffs[i*pad*pad1] += coeffs1[i] * base1;
                        }

                        Exact(minpow, coeffs).reduce()
                    },
                    None => Float(self.float_value() + self.float_value())
                }
            },
        }
    }
}

// These 3 variations take ownership of one or both args
impl<T: Coeffs> std::ops::Add<Scalar<T>> for Scalar<T> {
    type Output = Scalar<T>;
    fn add(self, rhs: Scalar<T>) -> Self::Output { &self + &rhs }
}

impl<'a, T: Coeffs> std::ops::Add<Scalar<T>> for &'a Scalar<T> {
    type Output = Scalar<T>;
    fn add(self, rhs: Scalar<T>) -> Self::Output { self + &rhs }
}

impl<'a, T: Coeffs> std::ops::Add<&'a Scalar<T>> for Scalar<T> {
    type Output = Scalar<T>;
    fn add(self, rhs: &Scalar<T>) -> Self::Output { &self + rhs }
}

impl<T: Coeffs> FromScalar<Scalar<T>> for Complex<f64> {
    fn from_scalar(s: &Scalar<T>) -> Complex<f64> {
        s.float_value()
    }
}

impl<S: Coeffs, T: Coeffs> FromScalar<Scalar<T>> for Scalar<S> {
    fn from_scalar(s: &Scalar<T>) -> Scalar<S> {
        match s {
            Exact(pow, coeffs) => {
                match S::new(coeffs.len()) {
                    Some((mut coeffs1, pad)) => {
                        for i in 0..coeffs.len() {
                            coeffs1[i*pad] = coeffs[i];
                        }
                        Exact(*pow, coeffs1)
                    },
                    None => Float(s.float_value()),
                }
            },
            Float(c) => Float(*c)
        }
    }
}


impl<T: Coeffs> AbsDiffEq<Scalar<T>> for Scalar<T> {
    type Epsilon = <f64 as AbsDiffEq>::Epsilon;

    // since this is mainly used for testing, we allow rounding errors much bigger than
    // machine-epsilon
    fn default_epsilon() -> Self::Epsilon {
        1e-6f64
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let c1 = self.float_value();
        let c2 = other.float_value();
        f64::abs_diff_eq(&c1.re, &c2.re, epsilon) &&
        f64::abs_diff_eq(&c1.im, &c2.im, epsilon)
    }
}

impl<T: Coeffs> PartialEq for Scalar<T> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Float(c0), Float(c1)) => c0 == c1,
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                // since scalars are reduced via Scalar::reduce(), equal scalars
                // must have the same power of 2.
                if pow0 != pow1 { return false; }
                let (lcm, pad0, pad1) = lcm_with_padding(coeffs0.len(), coeffs1.len());

                let mut all_eq = true;
                for i in 0..lcm {
                    let c0 = if i % pad0 == 0 { coeffs0[i/pad0] } else { 0 };
                    let c1 = if i % pad1 == 0 { coeffs1[i/pad1] } else { 0 };
                    all_eq = all_eq && c0 == c1;
                }

                all_eq
            },
            _ => false
        }
    }
}

/// Implements Coeffs for an array of fixed size $n, and defines
/// the associated scalar type.
macro_rules! fixed_size_scalar {
    ( $name:ident, $n:expr ) => {
        impl Coeffs for [isize;$n] {
            fn len(&self) -> usize { $n }
            fn zero() -> Self { [0;$n] }
            fn one() -> Self {
                let mut a = [0;$n];
                a[0] = 1;
                a
            }
            fn new(sz: usize) -> Option<(Self,usize)> {
                if $n.is_multiple_of(&sz) {
                    Some(([0;$n], $n/sz))
                } else {
                    None
                }
            }
        }

        pub type $name = Scalar<[isize;$n]>;
        impl ndarray::ScalarOperand for $name { }
    }
}

fixed_size_scalar!(Scalar1, 1);
fixed_size_scalar!(Scalar2, 2);
fixed_size_scalar!(Scalar3, 3);
fixed_size_scalar!(Scalar4, 4);
fixed_size_scalar!(Scalar5, 5);
fixed_size_scalar!(Scalar6, 6);
fixed_size_scalar!(Scalar7, 7);
fixed_size_scalar!(Scalar8, 8);

impl Coeffs for Vec<isize> {
    fn len(&self) -> usize { self.len() }
    fn zero() -> Self { vec![0] }
    fn one() -> Self { vec![1] }
    fn new(sz: usize) -> Option<(Self,usize)> {
        Some((vec![0; sz],1))
    }
}

pub type ScalarN = Scalar<Vec<isize>>;

/// tests {{{
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn approx_mul() {
        let s: Scalar4 = Scalar::real(f64::sqrt(0.3) * f64::sqrt(0.3) - 0.3);
        let t: Scalar4 = Scalar::zero();
        assert_ne!(s, t);
        assert_abs_diff_eq!(s, t);
    }

    #[test]
    fn sqrt_i() {
        let s = Scalar4::from_int_coeffs(&[0, 1, 0, 0]);
        assert_abs_diff_eq!(s.to_float(), Scalar::complex(1.0 / f64::sqrt(2.0), 1.0 / f64::sqrt(2.0)));
    }

    #[test]
    fn mul_same_base() {
        let s = Scalar4::from_int_coeffs(&[1, 2, 3, 4]);
        let t = Scalar4::from_int_coeffs(&[4, 5, 6, 7]);
        let st = &s * &t;
        assert!(match st { Exact(_,_) => true, _ => false });
        assert_abs_diff_eq!(st.to_float(), s.to_float() * t.to_float());
    }

    #[test]
    fn phases() {
        let s: ScalarN = ScalarN::from_phase(Rational::new(4,3)) * ScalarN::from_phase(Rational::new(2,5));
        let t: ScalarN = ScalarN::from_phase(Rational::new(4,3) + Rational::new(2,5));
        assert_abs_diff_eq!(s,t);

        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(0,1)),  Scalar4::one());
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(1,1)),  Scalar4::real(-1.0));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(1,2)),  Scalar4::complex(0.0, 1.0));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(-1,2)), Scalar4::complex(0.0, -1.0));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(1,4)),  Scalar4::from_int_coeffs(&[0,1,0,0]));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(3,4)),  Scalar4::from_int_coeffs(&[0,0,0,1]));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(7,4)),  Scalar4::from_int_coeffs(&[0,0,0,-1]));
    }

    #[test]
    fn additions() {
        let s = ScalarN::from_int_coeffs(&[1,2,3,4]);
        let t = ScalarN::from_int_coeffs(&[2,3,4,5]);
        let st = ScalarN::from_int_coeffs(&[3,5,7,9]);
        assert_eq!(s + t, st);
    }

    #[test]
    fn sqrt2_powers() {
        let s = Scalar4::sqrt2_pow(0);
        assert_eq!(s, Scalar4::one());
        let s = Scalar4::sqrt2_pow(2);
        assert_eq!(s, Scalar4::from_int_coeffs(&[2]));
        let s = Scalar4::sqrt2_pow(1);
        assert_abs_diff_eq!(s, Scalar4::real(f64::sqrt(2f64)));
        let s = Scalar4::sqrt2_pow(-1);
        assert_abs_diff_eq!(s, Scalar4::real(1.0 / f64::sqrt(2f64)));

        for p in -7..7 {
            let s = Scalar4::sqrt2_pow(p);
            assert_abs_diff_eq!(s, Scalar4::real(f64::sqrt(2f64).powi(p)));
        }
    }

    #[test]
    fn one_plus_phases() {
        assert_abs_diff_eq!(ScalarN::one_plus_phase(Rational::new(1,1)), ScalarN::zero());

        let plus = ScalarN::one_plus_phase(Rational::new(1,2));
        let minus = ScalarN::one_plus_phase(Rational::new(-1,2));
        assert_abs_diff_eq!(plus * minus, Scalar::real(2.0));
    }

    #[test]
    fn mul_large_power_2() {
        let p1 = Scalar4::sqrt2_pow(200);
        let p2 = Scalar4::sqrt2_pow(-200);
        // multiplying small, large, and/or very different powers of 2 is ok
        let p3 = &p1 * &p2;
        assert_eq!(p3, Scalar4::one());
    }

    #[test]
    fn add_large_power_2() {
        let p1 = Scalar4::sqrt2_pow(200);
        let p2 = Scalar4::sqrt2_pow(210);
        // adding large or small powers of 2 is ok, as long as they are fairly
        // close
        let p3 = &p1 + &p2;

        let q1 = Scalar4::one();
        let q2 = Scalar4::sqrt2_pow(10);
        let q3 = Scalar4::sqrt2_pow(200) * (&q1 + &q2);

        assert_eq!(p3, q3);
    }

    #[test]
    #[should_panic(expected="attempt to multiply with overflow")]
    fn add_diff_power_2() {
        let p1 = Scalar4::sqrt2_pow(200);
        let p2 = Scalar4::sqrt2_pow(-200);
        // adding very different powers of 2 will panic
        let p3 = &p1 + &p2;
        assert_eq!(p3, Scalar4::one());
    }

    #[test]
    fn conjugates() {
        let ps = vec![
            Scalar4::Exact(0, [1, 1, 0, 0]),
            Scalar4::Exact(0, [1, 2, 0, 5]),
            Scalar4::Exact(10, [1, 1, 0, 0]),
            Scalar4::Exact(-3, [1, 1, 1, 1]),
        ];

        for p in ps {
            let p_conj = p.conj();

            let lhs = p.float_value().conj();
            let rhs = p_conj.float_value();
            assert_abs_diff_eq!(lhs.re, rhs.re, epsilon = 0.00001);
            assert_abs_diff_eq!(lhs.im, rhs.im, epsilon = 0.00001);

            let abs = &p * &p_conj;
            let absf = abs.float_value();
            println!("p = {:?}", p);
            println!("p_conj = {:?}", p_conj);
            println!("abs = {:?}", abs);
            assert_abs_diff_eq!(absf.im, 0.0, epsilon = 0.00001);
            assert!(absf.re > 0.0);
        }
    }
}
// }}}
// vim:foldlevel=0:

