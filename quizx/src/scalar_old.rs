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
use approx::AbsDiffEq;

/// A type for exact and approximate representation of complex
/// numbers.
///
/// The [Exact] representation of a scalar is given as a power of
/// sqrt(2) and an element of Z\[omega\], where omega is the 2N-th
/// root of unity, represented by its first N coefficients. Addition
/// for this type is O(N) and multiplication O(N^2).
///
/// The type of the coefficient list is given as a type parameter
/// implementing a trait [Coeffs]. This is to allow fixed N (with
/// an array) or variable N (with a [Vec]).  Only the former is
/// allowed to implement the [Copy] trait, needed for tensor/matrix
/// elements.
///
/// The [Float] representation of a scalar is given as a 64-bit
/// floating point [Complex] number.
///
/// TODO: Use a custom implementation of PartialEq to handle
/// scalars of different, compatible orders.
#[derive(Debug,Clone)]
pub enum Scalar<T: Coeffs> {
    Exact(i32, T),
    Float(Complex<f64>),
}

/// Adds the ability to take non-integer types modulo 2.
pub trait Mod2 {
    fn mod2(&self) -> Self;
}

impl Mod2 for Rational {
    fn mod2(&self) -> Rational {
       Rational::new(*self.numer() % (2 * *self.denom()), *self.denom())
    }
}

/// Produce a number from rational root of -1.
pub trait FromPhase {
    fn from_phase(p: Rational) -> Self;
}

/// Contains the numbers sqrt(2) and 1/sqrt(2), often used for renormalisation of
/// qubit tensors and matrices.
pub trait Sqrt2: Sized {
    fn sqrt2() -> Self { Self::sqrt2_pow(1) }
    fn one_over_sqrt2() -> Self { Self::sqrt2_pow(-1) }
    fn sqrt2_pow(p: i32) -> Self;
}

/// A list of coefficients. We give this as a parameter to allow either
/// fixed-size lists (e.g. [i32;4]) or dynamic ones (e.g. [Vec]\<i32\>). Only
/// the former can be used in tensors and matrices, because they have to
/// implement Copy (i.e. size must be known at compile time).
pub trait Coeffs: Clone + std::ops::IndexMut<usize,Output=i32> {
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
/// when converting Scalar\<S\> to Scalar\<T\>, which is already implemented
/// as a noop for [From] when S = T.
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
                let rt2_pow = Complex::new(2f64,0f64).powf((*pow as f64) / 2f64);

                let mut num = Complex::new(0f64, 0f64);
                for i in 0..coeffs.len() {
                    num += (coeffs[i] as f64) * omega.powu(i as u32) * rt2_pow;
                }
                num
            },
            Float(c) => *c
        }
    }

    /// If zero, make sqrt(2) power 0, otherwise make it as big as possible
    pub fn reduced(&self) -> Scalar<T> {
        if let Exact(mut pow,mut coeffs) = self.clone() {
            let mut all_zero = true;
            for i in 0..coeffs.len() {
                all_zero = all_zero && coeffs[i] == 0;
            }

            if all_zero { return Exact(0, coeffs); }

            let mut red = true;
            while red {
                for i in 0..coeffs.len() { red = red && coeffs[i].is_multiple_of(&2); }

                if red {
                    for i in 0..coeffs.len() {
                        coeffs[i] = coeffs[i] / 2;
                    }
                    pow += 2;
                }
            }
            Exact(pow, coeffs)
        } else {
            self.clone()
        }
    }

    pub fn mul_sqrt2_pow(&mut self, p: i32) {
        match self {
            Exact(pow, _) => *pow += p,
            Float(c) => {
                let rt2_pow = f64::powi(f64::sqrt(2.0), p);
                c.re *= rt2_pow;
                c.im *= rt2_pow;
            }
        }
    }

    pub fn mul_phase(&mut self, phase: Rational) {
        (*self) *= Scalar::from_phase(phase);
    }

    pub fn to_float(&self) -> Scalar<T> {
        Float(self.float_value())
    }

    pub fn one_plus_phase(p: Rational) -> Scalar<T> {
        Scalar::one() + Scalar::from_phase(p)
    }
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
        Scalar::Exact(p, T::one())
    }
}

impl<T: Coeffs> FromPhase for Scalar<T> {
    fn from_phase(p: Rational) -> Scalar<T> {
        let mut rnumer = *p.numer();
        let mut rdenom = *p.denom();

        if rdenom < 0 {
            rnumer *= -1;
            rdenom *= -1;
        }

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
                Scalar::Exact(0, coeffs)
            },
            None => {
                let f = (*p.numer() as f64) / (*p.denom() as f64);
                Scalar::Float(Complex::new(-1.0,0.0).powf(f))
            }
        }
    }
}


impl<T: Coeffs> fmt::Display for Scalar<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Exact(pow,coeffs) => {
                write!(f, "rt(2)^{} (", pow)?;

                for i in 0..coeffs.len() {
                    if i == 0 {
                        write!(f, "{}", coeffs[0])?;
                    } else {
                        write!(f, " + {} * om^{}", coeffs[i], i)?;
                    }
                }

                write!(f, ")")
            },
            Float(c) => write!(f, "{}", c),
        }
    }
}

// The main implementation of the Mul trait uses references, so we don't need
// to make a copy of the scalars to multiply them.
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
                                    coeffs[pos - lcm] -= coeffs0[i] * coeffs1[j];
                                }
                            }
                        }

                        // TODO: we may not want to reduce at every multiplication
                        Exact(pow0 + pow1, coeffs).reduced()
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
    fn mul(self, rhs: Scalar<T>) -> Self::Output { &self * &rhs }
}

impl<'a, T: Coeffs> std::ops::Mul<Scalar<T>> for &'a Scalar<T> {
    type Output = Scalar<T>;
    fn mul(self, rhs: Scalar<T>) -> Self::Output { self * &rhs }
}

impl<'a, T: Coeffs> std::ops::Mul<&'a Scalar<T>> for Scalar<T> {
    type Output = Scalar<T>;
    fn mul(self, rhs: &Scalar<T>) -> Self::Output { &self * rhs }
}

impl<'a, T: Coeffs> std::ops::MulAssign<Scalar<T>> for Scalar<T> {
    fn mul_assign(&mut self, rhs: Scalar<T>) {
        *self = &*self * &rhs;
    }
}

impl<'a, T: Coeffs> std::ops::MulAssign<&'a Scalar<T>> for Scalar<T> {
    fn mul_assign(&mut self, rhs: &Scalar<T>) {
        *self = &*self * rhs;
    }
}

// The main implementation of the Add trait uses references, so we don't need
// to make a copy of the scalars to multiply them.
impl<'a, 'b, T: Coeffs> std::ops::Add<&'b Scalar<T>> for &'a Scalar<T> {
    type Output = Scalar<T>;

    fn add(self, rhs: &Scalar<T>) -> Self::Output {
        match (self,rhs) {
            (Float(c), x) => Float(c + x.float_value()),
            (x, Float(c)) => Float(x.float_value() + c),
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                if (pow0 - pow1) % 2 != 0 {
                    // if sqrt(2) powers don't have the same parity, we have to fall back on float repr
                    Float(self.float_value() + rhs.float_value())
                } else {
                    let (lcm, pad0, pad1) = lcm_with_padding(coeffs0.len(), coeffs1.len());

                    // if the powers of sqrt(2) are different by an even number, we need to rescale
                    // the coefficients of the scalar with the larger power
                    let (new_pow, scale0, scale1) =
                        if *pow0 > *pow1 { (*pow1, (1 as i32) << (((pow0 - pow1)/2) as u32), 1) }
                        else if *pow0 < *pow1 { (*pow0, 1, (1 as i32) << (((pow1 - pow0)/2) as u32)) }
                        else { (*pow0, 1, 1) };

                    match T::new(lcm) {
                        Some((mut coeffs, pad)) => {
                            for i in 0..coeffs0.len() {
                                coeffs[i*pad*pad0] += scale0*coeffs0[i];
                            }

                            for i in 0..coeffs1.len() {
                                coeffs[i*pad*pad1] += scale1*coeffs1[i];
                            }

                            // TODO: we may not want to reduce at every addition
                            Exact(new_pow, coeffs).reduced()
                        },
                        None => Float(self.float_value() + self.float_value())
                    }
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
        0.0000000001f64 //f64::default_epsilon()
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
        match (self.reduced(), other.reduced()) {
            (Float(c0), Float(c1)) => c0 == c1,
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
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
        impl Coeffs for [i32;$n] {
            fn len(&self) -> usize { $n }
            fn zero() -> [i32;$n] { [0;$n] }
            fn one() -> [i32;$n] { let mut a = [0;$n]; a[0] = 1; a }
            fn new(sz: usize) -> Option<([i32;$n],usize)> {
                if $n.is_multiple_of(&sz) {
                    Some(([0;$n], $n/sz))
                } else {
                    None
                }
            }
        }

        pub type $name = Scalar<[i32;$n]>;
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

// impl Coeffs for [i32;4] {
//     fn len(&self) -> usize { 4 }
//     fn zero() -> [i32;4] { [0;4] }
//     fn one() -> [i32;4] { let mut a = [0;4]; a[0] = 1; a }
//     fn new(sz: usize) -> Option<([i32;4],usize)> {
//         if (sz as i32).divides(&4) {
//             Some(([0;4], 4/sz))
//         } else {
//             None
//         }
//     }
// }

impl Coeffs for Vec<i32> {
    fn len(&self) -> usize { self.len() }
    fn zero() -> Vec<i32> { vec![0] }
    fn one() -> Vec<i32> { vec![1] }
    fn new(sz: usize) -> Option<(Vec<i32>,usize)> {
        Some((vec![0; sz],1))
    }
}

pub type ScalarN = Scalar<Vec<i32>>;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn approx_mul() {
        let s: Scalar<[i32;4]> = Scalar::real(f64::sqrt(0.3) * f64::sqrt(0.3) - 0.3);
        let t: Scalar<[i32;4]> = Scalar::zero();
        assert_ne!(s, t);
        assert_abs_diff_eq!(s, t);
    }

    #[test]
    fn sqrt_i() {
        let s = Scalar::Exact(0, [0, 1, 0, 0]);
        assert_abs_diff_eq!(s.to_float(), Scalar::complex(1.0 / f64::sqrt(2.0), 1.0 / f64::sqrt(2.0)));
    }

    #[test]
    fn mul_same_base() {
        let s = Scalar::Exact(0, [1, 2, 3, 4]);
        let t = Scalar::Exact(0, [4, 5, 6, 7]);
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
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(1,4)),  Scalar4::Exact(0, [0,1,0,0]));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(3,4)),  Scalar4::Exact(0, [0,0,0,1]));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(7,4)),  Scalar4::Exact(0, [0,0,0,-1]));
    }

    #[test]
    fn additions() {
        let s: ScalarN = Exact(0, vec![1,2,3,4]);
        let t: ScalarN = Exact(0, vec![2,3,4,5]);
        let st: ScalarN = Exact(0, vec![3,5,7,9]);
        assert_eq!(s + t, st);
    }

    #[test]
    fn additions_diff_base() {
        let s: ScalarN = Exact(2, vec![1,2,3,4]);
        let t: ScalarN = Exact(0, vec![2,3,4,5]);
        let st: ScalarN = Exact(0, vec![4,7,10,13]);
        assert_eq!(s + t, st);
    }

    #[test]
    fn reductions() {
        let s: ScalarN = Exact(-2, vec![2]);
        assert_eq!(s,s);
        assert_eq!(s, ScalarN::one());
    }

    // #[test]
    // fn one_plus_phases() {
    //     assert_abs_diff_eq!(Scalar::one_plus_phase(Rational::new(1,1)), Scalar::zero());

    //     let plus = Scalar::one_plus_phase(Rational::new(1,2));
    //     let minus = Scalar::one_plus_phase(Rational::new(-1,2));
    //     assert_abs_diff_eq!(plus * minus, Scalar::real(2.0));
    // }
}
