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

use approx::AbsDiffEq;
use num::complex::Complex;
pub use num::traits::identities::{One, Zero};
use num::{Rational64, ToPrimitive};
use std::f64::consts::{PI, SQRT_2};
use std::fmt;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use crate::phase::Phase;
use crate::scalar::dyadic::DyadicExponentOverflowError;
pub use crate::scalar_traits::{FromPhase, Sqrt2};
pub mod dyadic;
pub use dyadic::Dyadic;

/// This is the main representation of scalars used in QuiZX. It is a wrapper around
/// four [`Dyadic`] rationals, used to represent the coefficients in a complex number of
/// the form:
///
///   a + b ω + c ω² + d ω³
///
/// where ω is the 4th root of -1, i.e. exp(i π/4).
///
/// When scalars come from a Clifford+T circuit or ZX diagram, i.e. where all phase
/// angles are integer multiples of π/4, scalars are guaranteed to be of the form
/// above, with rational coefficients of the form x/2^y, for integers x, y.
///
/// Note that ω² = i, and [`Dyadic`] can store any [`f64`] losslessly. Hence, for all
/// other complex numbers, [`Scalar4`] gives an approximate representation of that number
/// as a + c ω² = a + i c. This allows easy conversions to and from [`Complex<f64>`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scalar4([Dyadic; 4]);

impl Scalar4 {
    /// Constructs a scalar from 4 integer coefficients and a power of 2. For
    /// `coeffs := [a,b,c,d]`, the resulting scalar is: `2^pow * (a + b ω + c ω² + d ω³)`.
    pub fn new(coeffs: [i64; 4], pow: i32) -> Self {
        Scalar4(coeffs.map(|coeff| Dyadic::new(coeff, pow)))
    }

    /// Constructs a scalar as a real number with the given float value
    pub fn real(r: f64) -> Self {
        Scalar4([r.into(), 0.into(), 0.into(), 0.into()])
    }

    /// Constructs a scalar as a complex number `re + im * i`.
    pub fn complex(re: f64, im: f64) -> Self {
        Scalar4([re.into(), 0.into(), im.into(), 0.into()])
    }

    /// Constructs a scalar `1 + e^(i π * phase)` for the given `Phase`
    pub fn one_plus_phase(phase: impl Into<Phase>) -> Self {
        Scalar4::one() + Scalar4::from_phase(phase)
    }

    /// Convience method for multipling by the given power of sqrt(2)
    pub fn mul_sqrt2_pow(&mut self, p: i32) {
        *self *= Scalar4::sqrt2_pow(p);
    }

    /// Convience method for multipling by `e^(i π * phase)` for the given `Phase`
    pub fn mul_phase(&mut self, phase: impl Into<Phase>) {
        *self *= Scalar4::from_phase(phase);
    }

    /// Convience method for multipling by `1 + e^(i π * phase)` for the given `Phase`
    pub fn mul_one_plus_phase(&mut self, phase: impl Into<Phase>) {
        *self *= Scalar4::one_plus_phase(phase)
    }

    /// Returns the complex conjugate of the scalar
    pub fn conj(&self) -> Self {
        Scalar4([self.0[0], -self.0[3], -self.0[2], -self.0[1]])
    }

    /// Converts `Scalar4` into a `Complex<f64>`
    pub fn complex_value(&self) -> Complex<f64> {
        self.try_into().unwrap()
    }

    pub fn approx(&self) -> bool {
        self.0.iter().any(|c| c.approx())
    }

    pub fn set_approx(&mut self, approx: bool) {
        self.0.iter_mut().for_each(|c| c.set_approx(approx));
    }

    /// The number of non-zero coefficients
    fn num_coeffs(&self) -> usize {
        self.0
            .iter()
            .fold(0, |n, &co| if !co.is_zero() { n + 1 } else { n })
    }

    /// If the scalar is an exact representation of a number in the form `exp(i k π / 4) * sqrt(2)^p`,
    /// return `Some((k/4, p))`, otherwise return `None`.
    pub fn exact_phase_and_sqrt2_pow(&self) -> Option<(Phase, i32)> {
        let mut s: Scalar4 = *self;

        let p;
        if s.num_coeffs() != 1 {
            p = -1;
            s *= Self::sqrt2();
            if s.num_coeffs() != 1 {
                return None;
            }
        } else {
            p = 0;
        }

        let (i, c) = s.0.iter().enumerate().find(|(_, c)| !c.is_zero()).unwrap();

        if c.val() == 1 {
            Some((Phase::new(Rational64::new(i as i64, 4)), c.exp() * 2 + p))
        } else if c.val() == -1 {
            Some((
                Phase::new(Rational64::new((i + 4) as i64, 4)),
                c.exp() * 2 + p,
            ))
        } else {
            None
        }
    }
}

impl Default for Scalar4 {
    fn default() -> Self {
        Self::zero()
    }
}

impl Sqrt2 for Scalar4 {
    fn sqrt2_pow(p: i32) -> Self {
        Scalar4(if p % 2 == 0 {
            [Dyadic::new(1, p / 2), 0.into(), 0.into(), 0.into()]
        } else {
            let d = Dyadic::new(1, (p - 1) / 2);
            [0.into(), d, 0.into(), -d]
        })
    }
}

impl FromPhase for Scalar4 {
    fn from_phase(phase: impl Into<Phase>) -> Self {
        let p: Phase = phase.into();
        p.into()
    }

    fn minus_one() -> Self {
        Scalar4([-1, 0, 0, 0].map(|c| c.into()))
    }
}

impl ndarray::ScalarOperand for Scalar4 {}

impl fmt::Display for Scalar4 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut fst = true;

        for (i, dy) in self.0.iter().enumerate() {
            let (mut v, mut e) = dy.val_and_exp();
            if v > -1024 && v < 1024 && e > 0 && e <= 10 {
                v *= 2i64.pow(e as u32);
                e = 0;
            }

            if v != 0 {
                if fst {
                    write!(f, "{v}")?;
                    fst = false;
                } else if v > 0 {
                    write!(f, " + {v}")?;
                } else {
                    write!(f, " - {}", -v)?;
                }

                if e != 0 {
                    write!(f, "e{e}")?;
                }
                if i == 1 {
                    write!(f, " ω")?;
                }
                if i == 2 {
                    write!(f, " ω²")?;
                }
                if i == 3 {
                    write!(f, " ω³")?;
                }
            }
        }

        if fst {
            write!(f, "0")?;
        }

        Ok(())
    }
}

impl Add<&Scalar4> for &Scalar4 {
    type Output = Scalar4;
    fn add(self, rhs: &Scalar4) -> Self::Output {
        Scalar4([
            self.0[0] + rhs.0[0],
            self.0[1] + rhs.0[1],
            self.0[2] + rhs.0[2],
            self.0[3] + rhs.0[3],
        ])
    }
}

impl Sum<Scalar4> for Scalar4 {
    fn sum<I: Iterator<Item = Scalar4>>(iter: I) -> Self {
        iter.fold(Scalar4::zero(), |accumulator, current_fscalar| {
            accumulator + current_fscalar // Uses Scalar4 + Scalar4
        })
    }
}

// These 3 variations take ownership of one or both args
impl Add<Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn add(self, rhs: Scalar4) -> Self::Output {
        let rself = &self;
        let rrhs = &rhs;
        rself + rrhs
    }
}

impl Add<Scalar4> for &Scalar4 {
    type Output = Scalar4;
    fn add(self, rhs: Scalar4) -> Self::Output {
        let rrhs = &rhs;
        self + rrhs
    }
}

impl Add<&Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn add(self, rhs: &Scalar4) -> Self::Output {
        let rself = &self;
        rself + rhs
    }
}

impl AddAssign for Scalar4 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl AddAssign<&Scalar4> for Scalar4 {
    fn add_assign(&mut self, rhs: &Scalar4) {
        *self = *self + rhs;
    }
}

impl Sub<&Scalar4> for &Scalar4 {
    type Output = Scalar4;
    fn sub(self, rhs: &Scalar4) -> Self::Output {
        Scalar4([
            self.0[0] - rhs.0[0],
            self.0[1] - rhs.0[1],
            self.0[2] - rhs.0[2],
            self.0[3] - rhs.0[3],
        ])
    }
}

impl Sub<Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn sub(self, rhs: Scalar4) -> Self::Output {
        &self - &rhs
    }
}

impl Sub<Scalar4> for &Scalar4 {
    type Output = Scalar4;
    fn sub(self, rhs: Scalar4) -> Self::Output {
        let rrhs = &rhs;
        self - rrhs
    }
}

impl Sub<&Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn sub(self, rhs: &Scalar4) -> Self::Output {
        let rself = &self;
        rself - rhs
    }
}

impl SubAssign for Scalar4 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl SubAssign<&Scalar4> for Scalar4 {
    fn sub_assign(&mut self, rhs: &Scalar4) {
        *self = *self - rhs;
    }
}

impl Mul<&Scalar4> for &Scalar4 {
    type Output = Scalar4;
    fn mul(self, rhs: &Scalar4) -> Self::Output {
        let mut s = Scalar4::zero();
        for i in 0..4 {
            if !self.0[i].is_zero() {
                for j in 0..4 {
                    let pos = (i + j) % 8;
                    if pos < 4 {
                        s.0[pos] += self.0[i] * rhs.0[j];
                    } else {
                        s.0[pos - 4] += -self.0[i] * rhs.0[j];
                    }
                }
            }
        }
        s
    }
}

// These 3 variations take ownership of one or both args
impl Mul<Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn mul(self, rhs: Scalar4) -> Self::Output {
        let rself = &self;
        let rrhs = &rhs;
        rself * rrhs
    }
}
impl Mul<Scalar4> for &Scalar4 {
    type Output = Scalar4;
    fn mul(self, rhs: Scalar4) -> Self::Output {
        let rrhs = &rhs;
        self * rrhs
    }
}
impl Mul<&Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn mul(self, rhs: &Scalar4) -> Self::Output {
        let rself = &self;
        rself * rhs
    }
}

/// Implements *=
impl MulAssign<Scalar4> for Scalar4 {
    fn mul_assign(&mut self, rhs: Scalar4) {
        *self = *self * rhs;
    }
}

// Variation takes ownership of rhs
impl MulAssign<&Scalar4> for Scalar4 {
    fn mul_assign(&mut self, rhs: &Scalar4) {
        *self = *self * rhs;
    }
}

impl Product<Scalar4> for Scalar4 {
    fn product<I: Iterator<Item = Scalar4>>(iter: I) -> Self {
        iter.fold(Scalar4::one(), |accumulator, current_scalar| {
            accumulator * current_scalar // Uses Scalar4 * Scalar4
        })
    }
}

impl Zero for Scalar4 {
    fn zero() -> Self {
        [0, 0, 0, 0].into()
    }

    fn is_zero(&self) -> bool {
        self.0.iter().all(|c| c.is_zero())
    }
}

impl One for Scalar4 {
    fn one() -> Self {
        [1, 0, 0, 0].into()
    }

    fn is_one(&self) -> bool {
        *self == Self::one()
    }
}

impl AbsDiffEq<Scalar4> for Scalar4 {
    type Epsilon = <f64 as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        1e-10f64
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let c1: Complex<f64> = self.try_into().unwrap();
        let c2: Complex<f64> = other.try_into().unwrap();
        f64::abs_diff_eq(&c1.re, &c2.re, epsilon) && f64::abs_diff_eq(&c1.im, &c2.im, epsilon)
    }
}

impl From<[i64; 4]> for Scalar4 {
    fn from(value: [i64; 4]) -> Self {
        Scalar4(value.map(|c| c.into()))
    }
}

impl From<[f64; 4]> for Scalar4 {
    fn from(value: [f64; 4]) -> Self {
        Scalar4(value.map(|c| c.into()))
    }
}

impl From<i64> for Scalar4 {
    fn from(value: i64) -> Self {
        [value, 0, 0, 0].into()
    }
}

impl From<f64> for Scalar4 {
    fn from(value: f64) -> Self {
        [value, 0.0, 0.0, 0.0].into()
    }
}

impl From<Complex<f64>> for Scalar4 {
    fn from(value: Complex<f64>) -> Self {
        [value.re, 0.0, value.im, 0.0].into()
    }
}

impl TryFrom<&Scalar4> for Complex<f64> {
    type Error = DyadicExponentOverflowError;
    fn try_from(value: &Scalar4) -> Result<Self, Self::Error> {
        Ok(Complex {
            re: f64::try_from(value.0[0])? + f64::try_from(value.0[1] - value.0[3])? * 0.5 * SQRT_2,
            im: f64::try_from(value.0[2])? + f64::try_from(value.0[1] + value.0[3])? * 0.5 * SQRT_2,
        })
    }
}

impl TryFrom<Scalar4> for Complex<f64> {
    type Error = DyadicExponentOverflowError;
    fn try_from(value: Scalar4) -> Result<Self, Self::Error> {
        Ok(Complex {
            re: f64::try_from(value.0[0])? + f64::try_from(value.0[1] - value.0[3])? * 0.5 * SQRT_2,
            im: f64::try_from(value.0[2])? + f64::try_from(value.0[1] + value.0[3])? * 0.5 * SQRT_2,
        })
    }
}

impl From<Phase> for Scalar4 {
    fn from(value: Phase) -> Self {
        let r: Rational64 = value.into();
        if 4 % r.denom() == 0 {
            let pos = (r.numer() * (4 / r.denom())).rem_euclid(8) as usize;
            let mut c: [i64; 4] = [0, 0, 0, 0];
            if pos >= 4 {
                c[pos - 4] = -1;
            } else {
                c[pos] = 1;
            }
            c.into()
        } else {
            let angle = PI * r.to_f64().unwrap();
            [f64::cos(angle), 0.0, f64::sin(angle), 0.0].into()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use rstest::rstest;

    #[test]
    fn display() {
        let s = Scalar4::zero();
        assert_eq!(format!("{s}"), "0");

        let s: Scalar4 = [1, 2, 3, 4].into();
        assert_eq!(format!("{s}"), "1 + 2 ω + 3 ω² + 4 ω³");

        let s: Scalar4 = [-1, -2, -3, -4].into();
        assert_eq!(format!("{s}"), "-1 - 2 ω - 3 ω² - 4 ω³");

        let s: Scalar4 = [0, 2, 0, 4].into();
        assert_eq!(format!("{s}"), "2 ω + 4 ω³");

        let s: Scalar4 = [0.5, 0.25, 0.125, 0.0625].into();
        assert_eq!(format!("{s}"), "1e-1 + 1e-2 ω + 1e-3 ω² + 1e-4 ω³");

        let s: Scalar4 = [
            2.0f64.powi(11),
            2.0f64.powi(21),
            2.0f64.powi(31),
            2.0f64.powi(41),
        ]
        .into();
        assert_eq!(format!("{s}"), "1e11 + 1e21 ω + 1e31 ω² + 1e41 ω³");
    }

    #[test]
    fn int_arith() {
        let s4: Scalar4 = 4.into();
        let s10: Scalar4 = 10.into();
        let s14: Scalar4 = 14.into();
        let sm1: Scalar4 = (-1).into();
        let s40: Scalar4 = 40.into();
        let sm14: Scalar4 = (-14).into();
        let sm40: Scalar4 = (-40).into();

        assert_eq!(s4 + s10, s14);
        assert_eq!(s4 * s10, s40);
        assert_eq!(sm1 * sm1, Scalar4::one());
        assert_eq!(sm1 * s40, sm40);
        assert_eq!(sm1 * (s4 + s10), sm14);
        assert_eq!(sm1 * s4 + sm1 * s10, sm14);
    }

    #[test]
    fn real_arith() {
        let a = 4.3;
        let b = 2e-11;
        let c = 0.3333333;
        let d = -1000000.0;
        let sa: Scalar4 = a.into();
        let sb: Scalar4 = b.into();
        let sc: Scalar4 = c.into();
        let sd: Scalar4 = d.into();

        assert_abs_diff_eq!(sa * sa, (a * a).into());
        assert_abs_diff_eq!(sa * sb, (a * b).into());
        assert_abs_diff_eq!(sc + (sa * sb), (c + (a * b)).into());
        assert_abs_diff_eq!(sd - (sa * sb), (d - (a * b)).into());
    }

    #[test]
    fn complex_arith() {
        let one: Scalar4 = 1.into();
        let i: Scalar4 = [0, 0, 1, 0].into();
        let om: Scalar4 = [0, 1, 0, 0].into();
        let sqrt2 = Scalar4::sqrt2_pow(1);
        assert_eq!(om * om, i);
        assert_eq!((one + i) * (one + i).conj(), one + one);
        assert_eq!(om + om.conj(), sqrt2);

        let c1: Complex<f64> = (one + i + i).complex_value();
        let c2: Complex<f64> = Complex::new(1.0, 2.0);
        assert_eq!(c1, c2);

        let c1: Complex<f64> = (Scalar4::sqrt2_pow(3) + i * sqrt2).complex_value();
        let c2: Complex<f64> = Complex::new(2.0 * SQRT_2, SQRT_2);
        assert_abs_diff_eq!(c1.re, c2.re);
        assert_abs_diff_eq!(c1.im, c2.im);
    }

    #[test]
    fn sqrt2() {
        // n.b. exact equality in the sqrt2_pow tests. This may break if conversion to Complex is
        // implemented differently.
        let sqrt2 = Scalar4::sqrt2_pow(1);
        let sqrt2_c: Complex<f64> = sqrt2.complex_value();
        assert_eq!(sqrt2_c.re, SQRT_2);
        assert_eq!(sqrt2_c.im, 0.0);

        let sqrt2_pow = Scalar4::sqrt2_pow(7);
        let sqrt2_pow_c: Complex<f64> = sqrt2_pow.complex_value();
        assert_eq!(sqrt2_pow_c.re, 2.0 * 2.0 * 2.0 * SQRT_2);
        assert_eq!(sqrt2_pow_c.im, 0.0);

        let two: Scalar4 = 2.into();
        let a = Scalar4::sqrt2_pow(10);
        let b = Scalar4::sqrt2_pow(11);
        let c: Scalar4 = 32.into();
        assert_eq!(two, sqrt2 * sqrt2);
        assert_eq!(a * sqrt2, b);
        assert_eq!(a, c);
    }

    #[test]
    fn from_t_phase() {
        let s: Scalar4 = Phase::new(Rational64::new(0, 1)).into();
        let s1: Scalar4 = [1, 0, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(1, 4)).into();
        let s1: Scalar4 = [0, 1, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(1, 2)).into();
        let s1: Scalar4 = [0, 0, 1, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(3, 4)).into();
        let s1: Scalar4 = [0, 0, 0, 1].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(1, 1)).into();
        let s1: Scalar4 = [-1, 0, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(5, 4)).into();
        let s1: Scalar4 = [0, -1, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(3, 2)).into();
        let s1: Scalar4 = [0, 0, -1, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(7, 4)).into();
        let s1: Scalar4 = [0, 0, 0, -1].into();
        assert_eq!(s, s1);

        let s: Scalar4 = Phase::new(Rational64::new(0, 1)).into();
        let s1: Scalar4 = [1, 0, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(-7, 4)).into();
        let s1: Scalar4 = [0, 1, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(-3, 2)).into();
        let s1: Scalar4 = [0, 0, 1, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(-5, 4)).into();
        let s1: Scalar4 = [0, 0, 0, 1].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(-1, 1)).into();
        let s1: Scalar4 = [-1, 0, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(-3, 4)).into();
        let s1: Scalar4 = [0, -1, 0, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(-1, 2)).into();
        let s1: Scalar4 = [0, 0, -1, 0].into();
        assert_eq!(s, s1);
        let s: Scalar4 = Phase::new(Rational64::new(-1, 4)).into();
        let s1: Scalar4 = [0, 0, 0, -1].into();
        assert_eq!(s, s1);
    }

    #[test]
    fn from_gen_phase() {
        let r = Rational64::new(-5, 37);
        let s1: Scalar4 = Phase::new(r).into();
        let c = Complex::new(0.0, PI * r.to_f64().unwrap()).exp();
        let s2: Scalar4 = c.into();

        assert_eq!(s1, s2);
        let r = Rational64::new(12, 117);
        let s1: Scalar4 = Phase::new(r).into();
        let c = Complex::new(0.0, PI * r.to_f64().unwrap()).exp();
        let s2: Scalar4 = c.into();
        assert_eq!(s1, s2);
    }

    #[rstest]
    #[case(Scalar4::zero())]
    #[case(Scalar4::one())]
    #[case(Scalar4::from_phase(1))]
    #[case(Scalar4::from_phase((1,2)))]
    #[case(Scalar4::from_phase((-1,2)))]
    #[case(Scalar4::real(2.0))]
    #[case(Scalar4::complex(1.0, 1.0))]
    #[case(Scalar4::new([0, 1, 0, -1], 3))]
    #[case(Scalar4::new([0, 7, 0, 7], -2))]
    #[case(Scalar4::new([-2, 0, -2, 0], 0))]
    #[case(Scalar4::new([2, 0, -2, 0], 30))]
    #[case(Scalar4::new([2, 0, 0, 0], -10))]
    fn scalar4_roundtrip(#[case] scalar: Scalar4) {
        let complex = scalar.complex_value();
        println!("complex value = {}", complex);
        let scalar1 = Scalar4::from(complex);
        assert_abs_diff_eq!(scalar, scalar1);
    }

    #[test]
    fn exact_phases() {
        let sqrt2 = Scalar4::sqrt2();
        let phase = Scalar4::new([0, 1, 0, 0], 0);

        let s = phase;
        println!("{}", s);
        assert_eq!(
            s.exact_phase_and_sqrt2_pow(),
            Some((Rational64::new(1, 4).into(), 0))
        );

        let s = phase * sqrt2;
        println!("{}", s);
        assert_eq!(
            s.exact_phase_and_sqrt2_pow(),
            Some((Rational64::new(1, 4).into(), 1))
        );

        let s = phase * sqrt2 * sqrt2;
        println!("{}", s);
        assert_eq!(
            s.exact_phase_and_sqrt2_pow(),
            Some((Rational64::new(1, 4).into(), 2))
        );

        let s = phase * phase * sqrt2;
        println!("{}", s);
        assert_eq!(
            s.exact_phase_and_sqrt2_pow(),
            Some((Rational64::new(1, 2).into(), 1))
        );
    }
}
