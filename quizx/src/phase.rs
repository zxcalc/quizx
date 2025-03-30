//! Phase encoded as either rational or floating point number of half-turns.

pub mod utils;

use std::fmt::{self, Display};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num::{FromPrimitive, One, Rational64, ToPrimitive, Zero};

use utils::limit_denominator;

/// A phase, expressed in half-turns and encoded as a rational number.
///
/// The phase is always normalized to be in the range (-1,1].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Phase {
    r: Rational64,
}

impl Phase {
    /// Creates a new phase.
    ///
    /// Normalizes the phase to be in the range (-1,1].
    pub fn new(r: impl Into<Rational64>) -> Self {
        Self { r: r.into() }.normalize()
    }

    /// Returns the phase as a rational number.
    pub fn to_rational(&self) -> Rational64 {
        self.r
    }

    /// Creates a new phase from a floating point number of half-turns.
    ///
    /// Rounds the floating point number to a rational number and
    /// normalizes it to be in the range (-1,1].
    pub fn from_f64(f: f64) -> Self {
        Self::new(Rational64::from_f64(f).unwrap())
    }

    /// Returns the phase as a floating point number of half-turns.
    pub fn to_f64(&self) -> f64 {
        self.r.to_f64().unwrap()
    }

    /// Normalizes the phase to be in the range (-1,1] by adding or subtracting multiples of 2.
    pub fn normalize(&self) -> Phase {
        let denom = *self.r.denom();
        let mut num = *self.r.numer();
        if -denom < num && num <= denom {
            return *self;
        }
        num = num.rem_euclid(2 * denom);
        if num > *self.r.denom() {
            num -= 2 * denom;
        }
        Rational64::new(num, denom).into()
    }

    /// Returns `true` if the phase is a multiple of 1/2.
    pub fn is_clifford(&self) -> bool {
        self.r.denom().abs() <= 2
    }

    /// Returns `true` if the phase is either -1/2 or 1/2.
    pub fn is_proper_clifford(&self) -> bool {
        self.r == Rational64::new(1, 2) || self.r == Rational64::new(-1, 2)
    }

    /// Returns `true` if the phase is 0 or 1.
    pub fn is_pauli(&self) -> bool {
        self.is_zero() || self.is_one()
    }

    /// Returns `true` if the phase a non-clifford multiple of 1/4.
    pub fn is_t(&self) -> bool {
        self.r.denom().abs() == 4
    }

    /// Approximate a phase's fraction to a Rational64 number with a small denominator.
    ///
    /// # Panics
    ///
    /// Panics if `max_denom` is 0.
    pub fn limit_denominator(&self, max_denom: i64) -> Self {
        Self::new(limit_denominator(self.r, max_denom))
    }
}

impl Display for Phase {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.r)
    }
}

impl From<Rational64> for Phase {
    fn from(r: Rational64) -> Phase {
        Phase::new(r)
    }
}

impl From<f64> for Phase {
    fn from(f: f64) -> Phase {
        Phase::from_f64(f)
    }
}

impl From<Phase> for Rational64 {
    fn from(phase: Phase) -> Rational64 {
        phase.to_rational()
    }
}

impl From<Phase> for f64 {
    fn from(phase: Phase) -> f64 {
        phase.to_f64()
    }
}

impl From<i64> for Phase {
    fn from(i: i64) -> Phase {
        Phase::new(Rational64::from_i64(i).unwrap())
    }
}

impl From<(i64, i64)> for Phase {
    fn from(i: (i64, i64)) -> Phase {
        let r: Rational64 = i.into();
        Phase::new(r)
    }
}

impl Zero for Phase {
    fn zero() -> Self {
        Phase::new(Rational64::zero())
    }

    fn is_zero(&self) -> bool {
        self.r.is_zero()
    }
}

impl One for Phase {
    fn one() -> Self {
        Phase::new(Rational64::one())
    }

    fn is_one(&self) -> bool {
        self.r.is_one()
    }
}

impl Neg for Phase {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.r)
    }
}

impl Add for Phase {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(self.r + other.r)
    }
}

impl AddAssign for Phase {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for Phase {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::new(self.r - other.r)
    }
}

impl SubAssign for Phase {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Mul for Phase {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self::new(self.r * other.r)
    }
}

impl Mul<i64> for Phase {
    type Output = Self;

    fn mul(self, other: i64) -> Self {
        Self::new(self.r * other)
    }
}

impl MulAssign for Phase {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl MulAssign<i64> for Phase {
    fn mul_assign(&mut self, other: i64) {
        *self = *self * other;
    }
}

impl Div for Phase {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Self::new(self.r / other.r)
    }
}

impl Div<i64> for Phase {
    type Output = Self;

    fn div(self, other: i64) -> Self {
        Self::new(self.r / other)
    }
}

impl DivAssign for Phase {
    fn div_assign(&mut self, other: Self) {
        *self = *self / other;
    }
}

impl DivAssign<i64> for Phase {
    fn div_assign(&mut self, other: i64) {
        *self = *self / other;
    }
}
