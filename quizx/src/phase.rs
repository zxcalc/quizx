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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalization() {
        // Values in range (-1, 1] are unchanged.
        assert_eq!(Phase::from((1, 2)).to_rational(), Rational64::new(1, 2));
        assert_eq!(Phase::from((-1, 2)).to_rational(), Rational64::new(-1, 2));
        assert_eq!(Phase::from((1, 4)).to_rational(), Rational64::new(1, 4));
        assert_eq!(Phase::from((1, 1)).to_rational(), Rational64::new(1, 1));

        // Values outside range are normalized by adding/subtracting 2.
        assert_eq!(Phase::from((3, 2)).to_rational(), Rational64::new(-1, 2));
        assert_eq!(Phase::from((5, 4)).to_rational(), Rational64::new(-3, 4));
        assert_eq!(Phase::from((-3, 2)).to_rational(), Rational64::new(1, 2));

        // Boundary: 1 is included, but not -1 (maps to 1).
        assert_eq!(Phase::from((2, 2)).to_rational(), Rational64::new(1, 1));
        assert_eq!(Phase::from((-2, 2)).to_rational(), Rational64::new(1, 1));
    }

    #[test]
    fn from_conversions() {
        assert_eq!(Phase::from(0_i64), Phase::zero());
        assert_eq!(Phase::from(1_i64), Phase::one());
        assert_eq!(Phase::from(2_i64), Phase::zero()); // normalized

        let p: Phase = (1, 4).into();
        assert_eq!(p.to_rational(), Rational64::new(1, 4));
    }

    #[test]
    fn f64_conversion() {
        let p = Phase::from_f64(0.5);
        assert!((p.to_f64() - 0.5).abs() < 1e-10);

        let p = Phase::from_f64(0.25);
        assert!((p.to_f64() - 0.25).abs() < 1e-10);
    }

    #[test]
    fn classification_predicates() {
        // Clifford phases: multiples of 1/2 (denominator <= 2).
        assert!(Phase::zero().is_clifford());
        assert!(Phase::one().is_clifford());
        assert!(Phase::from((1, 2)).is_clifford());
        assert!(Phase::from((-1, 2)).is_clifford());
        assert!(!Phase::from((1, 4)).is_clifford());
        assert!(!Phase::from((1, 3)).is_clifford());

        // Proper Clifford: exactly +/- 1/2.
        assert!(Phase::from((1, 2)).is_proper_clifford());
        assert!(Phase::from((-1, 2)).is_proper_clifford());
        assert!(!Phase::zero().is_proper_clifford());
        assert!(!Phase::one().is_proper_clifford());
        assert!(!Phase::from((1, 4)).is_proper_clifford());

        // Pauli: 0 or 1.
        assert!(Phase::zero().is_pauli());
        assert!(Phase::one().is_pauli());
        assert!(!Phase::from((1, 2)).is_pauli());

        // T-gate phases: denominator exactly 4.
        assert!(Phase::from((1, 4)).is_t());
        assert!(Phase::from((3, 4)).is_t());
        assert!(Phase::from((-1, 4)).is_t());
        assert!(!Phase::from((1, 2)).is_t()); // denominator 2
        assert!(!Phase::from((1, 8)).is_t()); // denominator 8
    }

    #[test]
    fn arithmetic() {
        let a = Phase::from((1, 4));
        let b = Phase::from((1, 4));
        assert_eq!(a + b, Phase::from((1, 2)));

        let a = Phase::from((3, 4));
        let b = Phase::from((1, 2));
        assert_eq!(a + b, Phase::from((-3, 4))); // 5/4 normalizes to -3/4

        let a = Phase::from((1, 2));
        let b = Phase::from((1, 4));
        assert_eq!(a - b, Phase::from((1, 4)));

        assert_eq!(-Phase::from((1, 4)), Phase::from((-1, 4)));
        assert_eq!(-Phase::one(), Phase::one()); // -1 normalizes to 1
    }

    #[test]
    fn mul_div() {
        let a = Phase::from((1, 4));
        assert_eq!(a * 2, Phase::from((1, 2)));
        assert_eq!(a * 4, Phase::one());
        assert_eq!(a * 8, Phase::zero()); // 2 normalizes to 0

        let a = Phase::from((1, 2));
        assert_eq!(a / 2, Phase::from((1, 4)));
    }

    #[test]
    fn assign_ops() {
        let mut p = Phase::from((1, 4));
        p += Phase::from((1, 4));
        assert_eq!(p, Phase::from((1, 2)));

        let mut p = Phase::from((1, 2));
        p -= Phase::from((1, 4));
        assert_eq!(p, Phase::from((1, 4)));

        let mut p = Phase::from((1, 4));
        p *= 2;
        assert_eq!(p, Phase::from((1, 2)));

        let mut p = Phase::from((1, 2));
        p /= 2;
        assert_eq!(p, Phase::from((1, 4)));
    }

    #[test]
    fn zero_one_traits() {
        assert!(Phase::zero().is_zero());
        assert!(!Phase::one().is_zero());
        assert!(Phase::one().is_one());
        assert!(!Phase::zero().is_one());
    }

    #[test]
    fn limit_denominator_phase() {
        // High-precision phase approximated to simpler form.
        let p = Phase::from((355, 113)); // approximation of pi
        let approx = p.limit_denominator(10);
        assert!(approx.to_rational().denom().abs() <= 10);
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}", Phase::from((1, 2))), "1/2");
        assert_eq!(format!("{}", Phase::zero()), "0");
    }
}
