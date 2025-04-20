use approx::AbsDiffEq;
use num::complex::Complex;
pub use num::traits::identities::{One, Zero};
use num::{integer, Integer, Rational64};
use std::cmp::min;
use std::f64::consts::PI;
use std::fmt;
use std::ops::{Add, AddAssign, Mul};

#[derive(Debug, Clone, Copy)]
pub struct FScalar {
    c: [f64; 4]
}

impl Add<&FScalar> for &FScalar {
    type Output = FScalar;

    fn add(self, rhs: &FScalar) -> Self::Output {
        FScalar { c: [
            self.c[0] + rhs.c[0],
            self.c[1] + rhs.c[1],
            self.c[2] + rhs.c[2],
            self.c[3] + rhs.c[3]
            ] }
    }
}

// These 3 variations take ownership of one or both args
impl Add<FScalar> for FScalar {
    type Output = FScalar;
    fn add(self, rhs: FScalar) -> Self::Output {
        &self + &rhs
    }
}

impl Add<FScalar> for &FScalar {
    type Output = FScalar;
    fn add(self, rhs: FScalar) -> Self::Output {
        self + &rhs
    }
}

impl Add<&FScalar> for FScalar {
    type Output = FScalar;
    fn add(self, rhs: &FScalar) -> Self::Output {
        &self + rhs
    }
}

impl Mul<&FScalar> for &FScalar {
    type Output = FScalar;

    fn mul(self, rhs: &FScalar) -> Self::Output {
        let mut c = [0.0, 0.0, 0.0, 0.0];
        for i in 0..4 {
            for j in 0..4 {
                let pos = (i * j) % 8;
                if pos < 4 {
                    c[pos] += self.c[i] * rhs.c[j];
                } else {
                    c[pos - 4] += -self.c[i] * rhs.c[j];
                }
            }
        }
        FScalar { c }
    }
}

// These 3 variations take ownership of one or both args
impl Mul for FScalar {
    type Output = FScalar;
    fn mul(self, rhs: FScalar) -> Self::Output {
        &self * &rhs
    }
}
impl Mul<FScalar> for &FScalar {
    type Output = FScalar;
    fn mul(self, rhs: FScalar) -> Self::Output {
        self * &rhs
    }
}
impl Mul<&FScalar> for FScalar {
    type Output = FScalar;
    fn mul(self, rhs: &FScalar) -> Self::Output {
        &self * rhs
    }
}

impl Zero for FScalar {
    fn zero() -> Self {
        FScalar { c: [0.0, 0.0, 0.0, 0.0] }
    }

    fn is_zero(&self) -> bool {
        self.c == [0.0, 0.0, 0.0, 0.0]
    }
}

impl One for FScalar {
    fn one() -> Self {
        FScalar { c: [1.0, 0.0, 0.0, 0.0] }
    }

    fn is_one(&self) -> bool {
        self.c == [1.0, 0.0, 0.0, 0.0]
    }
}