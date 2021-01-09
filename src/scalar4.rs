use crate::scalar::*;
use num::complex::Complex;
use num::rational::Rational;
pub use num::traits::identities::{Zero,One};
use std::fmt;
use approx::AbsDiffEq;

/// A type for exact and approximate representation of Scalars.
/// Note that '==' is only reliable when the scalar is Exact
/// and N is a power of 2.

#[derive(Debug,Clone,Copy,PartialEq)]
pub enum Scalar4 {
    // An exact representation of a scalar, which is given as
    // a power of sqrt(2) and an element of Z[omega], where
    // omega is the 2N-th root of unity, represented by its
    // first N coefficients.
    Exact(i32, [i32; 4]),
    // A floating-point representation of a scalar. We should
    // fall back to this if N^2 gets too big.
    Float(Complex<f64>),
}

use Scalar4::{Exact,Float};

impl Scalar4 {
    pub fn complex(re: f64, im: f64) -> Scalar4 {
        Float(Complex::new(re, im))
    }

    pub fn real(re: f64) -> Scalar4 {
        Float(Complex::new(re, 0.0))
    }

    pub fn float_value(&self) -> Complex<f64> {
        match self {
            Exact(pow, coeffs) => {
                let omega = Complex::new(-1f64, 0f64).powf(1f64 / (coeffs.len() as f64));
                let rt2_pow = Complex::new(2f64,0f64).powf((*pow as f64) / 2f64);

                let mut num = Complex::new(0f64, 0f64);
                for (i, c) in coeffs.iter().enumerate() {
                    num += (*c as f64) * omega.powu(i as u32) * rt2_pow;
                }
                num
            },
            Float(c) => *c
        }
    }

    pub fn mul_rt2_pow(&mut self, p: i32) {
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
        (*self) *= Scalar4::from_phase(phase);
    }

    pub fn to_float(&self) -> Scalar4 {
        Float(self.float_value())
    }

    pub fn one_plus_phase(p: Rational) -> Scalar4 {
        let mut s = Scalar4::from_phase(p);
        if let Scalar4::Exact(_,ref mut coeffs) = s { coeffs[0] += 1; }
        s
    }

    pub fn rt2_pow(p: i32) -> Scalar4 {
        Scalar4::Exact(p, [1,0,0,0])
    }
}

impl ndarray::ScalarOperand for Scalar4 { }

impl Zero for Scalar4 {
    fn zero() -> Scalar4 {
        Exact(0, [0; 4])
    }

    fn is_zero(&self) -> bool {
        *self == Scalar4::zero()
    }
}

impl One for Scalar4 {
    fn one() -> Scalar4 {
        Exact(0, [1, 0, 0, 0])
    }

    fn is_one(&self) -> bool {
        *self == Scalar4::one()
    }
}

impl Sqrt2 for Scalar4 {
    fn sqrt2() -> Scalar4 { Scalar4::rt2_pow(1) }
    fn one_over_sqrt2() -> Scalar4 { Scalar4::rt2_pow(-1) }
}

impl FromPhase for Scalar4 {
    fn from_phase(p: Rational) -> Scalar4 {
        let mut rnumer = *p.numer();
        let mut rdenom = *p.denom();

        if rdenom < 0 {
            rnumer *= -1;
            rdenom *= -1;
        }

        if rdenom == 1 || rdenom == 2 || rdenom == 4 {
            let f = 4 / rdenom;
            rdenom *= f;
            rnumer *= f;
            rnumer = rnumer.rem_euclid(2 * rdenom);
            let sgn = if rnumer >= rdenom {
                rnumer = rnumer - rdenom;
                -1
            } else {
                1
            };

            let mut coeffs = [0; 4];
            coeffs[rnumer as usize] = sgn;

            Scalar4::Exact(0, coeffs)
        } else {
            // TODO
            Scalar4::Float(Complex::zero())
        }
    }
}


impl fmt::Display for Scalar4 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Exact(pow,coeffs) => {
                write!(f, "rt(2)^{} (", pow)?;

                for (i,c) in coeffs.iter().enumerate() {
                    if i == 0 {
                        write!(f, "{}", c)?;
                    } else {
                        write!(f, " + {} * om^{}", c, i)?;
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
impl<'a, 'b> std::ops::Mul<&'b Scalar4> for &'a Scalar4 {
    type Output = Scalar4;

    fn mul(self, rhs: &Scalar4) -> Self::Output {
        match (self,rhs) {
            (Float(c), x) => Float(c * x.float_value()),
            (x, Float(c)) => Float(x.float_value() * c),
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                let mut coeffs = [0; 4];

                for (i,x) in coeffs0.iter().enumerate() {
                    for (j,y) in coeffs1.iter().enumerate() {
                        let pos = (i + j).rem_euclid(8);
                        if pos < 4 {
                            coeffs[pos] += x * y;
                        } else {
                            coeffs[pos - 4] -= x * y;
                        }
                    }
                }

                Exact(pow0 + pow1, coeffs)
            },
        }
    }
}

// These 3 variations take ownership of one or both args
impl std::ops::Mul<Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn mul(self, rhs: Scalar4) -> Self::Output { &self * &rhs }
}

impl<'a> std::ops::Mul<Scalar4> for &'a Scalar4 {
    type Output = Scalar4;
    fn mul(self, rhs: Scalar4) -> Self::Output { self * &rhs }
}

impl<'a> std::ops::Mul<&'a Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn mul(self, rhs: &Scalar4) -> Self::Output { &self * rhs }
}

impl std::ops::MulAssign<Scalar4> for Scalar4 {
    fn mul_assign(&mut self, rhs: Scalar4) {
        *self = &*self * &rhs;
    }
}

impl<'a> std::ops::MulAssign<&'a Scalar4> for Scalar4 {
    fn mul_assign(&mut self, rhs: &Scalar4) {
        *self = &*self * rhs;
    }
}

// The main implementation of the Add trait uses references, so we don't need
// to make a copy of the scalars to multiply them.
impl<'a, 'b> std::ops::Add<&'b Scalar4> for &'a Scalar4 {
    type Output = Scalar4;

    fn add(self, rhs: &Scalar4) -> Self::Output {
        match (self,rhs) {
            (Float(c), x) => Float(c + x.float_value()),
            (x, Float(c)) => Float(x.float_value() + c),
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                if pow0 != pow1 {
                    // if rt2 powers don't match, we have to fall back on float repr
                    Float(self.float_value() + rhs.float_value())
                } else {
                    let mut coeffs = *coeffs0;

                    for (i,x) in coeffs1.iter().enumerate() {
                        coeffs[i] += x;
                    }

                    Exact(*pow0, coeffs)
                }
            },
        }
    }
}

// These 3 variations take ownership of one or both args
impl std::ops::Add<Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn add(self, rhs: Scalar4) -> Self::Output { &self * &rhs }
}

impl<'a> std::ops::Add<Scalar4> for &'a Scalar4 {
    type Output = Scalar4;
    fn add(self, rhs: Scalar4) -> Self::Output { self * &rhs }
}

impl<'a> std::ops::Add<&'a Scalar4> for Scalar4 {
    type Output = Scalar4;
    fn add(self, rhs: &Scalar4) -> Self::Output { &self * rhs }
}

impl AbsDiffEq for Scalar4 {
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

// impl RelativeEq for Scalar4 {
//     fn default_max_relative() -> f64 {
//         f64::default_max_relative()
//     }

//     fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
//         let c1 = self.float_value();
//         let c2 = other.float_value();
//         f64::relative_eq(&c1.re, &c2.re, epsilon, max_relative) &&
//         f64::relative_eq(&c1.im, &c2.im, epsilon, max_relative)
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn approx_mul() {
        let s = Scalar4::real(f64::sqrt(0.3) * f64::sqrt(0.3) - 0.3);
        let t = Scalar4::zero();
        assert_ne!(s, t);
        assert_abs_diff_eq!(s, t);
    }

    #[test]
    fn sqrt_i() {
        let s = Scalar4::Exact(0, [0, 1, 0, 0]);
        assert_abs_diff_eq!(s.to_float(), Scalar4::complex(1.0 / f64::sqrt(2.0), 1.0 / f64::sqrt(2.0)));
    }

    #[test]
    fn mul_same_base() {
        let s = Scalar4::Exact(0, [1, 2, 3, 4]);
        let t = Scalar4::Exact(0, [4, 5, 6, 7]);
        assert_abs_diff_eq!((&s * &t).to_float(), s.to_float() * t.to_float());
    }

    #[test]
    fn phases() {
        assert_abs_diff_eq!(
            Scalar4::from_phase(Rational::new(4,3)) * Scalar4::from_phase(Rational::new(2,5)),
            Scalar4::from_phase(Rational::new(4,3) + Rational::new(2,5))
        );

        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(0,1)), Scalar4::one());
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(1,1)), Scalar4::real(-1.0));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(1,2)), Scalar4::complex(0.0, 1.0));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(-1,2)), Scalar4::complex(0.0, -1.0));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(1,4)), Scalar4::Exact(0, [0,1,0,0]));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(3,4)), Scalar4::Exact(0, [0,0,0,1]));
        assert_abs_diff_eq!(Scalar4::from_phase(Rational::new(7,4)), Scalar4::Exact(0, [0,0,0,-1]));
    }

    #[test]
    fn one_plus_phases() {
        assert_abs_diff_eq!(Scalar4::one_plus_phase(Rational::new(1,1)), Scalar4::zero());

        let plus = Scalar4::one_plus_phase(Rational::new(1,2));
        let minus = Scalar4::one_plus_phase(Rational::new(-1,2));
        assert_abs_diff_eq!(plus * minus, Scalar4::real(2.0));
    }
}
