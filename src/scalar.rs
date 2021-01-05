use num::complex::Complex;
use num::rational::Rational;
use num::integer;
use std::fmt;
use approx::AbsDiffEq;

pub trait Phase {
    fn mod2(&self) -> Self;
}

impl Phase for Rational {
    fn mod2(&self) -> Rational {
       Rational::new(*self.numer() % (2 * *self.denom()), *self.denom())
    }
}

/// A type for exact and approximate representation of Scalars.
/// Note that '==' is only reliable when the scalar is Exact
/// and N is a power of 2.

#[derive(Debug,Clone,PartialEq)]
pub enum Scalar {
    // An exact representation of a scalar, which is given as
    // a power of sqrt(2) and an element of Z[omega], where
    // omega is the 2N-th root of unity, represented by its
    // first N coefficients.
    Exact(i32, Vec<i32>),
    // A floating-point representation of a scalar. We should
    // fall back to this if N^2 gets too big.
    Float(Complex<f64>),
}

use Scalar::{Exact,Float};

impl Scalar {
    pub fn zero() -> Scalar {
        Exact(0, vec![0])
    }

    pub fn one() -> Scalar {
        Exact(0, vec![1])
    }

    pub fn complex(re: f64, im: f64) -> Scalar {
        Float(Complex::new(re, im))
    }

    pub fn real(re: f64) -> Scalar {
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
        (*self) *= Scalar::phase(phase);
    }

    pub fn to_float(&self) -> Scalar {
        Float(self.float_value())
    }

    pub fn phase(p: Rational) -> Scalar {
        let mut rnumer = *p.numer();
        let mut rdenom = *p.denom();

        if rdenom < 0 {
            rnumer *= -1;
            rdenom *= -1;
        }

        rnumer = rnumer.rem_euclid(2 * rdenom);
        let sgn = if rnumer >= rdenom {
            rnumer = rnumer - rdenom;
            -1
        } else {
            1
        };

        let mut coeffs: Vec<i32> = vec![0; rdenom as usize];
        coeffs[rnumer as usize] = sgn;

        Scalar::Exact(0, coeffs)
    }

    pub fn one_plus_phase(p: Rational) -> Scalar {
        let mut s = Scalar::phase(p);
        if let Scalar::Exact(_,ref mut coeffs) = s { coeffs[0] += 1; }
        s
    }

    pub fn rt2_pow(p: i32) -> Scalar {
        Scalar::Exact(p, vec![1])
    }
}

impl fmt::Display for Scalar {
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
impl<'a, 'b> std::ops::Mul<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    fn mul(self, rhs: &Scalar) -> Self::Output {
        match (self,rhs) {
            (Float(c), x) => Float(c * x.float_value()),
            (x, Float(c)) => Float(x.float_value() * c),
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                let (lcm, pad0, pad1) = if coeffs0.len() == coeffs1.len() {
                    (coeffs0.len(), 1, 1)
                } else {
                    let lcm0 = integer::lcm(coeffs0.len(), coeffs1.len());
                    (lcm0, lcm0 / coeffs0.len(), lcm0 / coeffs1.len())
                };
                let mut coeffs = vec![0; lcm];

                for (i,x) in coeffs0.iter().enumerate() {
                    for (j,y) in coeffs1.iter().enumerate() {
                        let pos = (i*pad0 + j*pad1).rem_euclid(2 * lcm);
                        if pos < lcm {
                            coeffs[pos] += x * y;
                        } else {
                            coeffs[pos - lcm] -= x * y;
                        }
                    }
                }

                Exact(pow0 + pow1, coeffs)
            },
        }
    }
}

// These 3 variations take ownership of one or both args
impl std::ops::Mul<Scalar> for Scalar {
    type Output = Scalar;
    fn mul(self, rhs: Scalar) -> Self::Output { &self * &rhs }
}

impl<'a> std::ops::Mul<Scalar> for &'a Scalar {
    type Output = Scalar;
    fn mul(self, rhs: Scalar) -> Self::Output { self * &rhs }
}

impl<'a> std::ops::Mul<&'a Scalar> for Scalar {
    type Output = Scalar;
    fn mul(self, rhs: &Scalar) -> Self::Output { &self * rhs }
}

impl std::ops::MulAssign<Scalar> for Scalar {
    fn mul_assign(&mut self, rhs: Scalar) {
        *self = &*self * &rhs;
    }
}

impl<'a> std::ops::MulAssign<&'a Scalar> for Scalar {
    fn mul_assign(&mut self, rhs: &Scalar) {
        *self = &*self * rhs;
    }
}

impl AbsDiffEq for Scalar {
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

// impl RelativeEq for Scalar {
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
        let s = Scalar::real(f64::sqrt(0.3) * f64::sqrt(0.3) - 0.3);
        let t = Scalar::zero();
        assert_ne!(s, t);
        assert_abs_diff_eq!(s, t);
    }

    #[test]
    fn sqrt_i() {
        let s = Scalar::Exact(0, vec![0, 1, 0, 0]);
        assert_abs_diff_eq!(s.to_float(), Scalar::complex(1.0 / f64::sqrt(2.0), 1.0 / f64::sqrt(2.0)));
    }

    #[test]
    fn mul_same_base() {
        let s = Scalar::Exact(0, vec![1, 2, 3]);
        let t = Scalar::Exact(0, vec![4, 5, 6]);
        assert_abs_diff_eq!((&s * &t).to_float(), s.to_float() * t.to_float());

        let s = Scalar::Exact(2, vec![1, 2, 3, 4, 5]);
        let t = Scalar::Exact(2, vec![4, 5, 6, 7, 8]);
        assert_abs_diff_eq!((&s * &t).to_float(), s.to_float() * t.to_float());
    }

    #[test]
    fn mul_diff_base() {
        let s = Scalar::Exact(-3, vec![1, 2]);
        let t = Scalar::Exact(2, vec![3, 4, 5]);
        assert_abs_diff_eq!((&s * &t).to_float(), s.to_float() * t.to_float());

        let s = Scalar::Exact(10, vec![1, 2, 3, 4, 5, 6]);
        let t = Scalar::Exact(-1, vec![7, 8]);
        assert_abs_diff_eq!((&s * &t).to_float(), s.to_float() * t.to_float());
    }

    #[test]
    fn phases() {
        assert_abs_diff_eq!(
            Scalar::phase(Rational::new(4,3)) * Scalar::phase(Rational::new(2,5)),
            Scalar::phase(Rational::new(4,3) + Rational::new(2,5))
        );

        assert_abs_diff_eq!(Scalar::phase(Rational::new(0,1)), Scalar::one());
        assert_abs_diff_eq!(Scalar::phase(Rational::new(1,1)), Scalar::real(-1.0));
        assert_abs_diff_eq!(Scalar::phase(Rational::new(1,2)), Scalar::complex(0.0, 1.0));
        assert_abs_diff_eq!(Scalar::phase(Rational::new(-1,2)), Scalar::complex(0.0, -1.0));
        assert_abs_diff_eq!(Scalar::phase(Rational::new(1,4)), Scalar::Exact(0, vec![0,1,0,0]));
        assert_abs_diff_eq!(Scalar::phase(Rational::new(3,4)), Scalar::Exact(0, vec![0,0,0,1]));
        assert_abs_diff_eq!(Scalar::phase(Rational::new(7,4)), Scalar::Exact(0, vec![0,0,0,-1]));
    }

    #[test]
    fn one_plus_phases() {
        assert_abs_diff_eq!(Scalar::one_plus_phase(Rational::new(1,1)), Scalar::zero());

        let plus = Scalar::one_plus_phase(Rational::new(1,2));
        let minus = Scalar::one_plus_phase(Rational::new(-1,2));
        assert_abs_diff_eq!(plus * minus, Scalar::real(2.0));
    }

    #[test]
    fn redundant_rep() {
        // If order is not a power of 2, there is some redudancy in the Exact representation.
        // This could be fixed by implementing a function that reduces to the Zumbroich basis.
        // See e.g. https://github.com/CyclotomicFields/cyclotomic/blob/master/src/fields/dense/basis.rs
        let s = Scalar::Exact(0, vec![1,-1,1]);
        assert_ne!(s, Scalar::zero());
        assert_abs_diff_eq!(s, Scalar::zero());
    }
}

