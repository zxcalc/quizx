use num::{integer,Integer};
use num::complex::Complex;
use num::rational::Rational;
pub use num::traits::identities::{Zero,One};
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

pub trait FromPhase {
    fn from_phase(p: Rational) -> Self;
}

pub trait Sqrt2 {
    fn sqrt2() -> Self;
    fn one_over_sqrt2() -> Self;
}


/// A list of coefficients. We give this as a parameter to allow either
/// fixed-size lists (e.g. [i32;4]) or dynamic ones (e.g. Vec<i32>). Only
/// the former can be used in tensors and matrices, because they have to
/// implement Copy (i.e. size must be known at compile time).
pub trait Coeffs: PartialEq + Clone + std::ops::IndexMut<usize,Output=i32> {
    fn len(&self) -> usize;
    fn zero() -> Self;
    fn one() -> Self;
    fn new(sz: usize) -> Option<(Self,usize)>;
}

/// A type for exact and approximate representation of Scalars.
/// Note that '==' is only reliable when the scalar is Exact
/// and N is a power of 2.

#[derive(Debug,Clone,PartialEq)]
pub enum Scalar<T: Coeffs> {
    // An exact representation of a scalar, which is given as
    // a power of sqrt(2) and an element of Z[omega], where
    // omega is the 2N-th root of unity, represented by its
    // first N coefficients.
    Exact(i32, T),
    // A floating-point representation of a scalar. We should
    // fall back to this if N^2 gets too big.
    Float(Complex<f64>),
}

/// Implement Copy whenever our coefficient list allows us to.
impl<T: Coeffs + Copy> Copy for Scalar<T> {}

use Scalar::{Exact,Float};

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
        (*self) *= Scalar::from_phase(phase);
    }

    pub fn to_float(&self) -> Scalar<T> {
        Float(self.float_value())
    }

    pub fn one_plus_phase(p: Rational) -> Scalar<T> {
        Scalar::one() + Scalar::from_phase(p)
    }

    pub fn rt2_pow(p: i32) -> Scalar<T> {
        Scalar::Exact(p, T::one())
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
    fn sqrt2() -> Scalar<T> { Scalar::rt2_pow(1) }
    fn one_over_sqrt2() -> Scalar<T> { Scalar::rt2_pow(-1) }
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
                let sgn = if rnumer >= rdenom {
                    rnumer = rnumer - rdenom;
                    -1
                } else {
                    1
                };
                let i = (rnumer * pad as isize).rem_euclid(2 * rdenom);
                coeffs[i as usize] = sgn;
                Scalar::Exact(0, coeffs)
            },
            None => {
                // TODO
                Scalar::Float(Complex::zero())
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
                let (lcm, pad0, pad1) = if coeffs0.len() == coeffs1.len() {
                    (coeffs0.len(), 1, 1)
                } else {
                    let lcm0 = integer::lcm(coeffs0.len(), coeffs1.len());
                    (lcm0, lcm0 / coeffs0.len(), lcm0 / coeffs1.len())
                };

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

                        Exact(pow0 + pow1, coeffs)
                    },
                    None => {
                        // TODO
                        Float(Complex::zero())
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
                if pow0 != pow1 {
                    // if rt2 powers don't match, we have to fall back on float repr
                    Float(self.float_value() + rhs.float_value())
                } else {
                    let (lcm, pad0, pad1) = if coeffs0.len() == coeffs1.len() {
                        (coeffs0.len(), 1, 1)
                    } else {
                        let lcm0 = integer::lcm(coeffs0.len(), coeffs1.len());
                        (lcm0, lcm0 / coeffs0.len(), lcm0 / coeffs1.len())
                    };

                    match T::new(lcm) {
                        Some((mut coeffs, pad)) => {
                            for i in 0..coeffs0.len() {
                                coeffs[i*pad*pad0] += coeffs0[i];
                            }

                            for i in 0..coeffs1.len() {
                                coeffs[i*pad*pad1] += coeffs1[i];
                            }

                            Exact(*pow0, coeffs)
                        },
                        None => {
                            // TODO
                            Float(Complex::zero())
                        }
                    }

                }
            },
        }
    }
}

// These 3 variations take ownership of one or both args
impl<T: Coeffs> std::ops::Add<Scalar<T>> for Scalar<T> {
    type Output = Scalar<T>;
    fn add(self, rhs: Scalar<T>) -> Self::Output { &self * &rhs }
}

impl<'a, T: Coeffs> std::ops::Add<Scalar<T>> for &'a Scalar<T> {
    type Output = Scalar<T>;
    fn add(self, rhs: Scalar<T>) -> Self::Output { self * &rhs }
}

impl<'a, T: Coeffs> std::ops::Add<&'a Scalar<T>> for Scalar<T> {
    type Output = Scalar<T>;
    fn add(self, rhs: &Scalar<T>) -> Self::Output { &self * rhs }
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

/// Implements Coeffs for an array of fixed size $n, and defines
/// the associated scalar type.
macro_rules! fixed_size_scalar {
    ( $name:ident, $n:expr ) => {
        impl Coeffs for [i32;$n] {
            fn len(&self) -> usize { $n }
            fn zero() -> [i32;$n] { [0;$n] }
            fn one() -> [i32;$n] { let mut a = [0;$n]; a[0] = 1; a }
            fn new(sz: usize) -> Option<([i32;$n],usize)> {
                if (sz as i32).divides(&$n) {
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

        // assert_abs_diff_eq!(Scalar::from_phase(Rational::new(0,1)), Scalar::one());
        // assert_abs_diff_eq!(Scalar::from_phase(Rational::new(1,1)), Scalar::real(-1.0));
        // assert_abs_diff_eq!(Scalar::from_phase(Rational::new(1,2)), Scalar::complex(0.0, 1.0));
        // assert_abs_diff_eq!(Scalar::from_phase(Rational::new(-1,2)), Scalar::complex(0.0, -1.0));
        assert_abs_diff_eq!(Scalar::from_phase(Rational::new(1,4)), Scalar::Exact(0, [0,1,0,0]));
        assert_abs_diff_eq!(Scalar::from_phase(Rational::new(3,4)), Scalar::Exact(0, [0,0,0,1]));
        assert_abs_diff_eq!(Scalar::from_phase(Rational::new(7,4)), Scalar::Exact(0, [0,0,0,-1]));
    }

    // #[test]
    // fn one_plus_phases() {
    //     assert_abs_diff_eq!(Scalar::one_plus_phase(Rational::new(1,1)), Scalar::zero());

    //     let plus = Scalar::one_plus_phase(Rational::new(1,2));
    //     let minus = Scalar::one_plus_phase(Rational::new(-1,2));
    //     assert_abs_diff_eq!(plus * minus, Scalar::real(2.0));
    // }
}
