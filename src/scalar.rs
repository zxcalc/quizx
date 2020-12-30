use num::complex::Complex;
use num::rational::Rational;
use num::integer;
use std::fmt;

/// A type for exact and approximate representation of Scalars.
/// Note that '==' is only reliable when the scalar is Exact
/// and N is a power of 2.

#[derive(Clone,PartialEq)]
pub enum Scalar {
    // An exact representation of a scalar, which is given as
    // a power of sqrt(2) and an element of Z[omega], where
    // omega is the 2N-th root of unity, represented by its
    // first N coefficients.
    Exact(isize, Vec<isize>),
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

    // fn raise_order(&mut self, k: usize) {
    //     match self {
    //         Exact(_, coeffs) => {
    //             let mut new_coeffs = vec![0;k * coeffs.len()];
    //             for (i,x) in coeffs.iter().enumerate() {
    //                 new_coeffs[k*i] = *x;
    //             }
    //             *coeffs = new_coeffs;
    //         }
    //         _ => {}
    //     }
    // }

    pub fn as_float(&self) -> Complex<f64> {
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

    pub fn str(&self) -> String {
        match self {
            Exact(pow,coeffs) => {
                let mut s = format!("rt(2)^{} (", pow);

                for (i,c) in coeffs.iter().enumerate() {
                    if i == 0 {
                        s.push_str(&c.to_string());
                    } else {
                        s.push_str(&format!(" + {} * om^{}", c, i));
                    }
                }

                s.push_str(")");

                s
            },
            Float(c) => format!("{}", c),
        }
    }

    pub fn mult_with(&self, rhs: &Scalar) -> Scalar {
        match (self,rhs) {
            (Float(c), x) => Float(c * x.as_float()),
            (x, Float(c)) => Float(x.as_float() * c),
            (Exact(pow0, coeffs0), Exact(pow1, coeffs1)) => {
                let lcm = integer::lcm(coeffs0.len(), coeffs1.len());
                let (pad0, pad1) = (lcm / coeffs0.len(), lcm / coeffs1.len());
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

    pub fn phase(p: Rational) -> Scalar {
        let mut rnumer = *p.numer();
        let mut rdenom = *p.denom();

        if rdenom < 0 {
            rnumer *= -1;
            rdenom *= -1;
        }

        rnumer = rnumer.rem_euclid(2 * rdenom);
        let sgn = if rnumer > rdenom {
            rnumer = rnumer - rdenom;
            -1
        } else {
            1
        };

        let mut coeffs: Vec<isize> = vec![0; rdenom as usize];
        coeffs[rnumer as usize] = sgn;

        Scalar::Exact(0, coeffs)
    }

    pub fn one_plus_phase(p: Rational) -> Scalar {
        let mut s = Scalar::phase(p);
        if let Scalar::Exact(_,ref mut coeffs) = s { coeffs[0] = 1; }
        s
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.str())
    }
}

impl<'a, 'b> std::ops::Mul<&'a Scalar> for &'b Scalar {
    type Output = Scalar;

    fn mul(self, rhs: &Scalar) -> Self::Output {
        self.mult_with(rhs)
    }
}
