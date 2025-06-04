use approx::AbsDiffEq;
use num::complex::Complex;
pub use num::traits::identities::{One, Zero};
use num::{Float, Rational64, ToPrimitive};
use std::f64::consts::{PI, SQRT_2};
use std::fmt;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use crate::phase::Phase;
pub use crate::scalar_traits::{FromPhase, Sqrt2};
/// This is the main representation of scalars used in QuiZX. It is a wrapper around
/// four `f64` values, used to represent the coefficients in a complex number of
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
/// These correspond exactly to the signficand and exponent of a floating point
/// number (see e.g. [floats on Wikipedia](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)),
/// so as long as x and y for each coefficient are not too large to fit in an `f64`,
/// `FScalar` can exactly represent any Clifford+T scalar.
///
/// Note that ω² = i, so for all other complex numbers, `FScalar` gives an approximate
/// representation of that number as a + c ω² = a + i c. This allows easy conversions
/// to and from `Complex<f64>` provided in the standard library.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FScalar {
    c: [f64; 4],
}

impl FScalar {
    /// Constructs a scalar from a power of two and 4 integer coefficients. For
    /// `coeffs := [a,b,c,d]`, the resulting scalar is: `2^pow * (a + b ω + c ω² + d ω³)`.
    pub fn dyadic(pow: i32, coeffs: [i64; 4]) -> Self {
        let f = 2.0_f64.powi(pow);
        FScalar {
            c: coeffs.map(|coeff| f * (coeff as f64)),
        }
    }

    /// Constructs a scalar as a real number with the given float value
    pub fn real(r: f64) -> Self {
        FScalar {
            c: [r, 0.0, 0.0, 0.0],
        }
    }

    /// Constructs a scalar as a complex number `re + im * i`.
    pub fn complex(re: f64, im: f64) -> Self {
        FScalar {
            c: [re, 0.0, im, 0.0],
        }
    }

    /// Constructs a scalar `1 + e^(i π * phase)` for the given `Phase`
    pub fn one_plus_phase(phase: impl Into<Phase>) -> Self {
        FScalar::one() + FScalar::from_phase(phase)
    }

    /// Convience method for multipling by the given power of sqrt(2)
    pub fn mul_sqrt2_pow(&mut self, p: i32) {
        *self *= FScalar::sqrt2_pow(p);
    }

    /// Convience method for multipling by `e^(i π * phase)` for the given `Phase`
    pub fn mul_phase(&mut self, phase: impl Into<Phase>) {
        *self *= FScalar::from_phase(phase);
    }

    /// Convience method for multipling by `1 + e^(i π * phase)` for the given `Phase`
    pub fn mul_one_plus_phase(&mut self, phase: impl Into<Phase>) {
        *self *= FScalar::one_plus_phase(phase)
    }

    /// Returns the complex conjugate of the scalar
    pub fn conj(&self) -> Self {
        FScalar {
            c: [self.c[0], -self.c[3], -self.c[2], -self.c[1]],
        }
    }

    /// Converts `FScalar` into a `Complex<f64>`
    pub fn complex_value(&self) -> Complex<f64> {
        self.into()
    }

    /// Returns an array of 4 pairs giving each coefficient as a pair of integers
    /// `(a,p)`, meaning `a * 2^p`. This is convenient e.g. for the base-2 scientific
    /// notation used to output scalars.
    pub fn exact_dyadic_form(&self) -> [(i64, i16); 4] {
        self.c.map(|f| {
            let (mut m, mut e, s) = f.integer_decode();
            if m == 0 {
                (0, 0)
            } else {
                while m != 0 && m % 2 == 0 {
                    m /= 2;
                    e += 1;
                }

                ((m as i64) * (s as i64), e)
            }
        })
    }

    /// The number of non-zero coefficients
    fn num_coeffs(&self) -> usize {
        self.c
            .iter()
            .fold(0, |n, &co| if co != 0.0 { n + 1 } else { n })
    }

    /// If the scalar is an exact representation of a number in the form `exp(i k π / 4) * sqrt(2)^p`,
    /// return `Some((k/4, p))`, otherwise return `None`.
    pub fn exact_phase_and_sqrt2_pow(&self) -> Option<(Phase, i16)> {
        let mut p = 0;
        let mut s = *self;
        if self.num_coeffs() > 1 {
            p = -1;
            s *= Self::sqrt2()
        };

        let edf = s.exact_dyadic_form();

        for i in 0..4 {
            if edf[i].0 == 1
                && edf[(i + 1) % 4].0 == 0
                && edf[(i + 2) % 4].0 == 0
                && edf[(i + 3) % 4].0 == 0
            {
                return Some((Phase::new(Rational64::new(i as i64, 4)), edf[i].1 * 2 + p));
            } else if edf[i].0 == -1
                && edf[(i + 1) % 4].0 == 0
                && edf[(i + 2) % 4].0 == 0
                && edf[(i + 3) % 4].0 == 0
            {
                return Some((
                    Phase::new(Rational64::new((i + 4) as i64, 4)),
                    edf[i].1 * 2 + p,
                ));
            }
        }

        None
    }
}

impl Default for FScalar {
    fn default() -> Self {
        Self::zero()
    }
}

impl Sqrt2 for FScalar {
    fn sqrt2_pow(p: i32) -> Self {
        FScalar {
            c: if p % 2 == 0 {
                [2.0_f64.powi(p / 2), 0.0, 0.0, 0.0]
            } else {
                let f = 2.0_f64.powi((p - 1) / 2);
                [0.0, f, 0.0, -f]
            },
        }
    }
}

impl FromPhase for FScalar {
    fn from_phase(phase: impl Into<Phase>) -> Self {
        let p: Phase = phase.into();
        p.into()
    }

    fn minus_one() -> Self {
        FScalar {
            c: [-1.0, 0.0, 0.0, 0.0],
        }
    }
}

impl ndarray::ScalarOperand for FScalar {}

impl fmt::Display for FScalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut fst = true;

        let dy = self.exact_dyadic_form();
        for (i, (mut v, mut e)) in dy.iter().enumerate() {
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

impl Add<&FScalar> for &FScalar {
    type Output = FScalar;
    fn add(self, rhs: &FScalar) -> Self::Output {
        FScalar {
            c: [
                self.c[0] + rhs.c[0],
                self.c[1] + rhs.c[1],
                self.c[2] + rhs.c[2],
                self.c[3] + rhs.c[3],
            ],
        }
    }
}

impl Sum<FScalar> for FScalar {
    fn sum<I: Iterator<Item = FScalar>>(iter: I) -> Self {
        iter.fold(FScalar::zero(), |accumulator, current_fscalar| {
            accumulator + current_fscalar // Uses FScalar + FScalar
        })
    }
}

// These 3 variations take ownership of one or both args
impl Add<FScalar> for FScalar {
    type Output = FScalar;
    fn add(self, rhs: FScalar) -> Self::Output {
        let rself = &self;
        let rrhs = &rhs;
        rself + rrhs
    }
}

impl Add<FScalar> for &FScalar {
    type Output = FScalar;
    fn add(self, rhs: FScalar) -> Self::Output {
        let rrhs = &rhs;
        self + rrhs
    }
}

impl Add<&FScalar> for FScalar {
    type Output = FScalar;
    fn add(self, rhs: &FScalar) -> Self::Output {
        let rself = &self;
        rself + rhs
    }
}

impl AddAssign for FScalar {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl AddAssign<&FScalar> for FScalar {
    fn add_assign(&mut self, rhs: &FScalar) {
        *self = *self + rhs;
    }
}

impl Sub<&FScalar> for &FScalar {
    type Output = FScalar;
    fn sub(self, rhs: &FScalar) -> Self::Output {
        FScalar {
            c: [
                self.c[0] - rhs.c[0],
                self.c[1] - rhs.c[1],
                self.c[2] - rhs.c[2],
                self.c[3] - rhs.c[3],
            ],
        }
    }
}

impl Sub<FScalar> for FScalar {
    type Output = FScalar;
    fn sub(self, rhs: FScalar) -> Self::Output {
        &self - &rhs
    }
}

impl Sub<FScalar> for &FScalar {
    type Output = FScalar;
    fn sub(self, rhs: FScalar) -> Self::Output {
        let rrhs = &rhs;
        self - rrhs
    }
}

impl Sub<&FScalar> for FScalar {
    type Output = FScalar;
    fn sub(self, rhs: &FScalar) -> Self::Output {
        let rself = &self;
        rself - rhs
    }
}

impl SubAssign for FScalar {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl SubAssign<&FScalar> for FScalar {
    fn sub_assign(&mut self, rhs: &FScalar) {
        *self = *self - rhs;
    }
}

impl Mul<&FScalar> for &FScalar {
    type Output = FScalar;
    fn mul(self, rhs: &FScalar) -> Self::Output {
        let mut c = [0.0, 0.0, 0.0, 0.0];
        for i in 0..4 {
            if self.c[i] != 0.0 {
                for j in 0..4 {
                    let pos = (i + j) % 8;
                    if pos < 4 {
                        c[pos] += self.c[i] * rhs.c[j];
                    } else {
                        c[pos - 4] += -self.c[i] * rhs.c[j];
                    }
                }
            }
        }
        FScalar { c }
    }
}

// These 3 variations take ownership of one or both args
impl Mul<FScalar> for FScalar {
    type Output = FScalar;
    fn mul(self, rhs: FScalar) -> Self::Output {
        let rself = &self;
        let rrhs = &rhs;
        rself * rrhs
    }
}
impl Mul<FScalar> for &FScalar {
    type Output = FScalar;
    fn mul(self, rhs: FScalar) -> Self::Output {
        let rrhs = &rhs;
        self * rrhs
    }
}
impl Mul<&FScalar> for FScalar {
    type Output = FScalar;
    fn mul(self, rhs: &FScalar) -> Self::Output {
        let rself = &self;
        rself * rhs
    }
}

/// Implements *=
impl MulAssign<FScalar> for FScalar {
    fn mul_assign(&mut self, rhs: FScalar) {
        *self = *self * rhs;
    }
}

// Variation takes ownership of rhs
impl MulAssign<&FScalar> for FScalar {
    fn mul_assign(&mut self, rhs: &FScalar) {
        *self = *self * rhs;
    }
}

impl Product<FScalar> for FScalar {
    fn product<I: Iterator<Item = FScalar>>(iter: I) -> Self {
        iter.fold(FScalar::one(), |accumulator, current_fscalar| {
            accumulator * current_fscalar // Uses FScalar * FScalar
        })
    }
}

impl Zero for FScalar {
    fn zero() -> Self {
        FScalar {
            c: [0.0, 0.0, 0.0, 0.0],
        }
    }

    fn is_zero(&self) -> bool {
        self.c == [0.0, 0.0, 0.0, 0.0]
    }
}

impl One for FScalar {
    fn one() -> Self {
        FScalar {
            c: [1.0, 0.0, 0.0, 0.0],
        }
    }

    fn is_one(&self) -> bool {
        self.c == [1.0, 0.0, 0.0, 0.0]
    }
}

impl AbsDiffEq<FScalar> for FScalar {
    type Epsilon = <f64 as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        1e-10f64
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let c1: Complex<f64> = self.into();
        let c2: Complex<f64> = other.into();
        f64::abs_diff_eq(&c1.re, &c2.re, epsilon) && f64::abs_diff_eq(&c1.im, &c2.im, epsilon)
    }
}

impl From<f64> for FScalar {
    fn from(value: f64) -> Self {
        FScalar {
            c: [value, 0.0, 0.0, 0.0],
        }
    }
}

impl From<i32> for FScalar {
    fn from(value: i32) -> Self {
        FScalar {
            c: [value as f64, 0.0, 0.0, 0.0],
        }
    }
}

impl From<Complex<f64>> for FScalar {
    fn from(value: Complex<f64>) -> Self {
        FScalar {
            c: [value.re, 0.0, value.im, 0.0],
        }
    }
}

impl From<[i32; 4]> for FScalar {
    fn from(value: [i32; 4]) -> Self {
        FScalar {
            c: [
                value[0] as f64,
                value[1] as f64,
                value[2] as f64,
                value[3] as f64,
            ],
        }
    }
}

impl From<[f64; 4]> for FScalar {
    fn from(value: [f64; 4]) -> Self {
        FScalar { c: value }
    }
}

impl From<&FScalar> for Complex<f64> {
    fn from(value: &FScalar) -> Self {
        Complex {
            re: value.c[0] + (value.c[1] - value.c[3]) * 0.5 * SQRT_2,
            im: value.c[2] + (value.c[1] + value.c[3]) * 0.5 * SQRT_2,
        }
    }
}

impl From<FScalar> for Complex<f64> {
    fn from(value: FScalar) -> Self {
        Complex {
            re: value.c[0] + (value.c[1] - value.c[3]) * 0.5 * SQRT_2,
            im: value.c[2] + (value.c[1] + value.c[3]) * 0.5 * SQRT_2,
        }
    }
}

impl From<Phase> for FScalar {
    fn from(value: Phase) -> Self {
        let r: Rational64 = value.into();
        if 4 % r.denom() == 0 {
            let pos = (r.numer() * (4 / r.denom())).rem_euclid(8) as usize;
            let mut c: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            if pos >= 4 {
                c[pos - 4] = -1.0;
            } else {
                c[pos] = 1.0;
            }
            FScalar { c }
        } else {
            let angle = PI * r.to_f64().unwrap();
            FScalar {
                c: [f64::cos(angle), 0.0, f64::sin(angle), 0.0],
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn display() {
        let s = FScalar::zero();
        assert_eq!(format!("{s}"), "0");

        let s: FScalar = [1, 2, 3, 4].into();
        assert_eq!(format!("{s}"), "1 + 2 ω + 3 ω² + 4 ω³");

        let s: FScalar = [-1, -2, -3, -4].into();
        assert_eq!(format!("{s}"), "-1 - 2 ω - 3 ω² - 4 ω³");

        let s: FScalar = [0, 2, 0, 4].into();
        assert_eq!(format!("{s}"), "2 ω + 4 ω³");

        let s: FScalar = [0.5, 0.25, 0.125, 0.0625].into();
        assert_eq!(format!("{s}"), "1e-1 + 1e-2 ω + 1e-3 ω² + 1e-4 ω³");

        let s: FScalar = [
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
        let s4: FScalar = 4.into();
        let s10: FScalar = 10.into();
        let s14: FScalar = 14.into();
        let sm1: FScalar = (-1).into();
        let s40: FScalar = 40.into();
        let sm14: FScalar = (-14).into();
        let sm40: FScalar = (-40).into();

        assert_eq!(s4 + s10, s14);
        assert_eq!(s4 * s10, s40);
        assert_eq!(sm1 * sm1, FScalar::one());
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
        let sa: FScalar = a.into();
        let sb: FScalar = b.into();
        let sc: FScalar = c.into();
        let sd: FScalar = d.into();

        assert_abs_diff_eq!(sa * sa, (a * a).into());
        assert_abs_diff_eq!(sa * sb, (a * b).into());
        assert_abs_diff_eq!(sc + (sa * sb), (c + (a * b)).into());
        assert_abs_diff_eq!(sd - (sa * sb), (d - (a * b)).into());
    }

    #[test]
    fn complex_arith() {
        let one: FScalar = 1.into();
        let i: FScalar = [0, 0, 1, 0].into();
        let om: FScalar = [0, 1, 0, 0].into();
        let sqrt2 = FScalar::sqrt2_pow(1);
        assert_eq!(om * om, i);
        assert_eq!((one + i) * (one + i).conj(), one + one);
        assert_eq!(om + om.conj(), sqrt2);

        let c1: Complex<f64> = (one + i + i).into();
        let c2: Complex<f64> = Complex::new(1.0, 2.0);
        assert_eq!(c1, c2);

        let c1: Complex<f64> = (FScalar::sqrt2_pow(3) + i * sqrt2).into();
        let c2: Complex<f64> = Complex::new(2.0 * SQRT_2, SQRT_2);
        assert_abs_diff_eq!(c1.re, c2.re);
        assert_abs_diff_eq!(c1.im, c2.im);
    }

    #[test]
    fn sqrt2() {
        // n.b. exact equality in the sqrt2_pow tests. This may break if conversion to Complex is
        // implemented differently.
        let sqrt2 = FScalar::sqrt2_pow(1);
        let sqrt2_c: Complex<f64> = sqrt2.into();
        assert_eq!(sqrt2_c.re, SQRT_2);
        assert_eq!(sqrt2_c.im, 0.0);

        let sqrt2_pow = FScalar::sqrt2_pow(7);
        let sqrt2_pow_c: Complex<f64> = sqrt2_pow.into();
        assert_eq!(sqrt2_pow_c.re, 2.0 * 2.0 * 2.0 * SQRT_2);
        assert_eq!(sqrt2_pow_c.im, 0.0);

        let two: FScalar = 2.into();
        let a = FScalar::sqrt2_pow(10);
        let b = FScalar::sqrt2_pow(11);
        let c: FScalar = 32.into();
        assert_eq!(two, sqrt2 * sqrt2);
        assert_eq!(a * sqrt2, b);
        assert_eq!(a, c);
    }

    #[test]
    fn from_t_phase() {
        let s: FScalar = Phase::new(Rational64::new(0, 1)).into();
        assert_eq!(
            s,
            FScalar {
                c: [1.0, 0.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(1, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 1.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(1, 2)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, 1.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(3, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, 0.0, 1.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(1, 1)).into();
        assert_eq!(
            s,
            FScalar {
                c: [-1.0, 0.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(5, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, -1.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(3, 2)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, -1.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(7, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, 0.0, -1.0]
            }
        );

        let s: FScalar = Phase::new(Rational64::new(0, 1)).into();
        assert_eq!(
            s,
            FScalar {
                c: [1.0, 0.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(-7, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 1.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(-3, 2)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, 1.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(-5, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, 0.0, 1.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(-1, 1)).into();
        assert_eq!(
            s,
            FScalar {
                c: [-1.0, 0.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(-3, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, -1.0, 0.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(-1, 2)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, -1.0, 0.0]
            }
        );
        let s: FScalar = Phase::new(Rational64::new(-1, 4)).into();
        assert_eq!(
            s,
            FScalar {
                c: [0.0, 0.0, 0.0, -1.0]
            }
        );
    }

    #[test]
    fn from_gen_phase() {
        let r = Rational64::new(-5, 37);
        let s1: FScalar = Phase::new(r).into();
        let c = Complex::new(0.0, PI * r.to_f64().unwrap()).exp();
        let s2: FScalar = c.into();

        assert_eq!(s1, s2);
        let r = Rational64::new(12, 117);
        let s1: FScalar = Phase::new(r).into();
        let c = Complex::new(0.0, PI * r.to_f64().unwrap()).exp();
        let s2: FScalar = c.into();
        assert_eq!(s1, s2);
    }
}
