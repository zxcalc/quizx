use approx::AbsDiffEq;
use derive_more::derive::{Display, Error};
use num::{Float, Zero};
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

// use type aliases to make it easy to change the precision
pub type Mantissa = u64;
pub type SignedMantissa = i64;
pub type DoubleMantissa = u128;
pub type Exponent = i32;

#[derive(Default, Debug, Display, Error, PartialEq, Eq)]
#[display("exponent is too small or large for destination type")]
pub struct DyadicExponentOverflowError;

const SIGN: u8 = 0x01;
const SIGN_OFF: u8 = 0xff ^ SIGN;
const APPROX: u8 = 0x02;
const APPROX_OFF: u8 = 0xff ^ APPROX;

// A dyadic is essentially a floating point number, but we have more control over precision
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Dyadic {
    flags: u8,
    exp: Exponent,
    val: Mantissa,
}

impl Dyadic {
    pub fn new(val: SignedMantissa, exp: Exponent) -> Self {
        let mut d = if val.is_negative() {
            Dyadic {
                flags: SIGN,
                val: -val as Mantissa,
                exp,
            }
        } else {
            Dyadic {
                flags: 0,
                val: val as Mantissa,
                exp,
            }
        };
        d.normalize();
        d
    }

    #[inline]
    pub fn sign(&self) -> bool {
        self.flags & SIGN == SIGN
    }

    #[inline]
    pub fn approx(&self) -> bool {
        self.flags & APPROX == APPROX
    }

    #[inline]
    pub fn set_approx(&mut self, approx: bool) {
        if approx {
            self.flags |= APPROX;
        } else {
            self.flags &= APPROX_OFF;
        }
    }

    #[inline]
    pub fn val_and_exp(&self) -> (SignedMantissa, Exponent) {
        if self.is_zero() {
            (0, 0)
        } else {
            let shift = self.val.trailing_zeros();
            let v = self.val.wrapping_shr(shift) as SignedMantissa;
            if self.sign() {
                (-v, self.exp + (shift as Exponent))
            } else {
                (v, self.exp + (shift as Exponent))
            }
        }
    }

    #[inline]
    pub fn val(&self) -> SignedMantissa {
        let shift = self.val.trailing_zeros();
        let v = self.val.wrapping_shr(shift) as SignedMantissa;
        if self.sign() {
            -v
        } else {
            v
        }
    }

    #[inline]
    pub fn exp(&self) -> Exponent {
        if self.is_zero() {
            0
        } else {
            self.exp + (self.val.trailing_zeros() as Exponent)
        }
    }

    #[inline]
    fn normalize(&mut self) {
        if self.val == 0 {
            self.exp = 0;
            self.flags &= SIGN_OFF;
        } else {
            let head = self.val.leading_zeros();
            self.exp -= head as Exponent;
            self.val = self.val.wrapping_shl(head);
        }
    }

    #[inline]
    pub fn abs(&self) -> Self {
        let mut d = *self;
        d.flags &= SIGN_OFF;
        d
    }
}

impl fmt::Debug for Dyadic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)?;
        if self.approx() {
            write!(f, "~")?;
        }
        Ok(())
    }
}

impl fmt::Display for Dyadic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (mut v, mut e) = self.val_and_exp();

        let bnd: SignedMantissa = 1024;
        if v > -bnd && v < bnd && e.is_positive() && e < 10 {
            v *= (2 as SignedMantissa).pow(e as u32);
            e = 0;
        }

        write!(f, "{}", v)?;
        if e != 0 {
            write!(f, "e{}", e)?;
        }

        Ok(())
    }
}

impl PartialOrd for Dyadic {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Dyadic {
    fn cmp(&self, other: &Self) -> Ordering {
        // first compute the order as if "self" was positive
        let ord = if self.sign() != other.sign() {
            Ordering::Greater
        } else if self.exp == other.exp {
            self.val.cmp(&other.val)
        } else {
            self.exp.cmp(&other.exp)
        };

        // if "self" is actually negative, flip the order
        if self.sign() {
            ord.reverse()
        } else {
            ord
        }
    }
}

impl AbsDiffEq for Dyadic {
    type Epsilon = Dyadic;
    fn default_epsilon() -> Self::Epsilon {
        Dyadic::new(1, -100)
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        (*self - *other).abs() < epsilon
    }
}

impl Neg for Dyadic {
    type Output = Dyadic;

    #[inline]
    fn neg(mut self) -> Self::Output {
        if self.val != 0 {
            self.flags ^= SIGN;
        }
        self
    }
}

impl Add for Dyadic {
    type Output = Dyadic;

    fn add(mut self, mut rhs: Dyadic) -> Self::Output {
        if self.val == 0 {
            return rhs;
        } else if rhs.val == 0 {
            return self;
        }

        let shift;
        if self.exp > rhs.exp {
            shift = (self.exp - rhs.exp) as u32;
            if rhs.val.trailing_zeros() < shift {
                self.flags |= APPROX;
            }
            rhs.val = if shift >= (size_of::<Mantissa>() as u32) * 8 {
                0
            } else {
                rhs.val.wrapping_shr(shift)
            };
        } else if rhs.exp > self.exp {
            shift = (rhs.exp - self.exp) as u32;
            if self.val.trailing_zeros() < shift {
                self.flags |= APPROX;
            }
            self.val = if shift >= (size_of::<Mantissa>() as u32) * 8 {
                0
            } else {
                self.val.wrapping_shr(shift)
            };
            self.exp = rhs.exp;
        } else {
            shift = 0;
        }

        if self.sign() != rhs.sign() {
            if self.val > rhs.val {
                self.val -= rhs.val;
                self.flags |= rhs.flags & APPROX;
                self.normalize();
            } else if self.val < rhs.val {
                self.val = rhs.val - self.val;
                self.flags = rhs.flags | (self.flags & APPROX);
                self.normalize();
            } else {
                self.val = 0;
                self.exp = 0;
                self.flags = (self.flags & APPROX) | (rhs.flags & APPROX);
            }
        } else {
            let overflow = shift == 0
                || match self.val.checked_add(rhs.val) {
                    Some(v) => {
                        self.val = v;
                        false
                    }
                    None => true,
                };

            if overflow {
                if (self.val & 1 == 1) || (rhs.val & 1 == 1) {
                    self.flags |= APPROX;
                }

                self.val = self.val.wrapping_shr(1);
                rhs.val = rhs.val.wrapping_shr(1);
                self.val += rhs.val;
                self.exp += 1;
            }

            self.flags |= rhs.flags & APPROX;
        };

        self
    }
}

impl AddAssign for Dyadic {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl Sub for Dyadic {
    type Output = Dyadic;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl SubAssign for Dyadic {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl Mul for Dyadic {
    type Output = Dyadic;

    fn mul(mut self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            return self;
        } else if rhs.is_zero() {
            return rhs;
        }

        self.exp += rhs.exp;
        self.flags |= rhs.flags & APPROX;
        self.flags ^= rhs.flags & SIGN;

        let mut v = (self.val as DoubleMantissa) * (rhs.val as DoubleMantissa);
        let lead = v.leading_zeros();
        if lead < (size_of::<Mantissa>() as u32) * 8 {
            let shift = (size_of::<Mantissa>() as u32) * 8 - lead;
            if v.trailing_zeros() < shift {
                self.flags |= APPROX;
            }
            v = v.wrapping_shr(shift);
            self.exp += shift as Exponent;
        }

        self.val = v as Mantissa;
        self
    }
}

impl MulAssign for Dyadic {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl Zero for Dyadic {
    fn is_zero(&self) -> bool {
        self.val == 0
    }

    fn zero() -> Self {
        Dyadic {
            flags: 0,
            exp: 0,
            val: 0,
        }
    }
}

impl From<f64> for Dyadic {
    fn from(value: f64) -> Self {
        let (m, e, s) = value.integer_decode();
        let mut d = Dyadic {
            flags: if s == -1 { SIGN } else { 0 } | APPROX,
            exp: e as Exponent,
            val: m as Mantissa,
        };
        d.normalize();
        d
    }
}

impl TryFrom<Dyadic> for f64 {
    type Error = DyadicExponentOverflowError;
    fn try_from(value: Dyadic) -> Result<Self, Self::Error> {
        if value.exp >= f64::MIN_EXP && value.exp <= f64::MAX_EXP {
            let (v, e) = value.val_and_exp();
            Ok((v as f64) * 2.0f64.powi(e))
        } else {
            Err(DyadicExponentOverflowError)
        }
    }
}

impl From<SignedMantissa> for Dyadic {
    fn from(value: SignedMantissa) -> Self {
        Dyadic::new(value, 0)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use rstest::rstest;

    #[rstest]
    #[case(Dyadic::new(24, 0), "24")]
    #[case(Dyadic::new(-100, 0), "-100")]
    #[case(Dyadic::new(15, 100), "15e100")]
    #[case(Dyadic::new(2048, 0), "1e11")]
    #[case(Dyadic::new(-2048, 0), "-1e11")]
    #[case(Dyadic::new(1, -11), "1e-11")]
    #[case(Dyadic::new(-1, -11), "-1e-11")]
    fn str(#[case] d: Dyadic, #[case] s: &str) {
        assert_eq!(format!("{}", d), s);
    }

    #[rstest]
    #[case(Dyadic::new(3, 12), (3, 12))]
    #[case(Dyadic::new(12, 10), (3, 12))]
    #[case(Dyadic::new(-3, 4), (-3, 4))]
    fn repr(#[case] d: Dyadic, #[case] p: (SignedMantissa, Exponent)) {
        assert_eq!(d.val_and_exp(), p);
    }

    #[rstest]
    #[case(Dyadic::new(32, 0), Dyadic::new(12, 0), Dyadic::new(44, 0))]
    #[case(Dyadic::new(32, 0), Dyadic::new(-12, 0), Dyadic::new(20, 0))]
    #[case(Dyadic::new(12, 10), Dyadic::new(-3, 4), Dyadic::new(12*1024 - 3*16, 0))]
    #[case(Dyadic::new(5, -8), Dyadic::new(3, 2), Dyadic::new(5 + 3*1024, -8))]
    fn add(#[case] d1: Dyadic, #[case] d2: Dyadic, #[case] d3: Dyadic) {
        println!("{:?} + {:?} = {:?}", d1, d2, d3);
        assert_eq!(d1 + d2, d3);
    }

    #[rstest]
    #[case(Dyadic::new(32, 0), Dyadic::new(12, 0), Dyadic::new(32*12, 0))]
    #[case(Dyadic::new(32, 0), Dyadic::new(-12, 0), Dyadic::new(32 * -12, 0))]
    #[case(Dyadic::new(12, 10), Dyadic::new(-3, 4), Dyadic::new(12*1024 * -3*16, 0))]
    #[case(Dyadic::new(5, -8), Dyadic::new(3, 2), Dyadic::new(5 * 3, -6))]
    #[case(Dyadic::new(1, 0), Dyadic::new(0, 0), Dyadic::new(0, 0))]
    fn mul(#[case] d1: Dyadic, #[case] d2: Dyadic, #[case] d3: Dyadic) {
        println!("{:?} * {:?} = {:?}", d1, d2, d3);
        assert_eq!(d1 * d2, d3);
    }

    #[test]
    fn add_overflow() {
        let d1: Dyadic = Dyadic::new(5, 100);
        let d2 = Dyadic::new(5, 0);
        let d3 = d1 + d2;
        assert!(d3.approx());
        assert_ne!(d1, d3);
        assert_eq!(d3.val_and_exp(), (5, 100));
    }

    #[test]
    fn approx() {
        let d1: Dyadic = Dyadic::new(0, 0);
        let d2 = Dyadic::new(1, -200);
        let d3 = Dyadic::new(5, -200);

        assert_abs_diff_eq!(d1, d2);
        assert_abs_diff_eq!(d2, d3);
    }

    #[rstest]
    #[case(Dyadic::new(1, 0), (1,0))]
    #[case(Dyadic::new(0, 0), (0,0))]
    #[case(Dyadic::new(32, 0), (1, 5))]
    #[case(Dyadic::new(-5, 10), (-5, 10))]
    #[case(Dyadic::new(53, 100), (53, 100))]
    fn val_and_exp(#[case] s: Dyadic, #[case] ve: (SignedMantissa, Exponent)) {
        assert_eq!(ve, s.val_and_exp());
        assert_eq!(ve.0, s.val());
        assert_eq!(ve.1, s.exp());
    }

    #[rstest]
    #[case(Dyadic::new(1, 0))]
    #[case(Dyadic::new(0, 0))]
    #[case(Dyadic::new(32, 0))]
    #[case(Dyadic::new(-5, 10))]
    #[case(Dyadic::new(53, 100))]
    fn roundtrip1(#[case] s: Dyadic) {
        let mut s1 = Dyadic::from(f64::try_from(s).unwrap());
        s1.set_approx(false);
        assert_eq!(s, s1);
    }

    #[rstest]
    #[case(0.0)]
    #[case(1.0)]
    #[case(5.0)]
    #[case(0.3)]
    #[case(13e-60)]
    #[case(-55.13)]
    fn roundtrip2(#[case] f: f64) {
        let f1 = f64::try_from(Dyadic::from(f));
        assert_eq!(Ok(f), f1);
    }

    #[rstest]
    #[case(Dyadic::new(1, 1_000_000))]
    fn float_fail(#[case] s: Dyadic) {
        assert!(f64::try_from(s).is_err());
    }
}
