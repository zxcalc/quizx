use approx::AbsDiffEq;
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

pub type Mantissa = u32;
pub type DoubleMantissa = u64;
pub type Exponent = i32;
pub type SignedVal = i32;

const SIGN: u8 = 0x01;
const SIGN_OFF: u8 = 0xff ^ SIGN;
const APPROX: u8 = 0x02;
// const APPROX_OFF: u8 = 0xff ^ APPROX;

// A dyadic is essentially a floating point number, but we have more control over precision
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Dyadic {
    flags: u8,
    exp: Exponent,
    val: Mantissa,
}

impl Dyadic {
    pub fn new(val: SignedVal, exp: Exponent) -> Self {
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
    pub fn val_and_exp(&self) -> (SignedVal, Exponent) {
        let shift = self.val.trailing_zeros();
        let v = self.val.wrapping_shr(shift) as SignedVal;
        if self.sign() {
            (-v, self.exp + (shift as Exponent))
        } else {
            (v, self.exp + (shift as Exponent))
        }
    }

    #[inline]
    fn normalize(&mut self) {
        if self.val == 0 {
            self.exp = 0;
            self.flags &= SIGN_OFF;
        } else {
            let head = self.val.leading_zeros();
            self.exp = self.exp - (head as Exponent);
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

impl fmt::Display for Dyadic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (mut v, mut e) = self.val_and_exp();

        let bnd: SignedVal = 1024;
        if v > -bnd && v < bnd && e.is_positive() && e < 10 {
            v = v * (2 as SignedVal).pow(e as u32);
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
        // first compute the order as if "self" was positive
        let ord = if self.sign() != other.sign() {
            Ordering::Greater
        } else if self.exp == other.exp {
            self.val.cmp(&other.val)
        } else {
            self.exp.cmp(&other.exp)
        };

        // if "self" is actually negative, flip the order
        Some(if self.sign() { ord.reverse() } else { ord })
    }
}

impl Ord for Dyadic {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
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

impl Add<Dyadic> for Dyadic {
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
                self.val = self.val - rhs.val;
                self.flags = self.flags | (rhs.flags & APPROX);
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
                self.val = self.val + rhs.val;
                self.exp = self.exp + 1;
            }

            self.flags = self.flags | (rhs.flags & APPROX);
        };

        self
    }
}

impl Sub for Dyadic {
    type Output = Dyadic;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self + rhs.neg()
    }
}

impl Mul for Dyadic {
    type Output = Dyadic;

    fn mul(mut self, rhs: Self) -> Self::Output {
        if self.val == 0 {
            return rhs;
        } else if rhs.val == 0 {
            return self;
        }

        self.exp = self.exp + rhs.exp;
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
            self.exp = self.exp + (shift as Exponent);
        }

        self.val = v as Mantissa;
        self
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
    fn repr(#[case] d: Dyadic, #[case] p: (i32, i32)) {
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
}
