use num::Zero;
use std::ops::{Add, Index};

pub type Var = u16;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Parity(Vec<Var>);

impl Parity {
    pub fn single(v: Var) -> Self {
        Parity(vec![v])
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_one(&self) -> bool {
        self.len() == 1 && self[0] == 0
    }

    pub fn one() -> Self {
        Parity(vec![0])
    }
}

impl Index<usize> for Parity {
    type Output = Var;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl From<Vec<Var>> for Parity {
    fn from(value: Vec<Var>) -> Self {
        Parity(value)
    }
}

impl From<Parity> for Vec<Var> {
    fn from(value: Parity) -> Self {
        value.0
    }
}

impl Zero for Parity {
    fn is_zero(&self) -> bool {
        self.len() == 0
    }

    fn zero() -> Self {
        Parity(vec![])
    }
}

impl Add<&Parity> for &Parity {
    type Output = Parity;

    fn add(self, rhs: &Parity) -> Self::Output {
        let mut vars = Vec::with_capacity(self.len() + rhs.len());

        // merge the two vecs, keeping sorted and dropping duplicates
        let mut i = 0;
        let mut j = 0;
        loop {
            if i < self.len() {
                if j < rhs.len() {
                    if self[i] < rhs[j] {
                        vars.push(self[i]);
                        i += 1;
                    } else if self[i] > rhs[j] {
                        vars.push(rhs[j]);
                        j += 1;
                    } else {
                        i += 1;
                        j += 1;
                    }
                } else {
                    vars.push(self[i]);
                    i += 1;
                }
            } else {
                if j < rhs.len() {
                    vars.push(rhs[j]);
                    j += 1;
                } else {
                    break;
                }
            }
        }

        Parity(vars)
    }
}

impl Add<Parity> for Parity {
    type Output = Parity;
    fn add(self, rhs: Parity) -> Self::Output {
        &self + &rhs
    }
}

pub struct Factor {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parities() {
        let p1: Parity = vec![1, 3, 6, 8].into();
        let p2: Parity = vec![2, 4, 9].into();
        let p3: Parity = vec![1, 2, 3, 4, 6, 8, 9].into();
        let p2a: Parity = vec![2, 3, 4, 9].into();
        let p3a: Parity = vec![1, 2, 4, 6, 8, 9].into();
        assert_eq!(&p1 + &p2, p3);
        assert_eq!(&p1 + &p2a, p3a);
    }
}
