// QuiZX - Rust library for quantum circuit rewriting and optimisation
//         using the ZX-calculus
// Copyright (C) 2025 - Aleks Kissinger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use num::Zero;
use std::{
    cmp::Ordering,
    ops::{Add, Index},
};

pub type Var = u16;

/// A representation for an XOR of variables, represented as unsigned
/// integers. The variable "0" is reserved to mean the constant "1".
///
/// For example [0, 3, 4] means 1 ⊕ b3 ⊕ b4.
///
/// Variables are kept sorted to ensure uniqueness.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Parity(Vec<Var>);

/// A boolean expression, represented as a conjunction of XORs
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Expr(Vec<Parity>);

impl Parity {
    pub fn single(v: Var) -> Self {
        Parity(vec![v])
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn is_one(&self) -> bool {
        self.len() == 1 && self[0] == 0
    }

    pub fn one() -> Self {
        Parity(vec![0])
    }

    /// Returns of a copy of the parity negated
    pub fn negated(&self) -> Self {
        &Parity::one() + self
    }
}

impl Index<usize> for Parity {
    type Output = Var;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

/// Converts a vec into a Parity. This will sort the variables, but
/// it does not automatically cancel out duplicates.
impl From<Vec<Var>> for Parity {
    fn from(mut value: Vec<Var>) -> Self {
        value.sort();
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
        self.is_empty()
    }

    fn zero() -> Self {
        Parity(vec![])
    }
}
impl Add<&Parity> for &Parity {
    type Output = Parity;

    #[allow(clippy::collapsible_else_if)]
    fn add(self, rhs: &Parity) -> Self::Output {
        let mut vars = Vec::with_capacity(self.len() + rhs.len());

        // merge the two vecs, keeping sorted and dropping duplicates
        let mut i = 0;
        let mut j = 0;
        loop {
            if i < self.len() {
                if j < rhs.len() {
                    match self[i].cmp(&rhs[j]) {
                        Ordering::Less => {
                            vars.push(self[i]);
                            i += 1;
                        }
                        Ordering::Greater => {
                            vars.push(rhs[j]);
                            j += 1;
                        }
                        Ordering::Equal => {
                            i += 1;
                            j += 1;
                        }
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

impl Expr {
    pub fn linear(p: impl Into<Parity>) -> Self {
        Expr(vec![p.into()])
    }

    pub fn quadratic(p1: impl Into<Parity>, p2: impl Into<Parity>) -> Self {
        let mut p1 = p1.into();
        let mut p2 = p2.into();
        if p1 > p2 {
            (p1, p2) = (p2, p1);
        }

        if p1.is_one() || p1 == p2 {
            Expr(vec![p2])
        } else {
            Expr(vec![p1, p2])
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn is_linear(&self) -> bool {
        self.len() == 1
    }
}

impl Index<usize> for Expr {
    type Output = Parity;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
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
