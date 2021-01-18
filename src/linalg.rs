// QuiZX - Rust library for quantum circuit rewriting and optimisation
//         using the ZX-calculus
// Copyright (C) 2021 - Aleks Kissinger
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

//! Matrices and linear algebra over F2

use std::fmt;

/// A type for matrices over F2
pub struct MatF2 {
    d: Vec<Vec<i32>>
}

pub trait RowColOps {
    fn row_add(&mut self, r0: usize, r1: usize);
    fn col_add(&mut self, c0: usize, c1: usize);
    fn row_swap(&mut self, r0: usize, r1: usize);
    fn col_swap(&mut self, c0: usize, c1: usize);
}


impl MatF2 {
    pub fn zero(rows: usize, cols: usize) -> MatF2 {
        MatF2 { d: vec![vec![0; cols]; rows] }
    }

    pub fn id(dim: usize) -> MatF2 {
        MatF2 {
            d: (0..dim).map(|x| (0..dim).map(|y|
                   if x == y { 1 } else { 0 })
               .collect()).collect()
        }
    }

    pub fn num_rows(&self) -> usize {
        self.d.len()
    }

    pub fn num_cols(&self) -> usize {
        if self.d.len() > 0 { self.d[0].len() }
        else { 0 }
    }
}

impl RowColOps for MatF2 {
    fn row_add(&mut self, r0: usize, r1: usize) {
        for i in 0..self.num_cols() {
            self.d[r0][i] += self.d[r1][i];
        }
    }

    fn col_add(&mut self, c0: usize, c1: usize) {
        for i in 0..self.num_rows() {
            self.d[i][c0] += self.d[i][c1];
        }
    }

    fn row_swap(&mut self, r0: usize, r1: usize) {
        self.d.swap(r0, r1);
    }

    fn col_swap(&mut self, c0: usize, c1: usize) {
        for i in 0..self.num_rows() {
            self.d[i].swap(c0, c1);
        }
    }
}

impl fmt::Display for MatF2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in &self.d {
            write!(f, "[ ")?;
            for x in row { write!(f, "{} ", x)?; }
            write!(f, "]\n")?;
        }
        Ok(())
    }
}

impl std::ops::Index<(usize,usize)> for MatF2 {
    type Output = i32;
    fn index(&self, idx: (usize,usize)) -> &Self::Output { &self.d[idx.0][idx.1] }
}

impl std::ops::IndexMut<(usize,usize)> for MatF2 {
    fn index_mut(&mut self, idx: (usize,usize)) -> &mut Self::Output { &mut self.d[idx.0][idx.1] }
}

impl std::ops::Index<usize> for MatF2 {
    type Output = Vec<i32>;
    fn index(&self, idx: usize) -> &Self::Output { &self.d[idx] }
}

impl std::ops::IndexMut<usize> for MatF2 {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output { &mut self.d[idx] }
}
