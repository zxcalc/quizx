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

use num::Rational;
use crate::graph::*;
use crate::scalar::*;

pub struct Decomposer<G: GraphLike> {
    pub stack: Vec<G>,
}

impl<'a, G: GraphLike> Decomposer<G> {
    pub fn new(g: &G) -> Decomposer<G> {
        Decomposer { stack: vec![g.clone()] }
    }

    pub fn decomp_top(&mut self) -> &mut Self {
        let g = self.stack.pop().unwrap();
        let mut t = vec![];

        for v in g.vertices() {
            if *g.phase(v).denom() == 4 { t.push(v); }
            if t.len() == 6 { break; }
        }

        if t.len() == 6 { self.push_bss_decomp(&g, t) }
        else { panic!("decompose fewer than 6 T gates not implemented") }
    }

    /// Perform the Bravyi-Smith-Smolin decomposition of 6 T gates
    /// into a sum of 7 terms
    ///
    /// See Section IV of:
    /// https://journals.aps.org/prx/pdf/10.1103/PhysRevX.6.021043
    ///
    /// In particular, see the text below equation (10) and
    /// equation (11) itself.
    ///
    fn push_bss_decomp(&mut self, g: &G, verts: Vec<V>) -> &mut Self {
        self.stack.push(Decomposer::replace_b60(g, &verts));
        self.stack.push(Decomposer::replace_b66(g, &verts));
        self.stack.push(Decomposer::replace_e6(g, &verts));
        self.stack.push(Decomposer::replace_o6(g, &verts));
        self.stack.push(Decomposer::replace_k6(g, &verts));
        self.stack.push(Decomposer::replace_phi1(g, &verts));
        self.stack.push(Decomposer::replace_phi2(g, &verts));
        self
    }

    fn replace_b60(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![29, -1, -29, 0]);
        for &v in verts { g.add_to_phase(v, Rational::new(-1,4)); }
        g
    }

    fn replace_b66(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![169, 1, -169, 0]);
        for &v in verts { g.add_to_phase(v, Rational::new(3,4)); }
        g
    }

    fn replace_e6(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(2, vec![0, -7, 10, -7]);

        let w = g.add_vertex_with_phase(VType::Z, Rational::one());
        for &v in verts {
            g.add_to_phase(v, Rational::new(1,4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn replace_o6(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(0, vec![99, -70, -1, 70]);

        let w = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(1,4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn replace_k6(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![0, 0, 0, 99]);

        let w = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1,4));
            g.add_edge_with_type(v, w, EType::N);
        }

        g
    }

    fn replace_phi1(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(3, vec![-1, 140, 1, 0]);

        let ws: Vec<_> = verts.iter().map(|&v| {
            g.add_to_phase(v, Rational::new(-1,4));
            g.add_vertex_with_phase(VType::Z, Rational::one())
        }).collect();

        for i in 0..5 {
            g.add_edge_with_type(verts[i], ws[i], EType::H);
            g.add_edge_with_type(ws[i], ws[5], EType::H);
        }

        g.add_edge_with_type(ws[0], ws[2], EType::H);
        g.add_edge_with_type(ws[0], ws[3], EType::H);
        g.add_edge_with_type(ws[1], ws[3], EType::H);
        g.add_edge_with_type(ws[1], ws[4], EType::H);
        g.add_edge_with_type(ws[2], ws[4], EType::H);

        g
    }

    fn replace_phi2(g: &G, verts: &Vec<V>) -> G {
        Decomposer::replace_phi1(g, &vec![
                                 verts[0],
                                 verts[1],
                                 verts[3],
                                 verts[4],
                                 verts[5],
                                 verts[2]])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalars() {
        // this test is just to record where the BSS scalars
        // came from. These were all globals at the top of
        // simulate.py in PyZX, some of which are sometimes 
        // combined with a rt2 power and a phase in the
        // replace_XXXX function. We do this all at once.

        // 1/(2+2j) = (2-2j)/8
        let mut global = Scalar4::Exact(-3, [2, 0, -2, 0]);
        // times -(7 + 5*rt2)
        global *= Scalar4::Exact(0, [-7, -5, 0, -5]);

        // -1/8 * (-16 + 12*rt2)
        let b60 = Scalar4::Exact(-3, [-16, 12, 0, 12]) * &global;

        // -1/8 * (96 - 68*rt2)
        let b66 = Scalar4::Exact(-3, [-96, 68, 0, 68]) * &global;

        // 4i * (10 - 7 rt(2))
        let e6 = ScalarN::Exact(2, vec![0, -7, 10, -7]); // e6

        // 2i * (-14 + 10*rt2)
        let o6 = Scalar4::Exact(1, [0, 10, -14, 10]) * &global;

        // rt2^5 * (7 - 5*rt2)
        let mut k6 = Scalar4::Exact(0, [7, -5, 0, -5]) * &global;
        k6.mul_sqrt2_pow(5);

        // -i rt2^9 * (10 - 7*rt2)
        let mut phi = Scalar4::Exact(0, [10, -7, 0, -7]) * &global;
        phi.mul_sqrt2_pow(9);
        phi.mul_phase(Rational::new(3,2));

        println!("ScalarN::{:?}; // b60", b60);
        println!("ScalarN::{:?}; // b66", b66);
        println!("ScalarN::{:?}; // e6", e6);
        println!("ScalarN::{:?}; // o6", o6);
        println!("ScalarN::{:?}; // k6", k6);
        println!("ScalarN::{:?}; // phi", phi);

        // TODO: test these are indeed the scalars each term gets
    }
}
