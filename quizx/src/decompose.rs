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
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![29, -1, -29, 0]); // b60

        for &v in verts {
            g.add_to_phase(v, Rational::new(-1,4));
        }

        g
    }

    fn replace_b66(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![169, 1, -169, 0]); // b66

        for &v in verts {
            g.add_to_phase(v, Rational::new(3,4));
        }

        g
    }

    fn replace_e6(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        // 4i * (10 - 7 rt(2))
        *g.scalar_mut() *= ScalarN::Exact(2, vec![0, -7, 10, -7]);

        let w = g.add_vertex(VType::Z);
        g.set_phase(w, Rational::one());

        for &v in verts {
            g.add_to_phase(v, Rational::new(1,4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn replace_o6(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(0, vec![99, -70, -1, 70]); // o6

        let w = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(1,4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn replace_k6(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![0, 0, 0, 99]); // k6

        let w = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1,4));
            g.add_edge_with_type(v, w, EType::N);
        }

        g
    }

    fn replace_phi1(g: &G, verts: &Vec<V>) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(3, vec![-1, 140, 1, 0]); // phi

        let ws: Vec<_> = (0..6).map(|_| {
            let w = g.add_vertex(VType::Z);
            g.set_phase(w, Rational::one());
            w
        }).collect();

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
