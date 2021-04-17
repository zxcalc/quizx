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

/// Store the (partial) decomposition of a graph into stabilisers
pub struct Decomposer<G: GraphLike> {
    pub stack: Vec<G>,
    pub done: Vec<G>,
    simp: fn (&mut G),
}

impl<'a, G: GraphLike> Decomposer<G> {
    pub fn new(g: &G) -> Decomposer<G> {
        Decomposer {
            stack: vec![g.clone()],
            done: vec![],
            simp: Decomposer::trivial_simp,
        }
    }

    fn trivial_simp(_: &mut G) {}

    /// Decompose up to 6 T gates in the graph on the top of the
    /// stack.
    pub fn decomp_top(&mut self) -> &mut Self {
        let g = self.stack.pop().unwrap();
        let mut t = vec![];

        for v in g.vertices() {
            if *g.phase(v).denom() == 4 { t.push(v); }
            if t.len() == 6 { break; }
        }

        if t.len() == 6 { self.push_bss_decomp(&g, &t) }
        else if t.len() >= 2 { self.push_sym_decomp(&g, &t[0..2]) }
        else if t.len() == 1 { self.push_single_decomp(&g, &t) }
        else {
            println!("done");
            self.done.push(g);
            self
        }
    }

    fn push_decomp(&mut self, fs: &[fn (&G, &[V]) -> G], g: &G, verts: &[V]) -> &mut Self {
        for f in fs {
            let mut g = f(g, verts);
            (self.simp)(&mut g);
            self.stack.push(g);
        }

        self
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
    fn push_bss_decomp(&mut self, g: &G, verts: &[V]) -> &mut Self {
        self.push_decomp(&[
            Decomposer::replace_b60,
            Decomposer::replace_b66,
            Decomposer::replace_e6,
            Decomposer::replace_o6,
            Decomposer::replace_k6,
            Decomposer::replace_phi1,
            Decomposer::replace_phi2,
        ], g, verts)
    }

    /// Perform a decomposition of 2 T gates in the symmetric 2-qubit
    /// space spanned by stabilisers
    fn push_sym_decomp(&mut self, g: &G, verts: &[V]) -> &mut Self {
        self.push_decomp(&[
            Decomposer::replace_bell_s,
            Decomposer::replace_epr,
        ], g, verts)
    }

    /// Replace a single T gate with its decomposition
    fn push_single_decomp(&mut self, g: &G, verts: &[V]) -> &mut Self {
        self.push_decomp(&[
            Decomposer::replace_t0,
            Decomposer::replace_t1,
        ], g, verts)
    }

    fn replace_b60(g: &G, verts: &[V]) -> G {
        // println!("replace_b60");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![-1, 0, 1, 1]);
        for &v in &verts[0..6] { g.add_to_phase(v, Rational::new(-1,4)); }
        g
    }

    fn replace_b66(g: &G, verts: &[V]) -> G {
        // println!("replace_b66");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![-1, 0, 1, -1]);
        for &v in verts { g.add_to_phase(v, Rational::new(3,4)); }
        g
    }

    fn replace_e6(g: &G, verts: &[V]) -> G {
        // println!("replace_e6");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![0, -1, 0, 0]);

        let w = g.add_vertex_with_phase(VType::Z, Rational::one());
        for &v in verts {
            g.add_to_phase(v, Rational::new(1,4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn replace_o6(g: &G, verts: &[V]) -> G {
        // println!("replace_o6");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![-1, 0, -1, 0]);

        let w = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(1,4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn replace_k6(g: &G, verts: &[V]) -> G {
        // println!("replace_k6");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![1, 0, 0, 0]);

        let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1,2));
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1,4));
            g.add_edge_with_type(v, w, EType::N);
        }

        g
    }

    fn replace_phi1(g: &G, verts: &[V]) -> G {
        // println!("replace_phi1");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(3, vec![1, 0, 1, 0]);

        let mut ws = vec![];
        for i in 0..5 {
            let w = g.add_vertex(VType::Z);
            ws.push(w);
            g.add_edge_with_type(verts[i], ws[i], EType::H);
            g.add_edge_with_type(ws[i], verts[5], EType::H);
            g.add_to_phase(verts[i], Rational::new(-1,4));
        }

        g.add_to_phase(verts[5], Rational::new(3,4));

        g.add_edge_with_type(ws[0], ws[2], EType::H);
        g.add_edge_with_type(ws[0], ws[3], EType::H);
        g.add_edge_with_type(ws[1], ws[3], EType::H);
        g.add_edge_with_type(ws[1], ws[4], EType::H);
        g.add_edge_with_type(ws[2], ws[4], EType::H);

        g
    }

    fn replace_phi2(g: &G, verts: &[V]) -> G {
        // print!("replace_phi2 -> ");
        Decomposer::replace_phi1(g, &vec![
                                 verts[0],
                                 verts[1],
                                 verts[3],
                                 verts[4],
                                 verts[5],
                                 verts[2]])
    }

    fn replace_bell_s(g: &G, verts: &[V]) -> G {
        // println!("replace_bell_s");
        let mut g = g.clone();
        g.add_edge(verts[0], verts[1]);
        g.add_to_phase(verts[0], Rational::new(-1,4));
        g.add_to_phase(verts[1], Rational::new(1,4));

        g
    }

    fn replace_epr(g: &G, verts: &[V]) -> G {
        // println!("replace_epr");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::from_phase(Rational::new(1,4));
        let w = g.add_vertex_with_phase(VType::Z, Rational::one());
        for &v in verts {
            g.add_edge_with_type(v, w, EType::H);
            g.add_to_phase(v, Rational::new(-1,4));
        }

        g
    }

    fn replace_t0(g: &G, verts: &[V]) -> G {
        // println!("replace_t0");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![0,1,0,-1]);
        let w = g.add_vertex(VType::Z);
        g.add_edge_with_type(verts[0], w, EType::H);
        g.add_to_phase(verts[0], Rational::new(-1,4));
        g
    }

    fn replace_t1(g: &G, verts: &[V]) -> G {
        // println!("replace_t1");
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![1,0,1,0]);
        let w = g.add_vertex_with_phase(VType::Z, Rational::one());
        g.add_edge_with_type(verts[0], w, EType::H);
        g.add_to_phase(verts[0], Rational::new(-1,4));
        g
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::*;
    use crate::vec_graph::Graph;

    #[test]
    fn bss_scalars() {
        // this test is mainly to record how each of the exact
        // form scalars for the BSS decomposition were computed
        let one = ScalarN::one();
        let om = ScalarN::Exact(0, vec![0, 1, 0, 0]);
        let om2 = &om * &om;
        let om7 = ScalarN::Exact(0, vec![0, 0, 0, -1]);
        assert_eq!(&om * &om7, ScalarN::one());

        let minus = ScalarN::Exact(0, vec![-1, 0, 0, 0]);
        let onefourth = ScalarN::Exact(-2, vec![1, 0, 0, 0]);
        let two = &one + &one;
        let sqrt2 = ScalarN::sqrt2();
        let eight = &two * &two * &two;

        let k6 = &om7 * &two * &om;
        let phi = &om7 * &eight * &sqrt2 * &om2;
        let b60 = &om7 * &minus * &onefourth * (&one + &sqrt2);
        let b66 = &om7 * &onefourth * (&one + (&minus * &sqrt2));
        let o6 = &om7 * &minus * &two * &sqrt2 * &om2;
        let e6 = &om7 * &minus * &two * &om2;

        assert_eq!(b60, ScalarN::Exact(-2, vec![-1, 0, 1, 1]));
        assert_eq!(b66, ScalarN::Exact(-2, vec![-1, 0, 1, -1]));
        assert_eq!(e6, ScalarN::Exact(1, vec![0, -1, 0, 0]));
        assert_eq!(o6, ScalarN::Exact(1, vec![-1, 0, -1, 0]));
        assert_eq!(k6, ScalarN::Exact(1, vec![1, 0, 0, 0]));
        assert_eq!(phi, ScalarN::Exact(3, vec![1, 0, 1, 0]));
    }

    #[test]
    fn single_scalars() {
        let s0 = ScalarN::sqrt2_pow(-1);
        let s1 = ScalarN::from_phase(Rational::new(1,4)) * &s0;
        println!("s0 = {:?}\ns1 = {:?}", s0, s1);
        assert_eq!(s0, ScalarN::Exact(-1, vec![0,1,0,-1]));
        assert_eq!(s1, ScalarN::Exact(-1, vec![1,0,1,0]));
    }

    #[test]
    fn single() {
        let mut g = Graph::new();
        let v = g.add_vertex_with_phase(VType::Z, Rational::new(1,4));
        let w = g.add_vertex(VType::B);
        g.add_edge(v,w);
        g.set_outputs(vec![w]);

        let mut d = Decomposer::new(&g);
        d.decomp_top();
        assert_eq!(d.stack.len(), 2);

        let t = g.to_tensor4();
        let mut tsum = Tensor4::zeros(vec![2]);
        for h in &d.stack { tsum = tsum + h.to_tensor4(); }
        assert_eq!(t, tsum);
    }

    #[test]
    fn sym() {
        let mut g = Graph::new();
        let mut outs = vec![];
        for _ in 0..2 {
            let v = g.add_vertex_with_phase(VType::Z, Rational::new(1,4));
            let w = g.add_vertex(VType::B);
            outs.push(w);
            g.add_edge(v, w);
        }
        g.set_outputs(outs);

        let mut d = Decomposer::new(&g);
        d.decomp_top();
        assert_eq!(d.stack.len(), 2);

        let t = g.to_tensor4();
        let mut tsum = Tensor4::zeros(vec![2; 2]);
        for h in &d.stack { tsum = tsum + h.to_tensor4(); }
        assert_eq!(t, tsum);
    }

    #[test]
    fn bss() {
        let mut g = Graph::new();
        let mut outs = vec![];
        for _ in 0..6 {
            let v = g.add_vertex_with_phase(VType::Z, Rational::new(1,4));
            let w = g.add_vertex(VType::B);
            outs.push(w);
            g.add_edge(v, w);
        }
        g.set_outputs(outs);

        let mut d = Decomposer::new(&g);
        d.decomp_top();
        assert_eq!(d.stack.len(), 7);

        let t = g.to_tensor4();
        let mut tsum = Tensor4::zeros(vec![2; 6]);
        for h in &d.stack { tsum = tsum + h.to_tensor4(); }
        assert_eq!(t, tsum);
    }

    #[test]
    fn mixed() {
        let mut g = Graph::new();
        let mut outs = vec![];
        for _ in 0..9 {
            let v = g.add_vertex_with_phase(VType::Z, Rational::new(1,4));
            let w = g.add_vertex(VType::B);
            outs.push(w);
            g.add_edge(v, w);
        }
        g.set_outputs(outs);

        let mut d = Decomposer::new(&g);
        while d.stack.len() > 0 { d.decomp_top(); }

        assert_eq!(d.done.len(), 7*2*2);

        // thorough but SLOW
        // let t = g.to_tensor4();
        // let mut tsum = Tensor4::zeros(vec![2; 9]);
        // for h in &d.done { tsum = tsum + h.to_tensor4(); }
        // assert_eq!(t, tsum);
    }
}
