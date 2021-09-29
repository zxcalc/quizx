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
use std::collections::VecDeque;
use rand::{thread_rng, Rng};
use rayon::prelude::*;
use crate::graph::*;
use crate::scalar::*;

#[derive(Copy,Clone,PartialEq,Eq,Debug)]
pub enum SimpFunc {
    FullSimp,
    CliffordSimp,
    NoSimp,
}
use SimpFunc::*;

/// Store the (partial) decomposition of a graph into stabilisers
#[derive(Clone)]
pub struct Decomposer<G: GraphLike> {
    pub stack: VecDeque<(usize,G)>,
    pub done: Vec<G>,
    pub scalar: ScalarN,
    pub nterms: usize,
    simp_func: SimpFunc,
    random_t: bool,
    use_cats: bool,
    save: bool, // save graphs on 'done' stack
}

// impl<G: GraphLike> Send for Decomposer<G> {}

/// Gives upper bound for number of terms needed for BSS decomposition
///
/// Note this number can be very large. We use a float here to avoid overflows.
pub fn terms_for_tcount(tcount: usize) -> f64 {
    let mut t = tcount as i32;
    let mut count = 7f64.powi(t / 6i32);
    t = t % 6;
    count *= 2f64.powi(t / 2i32);
    if t % 2 == 1 { count *= 2.0; }
    return count;
}

impl<'a, G: GraphLike> Decomposer<G> {
    pub fn empty() -> Decomposer<G> {
        Decomposer {
            stack: VecDeque::new(),
            done: vec![],
            scalar: ScalarN::zero(),
            nterms: 0,
            simp_func: NoSimp,
            random_t: false,
            use_cats: false,
            save: false,
        }
    }

    pub fn new(g: &G) -> Decomposer<G> {
        let mut d = Decomposer::empty();
        d.stack.push_back((0, g.clone()));
        d
    }

    /// Split a Decomposer with N graphs on the stack into N Decomposers
    /// with 1 graph each.
    ///
    /// Used for parallelising. The last decomposer in the list keeps the
    /// current state (e.g. `nterms` and `scalar`).
    pub fn split(mut self) -> Vec<Decomposer<G>> {
        let mut ds = vec![];
        while self.stack.len() > 1 {
            let (_,g) = self.stack.pop_front().unwrap();
            let mut d1 = Decomposer::new(&g);
            d1.save(self.save)
              .random_t(self.random_t)
              .with_simp(self.simp_func);
            ds.push(d1);
        }
        ds.push(self);
        ds
    }

    /// Merge N decomposers into 1, adding scalars together
    pub fn merge(mut ds: Vec<Decomposer<G>>) -> Decomposer<G> {
        if let Some(mut d) = ds.pop() {
            while !ds.is_empty() {
                let d1 = ds.pop().unwrap();
                d.scalar = d.scalar + d1.scalar;
                d.nterms += d1.nterms;
                d.stack.extend(d1.stack);
                d.done.extend(d1.done);
            }
            d
        } else {
            Decomposer::empty()
        }
    }

    // pub fn seed(&mut self, seed: u64) -> &mut Self { self.rng = StdRng::seed_from_u64(seed); self }

    pub fn with_simp(&mut self, f: SimpFunc) -> &mut Self {
        self.simp_func = f;
        self
    }

    pub fn with_full_simp(&mut self) -> &mut Self { self.with_simp(FullSimp) }
    pub fn with_clifford_simp(&mut self) -> &mut Self { self.with_simp(CliffordSimp) }

    pub fn random_t(&mut self, b: bool) -> &mut Self {
        self.random_t = b;
        self
    }

    pub fn use_cats(&mut self, b: bool) -> &mut Self {
        self.use_cats = b;
        self
    }

    pub fn save(&mut self, b: bool) -> &mut Self {
        self.save = b;
        self
    }

    /// Computes `terms_for_tcount` for every graph on the stack
    pub fn max_terms(&self) -> f64 {
        let mut n = 0.0;
        for (_,g) in &self.stack {
            n += terms_for_tcount(g.tcount());
        }

        n
    }


    pub fn pop_graph(&mut self) -> G {
        let (_, g) = self.stack.pop_back().unwrap();
        g
    }

    /// Decompose the first <= 6 T gates in the graph on the top of the
    /// stack.
    pub fn decomp_top(&mut self) -> &mut Self {
        let (depth, g) = self.stack.pop_back().unwrap();
        if self.use_cats {
            let cat_nodes = Decomposer::cat_ts(&g);//gadget_ts(&g);
            //println!("{:?}", gadget_nodes);
            //let nts = cat_nodes.iter().fold(0, |acc, &x| if g.phase(x).denom() == &4 { acc + 1 } else { acc });
            if cat_nodes.len() > 0 {
                // println!("using cat!");
                return self.push_cat_decomp(depth+1, &g, &cat_nodes)                
            }
            let ts = Decomposer::first_ts(&g);
            if ts.len()>=5 {
                return self.push_magic5_from_cat_decomp(depth+1, &g, &ts[..5])
            }
        }
        let ts = if self.random_t { Decomposer::random_ts(&g, &mut thread_rng()) }
                 else { Decomposer::first_ts(&g) };
        self.decomp_ts(depth, g, &ts);
        self
    }

    /// Decompose until there are no T gates left
    pub fn decomp_all(&mut self) -> &mut Self {
        while self.stack.len() > 0 { self.decomp_top(); }
        self
    }

    /// Decompose breadth-first until the given depth
    pub fn decomp_until_depth(&mut self, depth: usize) -> &mut Self {
        while self.stack.len() > 0 {
            // pop from the bottom of the stack to work breadth-first
            let (d, g) = self.stack.pop_front().unwrap();
            if d >= depth {
                self.stack.push_front((d,g));
                break;
            } else {
                if self.use_cats {
                    let cat_nodes = Decomposer::cat_ts(&g);//gadget_ts(&g);
                    //println!("{:?}", gadget_nodes);
                    let nts = cat_nodes.iter().fold(0, |acc, &x| if g.phase(x).denom() == &4 { acc + 1 } else { acc });
                    if nts > 2 {
                        // println!("using cat!");
                        return self.push_cat_decomp(depth+1, &g, &cat_nodes)                
                    }
                }
                let ts = if self.random_t { Decomposer::random_ts(&g, &mut thread_rng()) }
                         else { Decomposer::first_ts(&g) };
                self.decomp_ts(d, g, &ts);
            }
        }
        self
    }

    /// Decompose in parallel, starting at the given depth
    pub fn decomp_parallel(mut self, depth: usize) -> Self {
        self.decomp_until_depth(depth);
        let ds = self.split();
        Decomposer::merge(ds.into_par_iter().map(|mut d| {
            d.decomp_all(); d
        }).collect())
    }

    pub fn decomp_ts(&mut self, depth: usize, g: G, ts: &[usize]) {
        if ts.len() == 6 { self.push_bss_decomp(depth+1, &g, ts); }
        else if ts.len() >= 2 { self.push_sym_decomp(depth+1, &g, &ts[0..2]); }
        else if ts.len() > 0 { self.push_single_decomp(depth+1, &g, ts); }
        else {
            // crate::simplify::full_simp(&mut g);
            self.scalar = &self.scalar + g.scalar();
            self.nterms += 1;
            if g.num_vertices() != 0 {
                println!("{}", g.to_dot());
                println!("WARNING: graph was not fully reduced");
                // println!("{}", g.to_dot());
            }
            if self.save { self.done.push(g); }
        }
    }

    /// Pick the first <= 6 T gates from the given graph
    pub fn first_ts(g: &G) -> Vec<V> {
        let mut t = vec![];

        for v in g.vertices() {
            if *g.phase(v).denom() == 4 { t.push(v); }
            if t.len() == 6 { break; }
        }

        t
    }

    /// Pick <= 6 T gates from the given graph, chosen at random
    pub fn random_ts(g: &G, rng: &mut impl Rng) -> Vec<V> {
        let mut all_t: Vec<_> = g.vertices().filter(|&v| *g.phase(v).denom() == 4).collect();
        let mut t = vec![];

        while t.len() < 6 && all_t.len() > 0 {
            let i = rng.gen_range(0..all_t.len());
            t.push(all_t.swap_remove(i));
        }

        t
    }

    /// Returns a best occurrence of a cat state
    /// The fist vertex in the result is the Clifford spider
    pub fn cat_ts(g: &G) -> Vec<V> {
        // the graph g is supposed to be completely simplified
        let prefered_order = [4,6,5,3];
        let mut res = vec![];
        let mut index = None;
        for v in g.vertices() {
            if g.phase(v).denom() == &1{
                let mut neigh = g.neighbor_vec(v);
                if neigh.len() <= 6 {
                    match prefered_order.iter().position(|&r| r == neigh.len()) {
                        Some(this_ind) => match index {
                            Some(ind) if this_ind < ind => {res = vec![v]; res.append(&mut neigh); index = Some(this_ind);},
                            None => {res = vec![v]; res.append(&mut neigh); index = Some(this_ind);},
                            _ => (),
                        },
                        _ => (),
                    }
                    if index == Some(0){break;}
                }
            }
        }
        res
    }

    fn push_decomp(&mut self, fs: &[fn (&G, &[V]) -> G], depth: usize, g: &G, verts: &[V]) -> &mut Self {
        for f in fs {
            let mut g = f(g, verts);
            match self.simp_func {
                FullSimp => { crate::simplify::full_simp(&mut g); },
                CliffordSimp => { crate::simplify::clifford_simp(&mut g); },
                _ => {}
            }

            // let comps = g.component_vertices();
            // if comps.len() > 1 {
            //     println!("GOT {} COMPONENTS ({})", comps.len(), comps.iter().map(|c| c.len()).format(","));
            // }
            self.stack.push_back((depth, g));
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
    fn push_bss_decomp(&mut self, depth: usize, g: &G, verts: &[V]) -> &mut Self {
        self.push_decomp(&[
            Decomposer::replace_b60,
            Decomposer::replace_b66,
            Decomposer::replace_e6,
            Decomposer::replace_o6,
            Decomposer::replace_k6,
            Decomposer::replace_phi1,
            Decomposer::replace_phi2,
        ], depth, g, verts)
    }

    /// Perform a decomposition of 2 T gates in the symmetric 2-qubit
    /// space spanned by stabilisers
    fn push_sym_decomp(&mut self, depth: usize, g: &G, verts: &[V]) -> &mut Self {
        self.push_decomp(&[
            Decomposer::replace_bell_s,
            Decomposer::replace_epr,
        ], depth, g, verts)
    }

    /// Replace a single T gate with its decomposition
    fn push_single_decomp(&mut self, depth: usize, g: &G, verts: &[V]) -> &mut Self {
        self.push_decomp(&[
            Decomposer::replace_t0,
            Decomposer::replace_t1,
        ], depth, g, verts)
    }

    /// Perform a decomposition of 5 T-spiders, with one remaining
    fn push_magic5_from_cat_decomp(&mut self, depth: usize, g: &G, verts: &[V]) -> &mut Self {
        //println!("magic5");
        self.push_decomp(&[
            Decomposer::replace_magic5_0,
            Decomposer::replace_magic5_1,
            Decomposer::replace_magic5_2,
        ], depth, &g, &verts)
    }

    /// Perform a decomposition of cat states
    fn push_cat_decomp(&mut self, depth: usize, g: &G, verts: &[V]) -> &mut Self {
        // verts[0] is a 0- or pi-spider, linked to all and only to vs in verts[1..] which are T-spiders
        let mut g = g.clone(); // that is annoying ... 
        let mut verts = Vec::from(verts);
        if g.phase(verts[0]).numer() == &1 {
            g.set_phase(verts[0], Rational::new(0,1));
            let mut neigh = g.neighbor_vec(verts[1]);
            neigh.retain(|&x| x != verts[0]);
            for &v in &neigh{
                g.add_to_phase(v, Rational::new(1,1));
            }
            let tmp = g.phase(verts[1]);
            *g.scalar_mut() *= ScalarN::from_phase(tmp);
            g.set_phase(verts[1], g.phase(verts[1])*Rational::new(-1,1));
        }
        if [3,5].contains(&verts[1..].len()) {
            let w = g.add_vertex(VType::Z);
            let v = g.add_vertex(VType::Z);
            g.add_edge_with_type(v, w, EType::H);
            g.add_edge_with_type(v, verts[0], EType::H);
            verts.push(v);     
        }
        if verts[1..].len() == 6 {
            self.push_decomp(&[
                Decomposer::replace_cat6_0,
                Decomposer::replace_cat6_1,
                Decomposer::replace_cat6_2,
            ], depth, &g, &verts)
        }else if verts[1..].len() == 4 {
            self.push_decomp(&[
                Decomposer::replace_cat4_0,
                Decomposer::replace_cat4_1,
            ], depth, &g, &verts)
        }else { println!("this shouldn't be printed"); self }
    }

    fn replace_cat6_0(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, 0, 0]);    
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1,4));
            g.set_edge_type(v, verts[0], EType::N);
        }
        g.set_phase(verts[0], Rational::new(-1,2));
        g
    }

    fn replace_cat6_1(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();  
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![-1, 0, 1, 0]);    
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1,4));
        }
        g
    }

    fn replace_cat6_2(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(7, vec![0, -1, 0, 0]);
        for i in 1..verts.len() {
            g.add_to_phase(verts[i], Rational::new(-1,4));
            for j in i+1..verts.len() {
                g.add_edge_smart(verts[i], verts[j], EType::H);
            }
        }
        g
    }

    fn replace_magic5_0(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![1, 0, 0, 0]);    
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1,4));
            g.add_edge_smart(v, verts[0], EType::N);
        }
        g.add_to_phase(verts[0], Rational::new(-3,4));
        g
    }

    fn replace_magic5_1(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();  
        *g.scalar_mut() *= ScalarN::Exact(1, vec![-1, 0, 1, 0]);
        let p = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1,4));
            g.add_edge_with_type(v, p, EType::H);
        }
        let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1,4));
        g.add_edge_with_type(w, p, EType::H);
        g
    }

    fn replace_magic5_2(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(9, vec![0, -1, 0, 0]);
        let p = g.add_vertex(VType::Z);
        let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1,4));
        g.add_edge_with_type(p, w, EType::H);
        for i in 0..verts.len() {
            g.add_to_phase(verts[i], Rational::new(-1,4));
            g.add_edge_with_type(verts[i], p, EType::H);
            g.add_edge_with_type(verts[i], w, EType::H);
            for j in i+1..verts.len() {
                g.add_edge_smart(verts[i], verts[j], EType::H);
            }
        }
        g
    }

    fn replace_cat4_0(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(0, vec![0, 0, 1, 0]);    
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1,4));
        }
        g
    }

    fn replace_cat4_1(g: &G, verts: &[V]) -> G {
        // same as replace_cat6_0, only with a different scalar
        let mut g = g.clone();  
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, -1, 0]);    
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1,4));
            g.set_edge_type(v, verts[0], EType::N);
        }
        g.set_phase(verts[0], Rational::new(-1,2));
        g
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
        g.add_edge_smart(verts[0], verts[1], EType::N);
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
        for (_,h) in &d.stack { tsum = tsum + h.to_tensor4(); }
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
        for (_,h) in &d.stack { tsum = tsum + h.to_tensor4(); }
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
        for (_,h) in &d.stack { tsum = tsum + h.to_tensor4(); }
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
        d.save(true);
        assert_eq!(d.max_terms(), 7.0*2.0*2.0);
        while d.stack.len() > 0 { d.decomp_top(); }

        assert_eq!(d.done.len(), 7*2*2);

        // thorough but SLOW
        // let t = g.to_tensor4();
        // let mut tsum = Tensor4::zeros(vec![2; 9]);
        // for h in &d.done { tsum = tsum + h.to_tensor4(); }
        // assert_eq!(t, tsum);
    }

    #[test]
    fn mixed_sc() {
        let mut g = Graph::new();
        for i in 0..11 {
            g.add_vertex_with_phase(VType::Z, Rational::new(1,4));

            for j in 0..i {
                g.add_edge_with_type(i, j, EType::H);
            }
            // let w = g.add_vertex(VType::Z);
            // g.add_edge(v, w);
        }


        let mut d = Decomposer::new(&g);
        d.with_full_simp();
        // assert_eq!(d.max_terms(), 7.0*2.0*2.0);
        d.decomp_all();
        // assert_eq!(d.nterms, 7*2*2);

        let sc = g.to_tensor4()[[]];
        assert_eq!(Scalar::from_scalar(&sc), d.scalar);
    }


    #[test]
    fn all_and_depth() {
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
        d.with_full_simp();
        d.save(true).decomp_all();
        assert_eq!(d.done.len(), 7*2*2);
        let mut d = Decomposer::new(&g);
        d.with_full_simp();
        d.decomp_until_depth(2);
        assert_eq!(d.stack.len(), 7*2);
    }

    #[test]
    fn full_simp() {
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
        d.with_full_simp()
         .save(true)
         .decomp_all();
        assert_eq!(d.done.len(), 7*2*2);
    }
}
