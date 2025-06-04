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

// use crate::decompose;
use crate::fscalar::*;
use crate::graph::*;
// use crate::hash_graph::Graph;
// use crate::tensor::Tensor;
// use itertools::Itertools;
// use itertools::Itertools;
use num::Rational64;
use rand::{thread_rng, Rng};
// use rand::rngs::StdRng;
use rayon::prelude::*;

/// Gives upper bound for number of terms needed for BSS decomposition
///
/// Note this number can be very large. We use a float here to avoid overflows.
pub fn terms_for_tcount(tcount: usize) -> f64 {
    let mut t = tcount as i32;
    let mut count = 7f64.powi(t / 6i32);
    t %= 6;
    count *= 2f64.powi(t / 2i32);
    if t % 2 == 1 {
        count *= 2.0;
    }
    count
}

/// Pick the first <= 6 T gates from the given graph
pub fn first_ts<G: GraphLike>(g: &G) -> Vec<V> {
    let mut t = vec![];

    for v in g.vertices() {
        if g.phase(v).is_t() {
            t.push(v);
        }
        if t.len() == 6 {
            break;
        }
    }
    t
}

/// Pick <= 6 T gates from the given graph, chosen at random
pub fn random_ts<G: GraphLike>(g: &G, rng: &mut impl Rng) -> Vec<V> {
    // the graph g is assumed to contain no X spiders
    let mut all_t: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();
    let mut t = vec![];

    while t.len() < 6 && !all_t.is_empty() {
        let i = rng.gen_range(0..all_t.len());
        t.push(all_t.swap_remove(i));
    }
    t
}

/// Returns a best occurrence of a cat state
/// The fist vertex in the result is the Pauli spider
pub fn cat_ts<G: GraphLike>(g: &G) -> Vec<V> {
    // the graph g is assumed to be graph-like
    let preferred_order = [4, 6, 5, 3];
    let mut res = vec![];
    let mut index = None;
    for v in g.vertices() {
        if g.vertex_type(v) == VType::Z && g.phase(v).is_pauli() {
            let mut neigh = g.neighbor_vec(v);
            if neigh.len() <= 6
                && neigh.iter().all(|&n| {
                    g.vertex_type(n) == VType::Z
                        && g.phase(n).is_t()
                        && g.edge_type(v, n) == EType::H
                })
            {
                if let Some(this_ind) = preferred_order.iter().position(|&r| r == neigh.len()) {
                    match index {
                        Some(ind) if this_ind < ind => {
                            res = vec![v];
                            res.append(&mut neigh);
                            index = Some(this_ind);
                        }
                        None => {
                            res = vec![v];
                            res.append(&mut neigh);
                            index = Some(this_ind);
                        }
                        _ => (),
                    }
                }
                if index == Some(0) {
                    break;
                }
            }
        }
    }
    res
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum SimpFunc {
    FullSimp,
    CliffordSimp,
    NoSimp,
}
use SimpFunc::*;

pub enum Decomp {
    CatDecomp(Vec<usize>),
    Magic5FromCat(Vec<usize>),
    TDecomp(Vec<usize>),
    BssDecomp(Vec<usize>),
    SymDecomp(Vec<usize>),
    SingleDecomp(Vec<usize>),
}
use Decomp::*;

#[derive(Clone, Debug, derive_more::Display)]
pub enum Driver {
    BssTOnly(bool),
    BssWithCats(bool),
}
use Driver::*;

impl Driver {
    fn choose_decomp(&self, g: &impl GraphLike) -> Decomp {
        match self {
            BssTOnly(random_t) => {
                let ts = if *random_t {
                    random_ts(g, &mut thread_rng())
                } else {
                    first_ts(g)
                };
                TDecomp(ts)
            }
            BssWithCats(random_t) => {
                let cat_nodes = cat_ts(g);
                if cat_nodes.len() > 3 {
                    // println!("using cat!");
                    CatDecomp(cat_nodes)
                } else {
                    let ts = if *random_t {
                        random_ts(g, &mut thread_rng())
                    } else {
                        first_ts(g)
                    };
                    if ts.len() >= 5 {
                        // println!("using M5!");
                        Magic5FromCat(ts[0..5].to_vec())
                    } else {
                        // println!("using Ts");
                        TDecomp(ts)
                    }
                }
            }
        }
    }
}

fn replace_cat6_0<G: GraphLike>(g: &G, verts: &[V]) -> G {
    let mut g = g.clone();
    g.scalar_mut().mul_sqrt2_pow(-2);
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational64::new(-1, 4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational64::new(-1, 2));
    g
}

fn replace_cat6_1<G: GraphLike>(g: &G, verts: &[V]) -> G {
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(-1, [-1, 0, 1, 0]);
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational64::new(-1, 4));
    }
    g
}

fn replace_cat6_2<G: GraphLike>(g: &G, verts: &[V]) -> G {
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(7, [0, -1, 0, 0]);
    for i in 1..verts.len() {
        g.add_to_phase(verts[i], Rational64::new(-1, 4));
        for j in i + 1..verts.len() {
            g.add_edge_smart(verts[i], verts[j], EType::H);
        }
    }
    g
}

fn replace_magic5_0<G: GraphLike>(g: &G, verts: &[V]) -> G {
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(1, [1, 0, 0, 0]);
    for &v in verts {
        g.add_to_phase(v, Rational64::new(-1, 4));
        g.add_edge_smart(v, verts[0], EType::N);
    }
    g.add_to_phase(verts[0], Rational64::new(-3, 4));
    g
}

fn replace_magic5_1<G: GraphLike>(g: &G, verts: &[V]) -> G {
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(1, [-1, 0, 1, 0]);
    let p = g.add_vertex(VType::Z);
    for &v in verts {
        g.add_to_phase(v, Rational64::new(-1, 4));
        g.add_edge_with_type(v, p, EType::H);
    }
    let w = g.add_vertex_with_phase(VType::Z, Rational64::new(-1, 4));
    g.add_edge_with_type(w, p, EType::H);
    g
}

fn replace_magic5_2<G: GraphLike>(g: &G, verts: &[V]) -> G {
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(9, [0, -1, 0, 0]);
    let p = g.add_vertex(VType::Z);
    let w = g.add_vertex_with_phase(VType::Z, Rational64::new(-1, 4));
    g.add_edge_with_type(p, w, EType::H);
    for i in 0..verts.len() {
        g.add_to_phase(verts[i], Rational64::new(-1, 4));
        g.add_edge_with_type(verts[i], p, EType::H);
        g.add_edge_with_type(verts[i], w, EType::H);
        for j in i + 1..verts.len() {
            g.add_edge_smart(verts[i], verts[j], EType::H);
        }
    }
    g
}

fn replace_cat4_0<G: GraphLike>(g: &G, verts: &[V]) -> G {
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(0, [0, 0, 1, 0]);
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational64::new(-1, 4));
    }
    g
}

fn replace_cat4_1<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // same as replace_cat6_0, only with a different scalar
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(-1, [1, 0, -1, 0]);
    for &v in &verts[1..] {
        g.add_to_phase(v, Rational64::new(-1, 4));
        g.set_edge_type(v, verts[0], EType::N);
    }
    g.set_phase(verts[0], Rational64::new(-1, 2));
    g
}

fn replace_b60<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_b60");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(-2, [-1, 0, 1, 1]);
    for &v in &verts[0..6] {
        g.add_to_phase(v, Rational64::new(-1, 4));
    }
    g
}

fn replace_b66<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_b66");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(-2, [-1, 0, 1, -1]);
    for &v in verts {
        g.add_to_phase(v, Rational64::new(3, 4));
    }
    g
}

fn replace_e6<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_e6");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(1, [0, -1, 0, 0]);

    let w = g.add_vertex_with_phase(VType::Z, Rational64::one());
    for &v in verts {
        g.add_to_phase(v, Rational64::new(1, 4));
        g.add_edge_with_type(v, w, EType::H);
    }

    g
}

fn replace_o6<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_o6");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(1, [-1, 0, -1, 0]);

    let w = g.add_vertex(VType::Z);
    for &v in verts {
        g.add_to_phase(v, Rational64::new(1, 4));
        g.add_edge_with_type(v, w, EType::H);
    }

    g
}

fn replace_k6<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_k6");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(1, [1, 0, 0, 0]);

    let w = g.add_vertex_with_phase(VType::Z, Rational64::new(-1, 2));
    for &v in verts {
        g.add_to_phase(v, Rational64::new(-1, 4));
        g.add_edge_with_type(v, w, EType::N);
    }

    g
}

fn replace_phi1<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_phi1");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(3, [1, 0, 1, 0]);

    let mut ws = vec![];
    for i in 0..5 {
        let w = g.add_vertex(VType::Z);
        ws.push(w);
        g.add_edge_with_type(verts[i], ws[i], EType::H);
        g.add_edge_with_type(ws[i], verts[5], EType::H);
        g.add_to_phase(verts[i], Rational64::new(-1, 4));
    }

    g.add_to_phase(verts[5], Rational64::new(3, 4));

    g.add_edge_with_type(ws[0], ws[2], EType::H);
    g.add_edge_with_type(ws[0], ws[3], EType::H);
    g.add_edge_with_type(ws[1], ws[3], EType::H);
    g.add_edge_with_type(ws[1], ws[4], EType::H);
    g.add_edge_with_type(ws[2], ws[4], EType::H);

    g
}

fn replace_phi2<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // print!("replace_phi2 -> ");
    replace_phi1(
        g,
        &[verts[0], verts[1], verts[3], verts[4], verts[5], verts[2]],
    )
}

fn replace_bell_s<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_bell_s");
    let mut g = g.clone();
    g.add_edge_smart(verts[0], verts[1], EType::N);
    g.add_to_phase(verts[0], Rational64::new(-1, 4));
    g.add_to_phase(verts[1], Rational64::new(1, 4));

    g
}

fn replace_epr<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_epr");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::from_phase(Rational64::new(1, 4));
    let w = g.add_vertex_with_phase(VType::Z, Rational64::one());
    for &v in verts {
        g.add_edge_with_type(v, w, EType::H);
        g.add_to_phase(v, Rational64::new(-1, 4));
    }

    g
}

fn replace_t0<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_t0");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(-1, [0, 1, 0, -1]);
    let w = g.add_vertex(VType::Z);
    g.add_edge_with_type(verts[0], w, EType::H);
    g.add_to_phase(verts[0], Rational64::new(-1, 4));
    g
}

fn replace_t1<G: GraphLike>(g: &G, verts: &[V]) -> G {
    // println!("replace_t1");
    let mut g = g.clone();
    *g.scalar_mut() *= FScalar::dyadic(-1, [1, 0, 1, 0]);
    let w = g.add_vertex_with_phase(VType::Z, Rational64::one());
    g.add_edge_with_type(verts[0], w, EType::H);
    g.add_to_phase(verts[0], Rational64::new(-1, 4));
    g
}

fn apply_ts_decomp<G: GraphLike>(g: &G, ts: &[V]) -> Vec<G> {
    if ts.len() == 6 {
        apply_bss_decomp(g, ts)
    } else if ts.len() >= 2 {
        apply_sym_decomp(g, &ts[0..2])
    } else if !ts.is_empty() {
        apply_single_decomp(g, ts)
    } else {
        panic!("No ts!")
    }
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
fn apply_bss_decomp<G: GraphLike>(g: &G, verts: &[V]) -> Vec<G> {
    vec![
        replace_b60(g, verts),
        replace_b66(g, verts),
        replace_e6(g, verts),
        replace_o6(g, verts),
        replace_k6(g, verts),
        replace_phi1(g, verts),
        replace_phi2(g, verts),
    ]
}

/// Perform a decomposition of 2 T gates in the symmetric 2-qubit
/// space spanned by stabilisers
fn apply_sym_decomp<G: GraphLike>(g: &G, verts: &[V]) -> Vec<G> {
    vec![replace_bell_s(g, verts), replace_epr(g, verts)]
}

/// Replace a single T gate with its decomposition
fn apply_single_decomp<G: GraphLike>(g: &G, verts: &[V]) -> Vec<G> {
    vec![replace_t0(g, verts), replace_t1(g, verts)]
}

/// Perform a decomposition of 5 T-spiders, with one remaining
fn apply_magic5_from_cat_decomp<G: GraphLike>(g: &G, verts: &[V]) -> Vec<G> {
    //println!("magic5");
    vec![
        replace_magic5_0(g, verts),
        replace_magic5_1(g, verts),
        replace_magic5_2(g, verts),
    ]
}

/// Perform a decomposition of cat states
fn apply_cat_decomp<G: GraphLike>(g: &G, verts: &[V]) -> Vec<G> {
    println!("{:?}", verts);
    // verts[0] is a 0- or pi-spider, linked to all and only to vs in verts[1..] which are T-spiders
    let mut g = g.clone(); // that is annoying ...
    let mut verts = Vec::from(verts);
    let mut num_verts = verts.len() - 1;
    if g.phase(verts[0]).is_one() {
        g.set_phase(verts[0], Rational64::new(0, 1));
        let mut neigh = g.neighbor_vec(verts[1]);
        neigh.retain(|&x| x != verts[0]);
        for &v in &neigh {
            g.add_to_phase(v, Rational64::new(1, 1));
        }
        let tmp = g.phase(verts[1]);
        g.scalar_mut().mul_phase(tmp);
        g.set_phase(verts[1], g.phase(verts[1]) * -1);
    }
    if num_verts == 3 || num_verts == 5 {
        let w = g.add_vertex(VType::Z);
        let v = g.add_vertex(VType::Z);
        g.add_edge_with_type(v, w, EType::H);
        g.add_edge_with_type(v, verts[0], EType::H);
        verts.push(v);
        num_verts += 1;
    }
    if num_verts == 6 {
        vec![
            replace_cat6_0(&g, &verts),
            replace_cat6_1(&g, &verts),
            replace_cat6_2(&g, &verts),
        ]
    } else if num_verts == 4 {
        vec![replace_cat4_0(&g, &verts), replace_cat4_1(&g, &verts)]
    } else {
        panic!("this shouldn't be printed")
    }
}

#[derive(Clone)]
enum ComputationNode<G: GraphLike> {
    Graph(G),
    Scalar(FScalar),
    Prod(Vec<ComputationNode<G>>),
    Sum(Vec<ComputationNode<G>>),
    None,
}

fn calc_max_terms(node: &ComputationNode<impl GraphLike>) -> f64 {
    match node {
        ComputationNode::None => 0f64,
        ComputationNode::Scalar(_) => 1f64,
        ComputationNode::Graph(g) => terms_for_tcount(g.tcount()),
        ComputationNode::Prod(terms) => terms.iter().map(|node| calc_max_terms(node)).sum(),
        ComputationNode::Sum(terms) => terms.iter().map(|node| calc_max_terms(node)).sum(),
    }
}

/// Store the (partial) decomposition of a graph into stabilisers
#[derive(Clone)]
pub struct Decomposer<G: GraphLike> {
    pub done: Vec<G>,
    pub nterms: usize,
    result: ComputationNode<G>,
    simp_func: SimpFunc,
    driver: Driver,
    split_graph_components: bool,
    save: bool, // save graphs on 'done' stack
}

impl<G: GraphLike> Decomposer<G> {
    pub fn empty() -> Decomposer<G> {
        Decomposer {
            result: ComputationNode::None,
            done: vec![],
            nterms: 0,
            simp_func: NoSimp,
            driver: BssTOnly(false),
            split_graph_components: false,
            save: false,
        }
    }

    pub fn new(g: &G) -> Decomposer<G> {
        Decomposer {
            result: ComputationNode::Graph(g.clone()),
            done: vec![],
            nterms: 0,
            simp_func: NoSimp,
            driver: BssTOnly(false),
            split_graph_components: false,
            save: false,
        }
    }

    pub fn scalar(&self) -> FScalar {
        match self.result {
            ComputationNode::Scalar(scalar) => scalar,
            ComputationNode::None => panic!("Not yet initialised!"),
            ComputationNode::Graph(_) => panic!("Not yet computed!"),
            ComputationNode::Prod(_) => panic!("Not yet computed!"), //TODO,
            ComputationNode::Sum(_) => panic!("Not yet computed!"),
        }
    }

    pub fn with_simp(&mut self, f: SimpFunc) -> &mut Self {
        self.simp_func = f;
        self
    }
    pub fn with_full_simp(&mut self) -> &mut Self {
        self.with_simp(FullSimp)
    }
    pub fn with_clifford_simp(&mut self) -> &mut Self {
        self.with_simp(CliffordSimp)
    }

    pub fn with_split_graphs_components(&mut self, b: bool) -> &mut Self {
        self.split_graph_components = b;
        self
    }

    pub fn with_save(&mut self, b: bool) -> &mut Self {
        self.save = b;
        self
    }

    pub fn with_driver(&mut self, d: Driver) -> &mut Self {
        self.driver = d;
        self
    }

    /// Computes the maximum number of terms that this decomposer will produce
    pub fn max_terms(&self) -> f64 {
        calc_max_terms(&self.result)
    }

    pub fn set_target(&mut self, g: G) -> &mut Self {
        self.result = ComputationNode::Graph(g.clone());
        self
    }

    #[deprecated = "This function is outdated and may not work as expected, please use with_driver instead"]
    pub fn use_cats(&mut self, b: bool) -> &mut Self {
        if b {
            self.driver = BssWithCats(false);
        } else {
            self.driver = BssTOnly(false);
        }
        self
    }

    #[deprecated = "This function is outdated and may not work as expected, please use with_driver instead"]
    pub fn random_t(&mut self, b: bool) -> &mut Self {
        self.driver = BssTOnly(b);
        self
    }

    #[deprecated = "This function is outdated, please use decompose_parallel instead, depth is not used"]
    pub fn decomp_parallel(mut self, _depth: usize) -> Self {
        self.decompose_parallel();
        self
    }

    pub fn decomp_until_depth(&mut self, depth: i64) -> &mut Self {
        self.result = self.decompose_graph(self.result.clone(), false, 0, depth, false);
        self
    }

    #[deprecated = "This function is outdated and will only decompose the inital graph for one step"]
    pub fn decomp_top(&mut self) -> &mut Self {
        self.decomp_until_depth(1)
    }

    /// Decompose until there are no T gates left
    pub fn decompose(&mut self) -> &mut Self {
        self.result = self.decompose_graph(self.result.clone(), false, 0, -1, true);
        self
    }

    pub fn decompose_parallel(&mut self) -> &mut Self {
        self.result = self.decompose_graph(self.result.clone(), true, 0, -1, true);
        self
    }

    #[deprecated = "This function is outdated, please use decompose instead"]
    pub fn decomp_all(&mut self) -> &mut Self {
        self.decompose()
    }

    fn node_to_scalar(&mut self, node: ComputationNode<G>) -> FScalar {
        if let ComputationNode::Scalar(scalar) = node {
            scalar
        } else {
            panic!("Not a Scalar!")
        }
    }

    fn decompose_graph(
        &mut self,
        node: ComputationNode<G>,
        parallel: bool,
        current_depth: i64,
        target_depth: i64,
        reduce_computation: bool,
    ) -> ComputationNode<G> {
        if reduce_computation && (target_depth != -1) {
            panic!("If reducing the computation the target_depth has to be -1")
        }
        match node {
            ComputationNode::None => panic!("Not yet initialised"),
            ComputationNode::Scalar(_) => node,
            ComputationNode::Sum(terms) => {
                let results: Vec<_> = terms
                    .into_iter()
                    .map(|term| {
                        self.decompose_graph(term, parallel, current_depth + 1, target_depth, true)
                    })
                    .collect();
                if reduce_computation {
                    ComputationNode::Scalar(
                        results
                            .into_iter()
                            .map(|node| self.node_to_scalar(node))
                            .sum(),
                    )
                } else {
                    ComputationNode::Sum(results)
                }
            }
            ComputationNode::Prod(terms) => {
                if reduce_computation {
                    let results: Vec<_> = terms
                        .into_iter()
                        .map(|term| {
                            self.decompose_graph(
                                term,
                                parallel,
                                current_depth + 1,
                                target_depth,
                                true,
                            )
                        })
                        .collect();
                    if reduce_computation {
                        ComputationNode::Scalar(
                            results
                                .into_iter()
                                .map(|node| self.node_to_scalar(node))
                                .product(),
                        )
                    } else {
                        ComputationNode::Prod(results)
                    }
                } else {
                    ComputationNode::Prod(
                        terms
                            .into_iter()
                            .map(|term| {
                                self.decompose_graph(
                                    term,
                                    parallel,
                                    current_depth + 1,
                                    target_depth,
                                    true,
                                )
                            })
                            .collect(),
                    )
                }
            }
            ComputationNode::Graph(mut g) => {
                if current_depth == target_depth {
                    ComputationNode::Graph(g)
                } else {
                    if self.split_graph_components {
                        let components = g.component_vertices();
                        if components.len() > 1 {
                            // println!("Number of components {}", components.len());
                            let subgraphs: Vec<G> = components
                                .into_iter()
                                .map(|component| {
                                    g.subgraph_from_vertices(component.into_iter().collect())
                                })
                                .collect();
                            let terms_vec: Vec<ComputationNode<G>> = if parallel {
                                subgraphs
                                    .into_par_iter()
                                    .map(|term| {
                                        let mut d = self.clone();
                                        d.decompose_graph(
                                            ComputationNode::Graph(term),
                                            parallel,
                                            current_depth + 1,
                                            target_depth,
                                            reduce_computation,
                                        )
                                    })
                                    .collect()
                            } else {
                                subgraphs
                                    .into_iter()
                                    .map(|term| {
                                        self.decompose_graph(
                                            ComputationNode::Graph(term),
                                            parallel,
                                            current_depth + 1,
                                            target_depth,
                                            reduce_computation,
                                        )
                                    })
                                    .collect()
                            };
                            if reduce_computation {
                                // println!("Component Scalars {:?}", (terms_vec.clone().into_iter().map(|node| self.node_to_scalar(node))).collect_vec());
                                // let prod: FScalar = terms_vec.clone()
                                // .into_iter()
                                // .map(|node| self.node_to_scalar(node))
                                // .product();
                                // println!("Product; {}", prod);
                                return ComputationNode::Scalar(
                                    terms_vec
                                        .into_iter()
                                        .map(|node| self.node_to_scalar(node))
                                        .product(),
                                );
                            } else {
                                return ComputationNode::Prod(terms_vec);
                            }
                        }
                    };
                    match self.simp_func {
                        FullSimp => {
                            crate::simplify::full_simp(&mut g);
                        }
                        CliffordSimp => {
                            crate::simplify::clifford_simp(&mut g);
                        }
                        _ => {}
                    }
                    //check if clifford
                    if g.tcount() == 0 {
                        crate::simplify::full_simp(&mut g);
                        self.nterms += 1;
                        if g.inputs().is_empty() && g.outputs().is_empty() && g.num_vertices() != 0
                        {
                            println!("{}", g.to_dot());
                            panic!("WARNING: graph was not fully reduced");
                        }
                        if self.save {
                            self.done.push(g.clone());
                        }
                        return ComputationNode::Scalar(*g.scalar());
                    }
                    let decomp = self.driver.choose_decomp(&g);
                    let terms = match decomp {
                        Magic5FromCat(vertices) => {
                            apply_magic5_from_cat_decomp(&g, &vertices[0..5])
                        }
                        TDecomp(vertices) => apply_ts_decomp(&g, &vertices),
                        CatDecomp(vertices) => apply_cat_decomp(&g, &vertices),
                        BssDecomp(vertices) => apply_bss_decomp(&g, &vertices),
                        SymDecomp(vertices) => apply_sym_decomp(&g, &vertices),
                        SingleDecomp(vertices) => apply_single_decomp(&g, &vertices),
                    };
                    let terms_vec: Vec<ComputationNode<G>> = if parallel {
                        terms
                            .into_par_iter()
                            .map(|term| {
                                let mut d = self.clone();
                                d.decompose_graph(
                                    ComputationNode::Graph(term),
                                    parallel,
                                    current_depth + 1,
                                    target_depth,
                                    reduce_computation,
                                )
                            })
                            .collect()
                    } else {
                        terms
                            .into_iter()
                            .map(|term| {
                                self.decompose_graph(
                                    ComputationNode::Graph(term),
                                    parallel,
                                    current_depth + 1,
                                    target_depth,
                                    reduce_computation,
                                )
                            })
                            .collect()
                    };
                    if reduce_computation {
                        ComputationNode::Scalar(
                            terms_vec
                                .into_iter()
                                .map(|node| self.node_to_scalar(node))
                                .sum(),
                        )
                    } else {
                        ComputationNode::Sum(terms_vec)
                    }
                }
            }
        }
    }
}
#[cfg(test)]
mod tests {
    // use num::rational::Ratio;

    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use super::*;
    use crate::tensor::*;
    use crate::vec_graph::Graph;
    // use itertools::Itertools;

    // Helper function to create a simple graph with T gates (no outputs)
    fn create_t_graph(n: usize) -> Graph {
        let mut g = Graph::new();
        for _ in 0..n {
            g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
        }
        g
    }

    fn create_graph(n: usize) -> Graph {
        let mut g = Graph::new();
        let mut rng = StdRng::seed_from_u64(42);

        // Create a random circuit with T gates and Hadamard edges
        for i in 0..n {
            let v = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
            // Randomly connect to previous vertices with 50% probability
            for j in 0..i {
                if rng.gen_bool(0.5) {
                    g.add_edge_with_type(v, j, EType::H);
                }
            }
        }
        g
    }

    // Helper function to create a cat state graph (no outputs)
    fn create_cat_graph(n: usize, pauli_phase: Rational64) -> Graph {
        let mut g = Graph::new();
        let z = g.add_vertex_with_phase(VType::Z, pauli_phase);
        for _ in 0..n {
            let t = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
            g.add_edge_with_type(z, t, EType::H);
        }
        g
    }

    // Test individual replacement functions by checking tensor equality
    #[test]
    fn test_single_t_replacements() {
        let mut g = Graph::new();
        let v = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
        let original_scalar = g.to_tensorf()[[]];

        // Test replace_t0
        let g0 = replace_t0(&g, &[v]);
        let s0 = g0.to_tensorf()[[]];

        // Test replace_t1
        let g1 = replace_t1(&g, &[v]);
        let s1 = g1.to_tensorf()[[]];

        // The sum should equal the original
        assert_eq!(original_scalar, s0 + s1);
    }

    #[test]
    fn test_symmetric_decomp() {
        let g = create_t_graph(2);
        let original_scalar = g.to_tensorf()[[]];
        let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();

        let decomp = apply_sym_decomp(&g, &ts[0..2]);
        assert_eq!(decomp.len(), 2);

        let sum: FScalar = decomp.iter().map(|g| g.to_tensorf()[[]]).sum();

        assert_eq!(original_scalar, sum);
    }

    #[test]
    fn test_bss_decomp() {
        let g = create_t_graph(6);
        let original_scalar = g.to_tensorf()[[]];
        let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();

        let decomp = apply_bss_decomp(&g, &ts[0..6]);
        assert_eq!(decomp.len(), 7);

        let sum: FScalar = decomp.iter().map(|g| g.to_tensorf()[[]]).sum();

        assert_eq!(original_scalar, sum);
    }

    #[test]
    fn test_cat3_decomp() {
        let g = create_cat_graph(3, Rational64::new(1, 1));
        let original_scalar = g.to_tensorf()[[]];
        let verts = cat_ts(&g);

        assert_eq!(verts.len(), 4); // 1 Pauli + 3 T

        let decomp = apply_cat_decomp(&g, &verts);
        assert_eq!(decomp.len(), 2);

        let sum: FScalar = decomp.iter().map(|g| g.to_tensorf()[[]]).sum();

        assert_eq!(original_scalar, sum);
    }
    #[test]
    fn test_cat4_decomp() {
        let g = create_cat_graph(4, Rational64::new(1, 1));
        let original_scalar = g.to_tensorf()[[]];
        let verts = cat_ts(&g);

        assert_eq!(verts.len(), 5); // 1 Pauli + 4 T

        let decomp = apply_cat_decomp(&g, &verts);
        assert_eq!(decomp.len(), 2);

        let sum: FScalar = decomp.iter().map(|g| g.to_tensorf()[[]]).sum();

        assert_eq!(original_scalar, sum);
    }

    #[test]
    fn test_cat5_decomp() {
        let g = create_cat_graph(5, Rational64::new(1, 1));
        let original_scalar = g.to_tensorf()[[]];
        let verts = cat_ts(&g);

        assert_eq!(verts.len(), 6); // 1 Pauli + 5 T

        let decomp = apply_cat_decomp(&g, &verts);
        assert_eq!(decomp.len(), 3);

        let sum: FScalar = decomp.iter().map(|g| g.to_tensorf()[[]]).sum();

        assert_eq!(original_scalar, sum);
    }

    #[test]
    fn test_cat6_decomp() {
        let g = create_cat_graph(6, Rational64::new(1, 1));
        let original_scalar = g.to_tensorf()[[]];
        let verts = cat_ts(&g);

        assert_eq!(verts.len(), 7); // 1 Pauli + 6 T

        let decomp = apply_cat_decomp(&g, &verts);
        assert_eq!(decomp.len(), 3);

        let sum: FScalar = decomp.iter().map(|g| g.to_tensorf()[[]]).sum();

        assert_eq!(original_scalar, sum);
    }
    // #[test]
    // fn test_cat_decomp_in_graph() {
    //     for size in 100..110{
    //         let mut g = create_graph(size);
    //         crate::simplify::full_simp(&mut g);
    //         // let original_scalar = g.to_tensorf()[[]];
    //         let verts = cat_ts(&g);
    //         if verts.len() < 3 {
    //             continue;
    //         }
    //         println!("Doing");
    //         let decomp = apply_cat_decomp(&g, &verts);
    //         let sum: FScalar = decomp.iter()
    //             .map(|g| g.to_tensorf()[[]])
    //             .sum();
    //         // assert_eq!(original_scalar, sum);
    //     }
    // }

    #[test]
    fn test_magic5_decomp() {
        let g = create_t_graph(5);
        let original_scalar = g.to_tensorf()[[]];
        let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();

        let decomp = apply_magic5_from_cat_decomp(&g, &ts[0..5]);
        assert_eq!(decomp.len(), 3);

        let sum: FScalar = decomp.iter().map(|g| g.to_tensorf()[[]]).sum();

        assert_eq!(original_scalar, sum);
    }

    // Test all configurations of the decomposer
    #[test]
    fn test_decomposer_all_configs() {
        // Test configurations
        let simp_funcs = vec![NoSimp, CliffordSimp, FullSimp];
        let drivers = vec![
            BssTOnly(false),
            BssTOnly(true),
            BssWithCats(false),
            BssWithCats(true),
        ];
        let split_components = vec![false, true];
        let parallel_modes = vec![false, true];
        let graph_generators = vec![create_t_graph, create_graph];

        // Test with different T-gate counts
        for graph_generator in graph_generators {
            for size in 0..=10 {
                let g = graph_generator(size);
                let expected_scalar = g.to_tensorf()[[]];

                for simp in &simp_funcs {
                    for driver in &drivers {
                        for &split in &split_components {
                            for &parallel in &parallel_modes {
                                let mut d = Decomposer::new(&g);
                                d.with_simp(*simp).with_driver(driver.clone());
                                d.split_graph_components = split;
                                // println!("X");
                                if parallel {
                                    d.decompose_parallel();
                                } else {
                                    d.decompose();
                                }

                                let result_scalar = d.scalar();
                                assert_eq!(
                                        expected_scalar, result_scalar,
                                        "Failed for t_count={}, simp={:?}, driver={:?}, split={}, parallel={}",
                                        size, simp, driver, split, parallel
                                    );
                            }
                        }
                    }
                }
            }
        }
    }

    // Test decomposer with cat states
    #[test]
    fn test_decomposer_with_cats() {
        for cat_size in [3, 4, 5, 6] {
            let g = create_cat_graph(cat_size, Rational64::new(0, 1));
            let expected_scalar = g.to_tensorf()[[]];

            let mut d = Decomposer::new(&g);
            d.with_full_simp()
                .with_driver(BssWithCats(false))
                .decompose();

            assert_eq!(
                expected_scalar,
                d.scalar(),
                "Cat state decomposition failed for size {}",
                cat_size
            );
        }
    }

    // Test mixed graphs with both regular T gates and cat states
    #[test]
    fn test_mixed_graph() {
        let mut g = Graph::new();

        // Add some regular T gates
        for _ in 0..3 {
            g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
        }

        // Add a cat state
        let z = g.add_vertex(VType::Z);
        for _ in 0..4 {
            let t = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
            g.add_edge_with_type(z, t, EType::H);
        }

        let expected_scalar = g.to_tensorf()[[]];

        // Test with cat-aware driver
        let mut d = Decomposer::new(&g);
        d.with_full_simp()
            .with_driver(BssWithCats(false))
            .decompose();

        assert_eq!(expected_scalar, d.scalar());
    }

    // Test edge cases
    #[test]
    fn test_empty_graph() {
        let g = Graph::new();
        let mut d = Decomposer::new(&g);
        d.decompose();
        assert_eq!(FScalar::one(), d.scalar());
    }

    #[test]
    fn test_clifford_only_graph() {
        let mut g = Graph::new();
        // Create a simple Clifford circuit that reduces to a scalar
        let v1 = g.add_vertex(VType::Z);
        let v2 = g.add_vertex(VType::Z);
        let v3 = g.add_vertex(VType::Z);
        g.add_edge_with_type(v1, v2, EType::H);
        g.add_edge_with_type(v2, v3, EType::H);
        g.add_edge_with_type(v1, v3, EType::H);

        let expected_scalar = g.to_tensorf()[[]];

        let mut d = Decomposer::new(&g);
        d.with_simp(CliffordSimp).decompose();

        // The decomposer should handle Clifford-only graphs
        assert_eq!(expected_scalar, d.scalar());
    }

    // Test specific replacement functions for correctness
    #[test]
    fn test_replacement_scalars() {
        // Test that each replacement has the correct scalar factor
        let g = create_t_graph(6);
        let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();

        // Check BSS replacement scalars
        let b60_g = replace_b60(&g, &ts);
        assert_eq!(
            b60_g.scalar(),
            &(g.scalar() * FScalar::dyadic(-2, [-1, 0, 1, 1]))
        );

        let b66_g = replace_b66(&g, &ts);
        assert_eq!(
            b66_g.scalar(),
            &(g.scalar() * FScalar::dyadic(-2, [-1, 0, 1, -1]))
        );

        let e6_g = replace_e6(&g, &ts);
        assert_eq!(
            e6_g.scalar(),
            &(g.scalar() * FScalar::dyadic(1, [0, -1, 0, 0]))
        );

        let o6_g = replace_o6(&g, &ts);
        assert_eq!(
            o6_g.scalar(),
            &(g.scalar() * FScalar::dyadic(1, [-1, 0, -1, 0]))
        );

        let k6_g = replace_k6(&g, &ts);
        assert_eq!(
            k6_g.scalar(),
            &(g.scalar() * FScalar::dyadic(1, [1, 0, 0, 0]))
        );

        let phi1_g = replace_phi1(&g, &ts);
        assert_eq!(
            phi1_g.scalar(),
            &(g.scalar() * FScalar::dyadic(3, [1, 0, 1, 0]))
        );
    }

    // Test phase handling in cat decomposition
    #[test]
    fn test_cat_with_pi_phase() {
        let g = create_cat_graph(4, Rational64::new(1, 1)); // Pi phase
        let expected_scalar = g.to_tensorf()[[]];

        let mut d = Decomposer::new(&g);
        d.with_driver(BssWithCats(false)).decompose();

        assert_eq!(expected_scalar, d.scalar());
    }

    // Test that all replacement functions preserve the sum
    #[test]
    fn test_all_replacements_sum_preservation() {
        // Test single T replacement
        {
            let g = create_t_graph(1);
            let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();
            let replacements = apply_single_decomp(&g, &ts);
            let sum: FScalar = replacements.iter().map(|g| g.to_tensorf()[[]]).sum();
            assert_eq!(g.to_tensorf()[[]], sum);
        }

        // Test symmetric (2 T) replacement
        {
            let g = create_t_graph(2);
            let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();
            let replacements = apply_sym_decomp(&g, &ts);
            let sum: FScalar = replacements.iter().map(|g| g.to_tensorf()[[]]).sum();
            assert_eq!(g.to_tensorf()[[]], sum);
        }

        // Test BSS (6 T) replacement
        {
            let g = create_t_graph(6);
            let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();
            let replacements = apply_bss_decomp(&g, &ts);
            let sum: FScalar = replacements.iter().map(|g| g.to_tensorf()[[]]).sum();
            assert_eq!(g.to_tensorf()[[]], sum);
        }

        // Test magic5 replacement
        {
            let g = create_t_graph(5);
            let ts: Vec<_> = g.vertices().filter(|&v| g.phase(v).is_t()).collect();
            let replacements = apply_magic5_from_cat_decomp(&g, &ts);
            let sum: FScalar = replacements.iter().map(|g| g.to_tensorf()[[]]).sum();
            assert_eq!(g.to_tensorf()[[]], sum);
        }
    }

    // Test split_graph_components functionality
    #[test]
    fn test_split_components() {
        // Create a graph with two disconnected components
        let mut g = Graph::new();

        // Component 1: 3 T gates
        for _ in 0..3 {
            g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
        }

        // Component 2: 2 T gates connected
        let t1 = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
        let t2 = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
        let z = g.add_vertex(VType::Z);
        g.add_edge_with_type(t1, z, EType::H);
        g.add_edge_with_type(t2, z, EType::H);

        let expected_scalar = g.to_tensorf()[[]];

        // Test with split_graph_components = true
        let mut d = Decomposer::new(&g);
        d.with_full_simp();
        d.split_graph_components = true;
        d.decompose();

        assert_eq!(expected_scalar, d.scalar());
    }

    // Existing tests from the original code (kept for compatibility)
    #[test]
    fn bss_scalars() {
        // ... (keep existing test)
        let one = FScalar::one();
        let om = FScalar::dyadic(0, [0, 1, 0, 0]);
        let om2 = om * om;
        let om7 = FScalar::dyadic(0, [0, 0, 0, -1]);
        assert_eq!(om * om7, FScalar::one());

        let minus = FScalar::dyadic(0, [-1, 0, 0, 0]);
        let onefourth = FScalar::dyadic(-2, [1, 0, 0, 0]);
        let two = one + one;
        let sqrt2 = FScalar::sqrt2();
        let eight = two * two * two;

        let k6 = om7 * two * om;
        let phi = om7 * eight * sqrt2 * om2;
        let b60 = om7 * minus * onefourth * (one + sqrt2);
        let b66 = om7 * onefourth * (one + (minus * sqrt2));
        let o6 = om7 * minus * two * sqrt2 * om2;
        let e6 = om7 * minus * two * om2;

        assert_eq!(b60, FScalar::dyadic(-2, [-1, 0, 1, 1]));
        assert_eq!(b66, FScalar::dyadic(-2, [-1, 0, 1, -1]));
        assert_eq!(e6, FScalar::dyadic(1, [0, -1, 0, 0]));
        assert_eq!(o6, FScalar::dyadic(1, [-1, 0, -1, 0]));
        assert_eq!(k6, FScalar::dyadic(1, [1, 0, 0, 0]));
        assert_eq!(phi, FScalar::dyadic(3, [1, 0, 1, 0]));
    }

    #[test]
    fn single_scalars() {
        let s0 = FScalar::sqrt2_pow(-1);
        let s1 = FScalar::from_phase(Rational64::new(1, 4)) * s0;
        // println!("s0 = {s0:?}\ns1 = {s1:?}");
        assert_eq!(s0, FScalar::dyadic(-1, [0, 1, 0, -1]));
        assert_eq!(s1, FScalar::dyadic(-1, [1, 0, 1, 0]));
    }

    #[test]
    fn mixed_sc() {
        let mut g = Graph::new();
        for i in 0..20 {
            g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));

            for j in 0..i {
                g.add_edge_with_type(i, j, EType::H);
            }
        }

        let mut d = Decomposer::new(&g);
        d.with_full_simp();
        d.decompose();

        let sc = g.to_tensorf()[[]];
        // println!("{}", d.nterms);
        assert_eq!(sc, d.scalar());
    }

    #[test]
    fn full_simp() {
        let mut g = Graph::new();
        let mut outs = vec![];
        for _ in 0..9 {
            let v = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
            let w = g.add_vertex(VType::B);
            outs.push(w);
            g.add_edge(v, w);
        }
        g.set_outputs(outs);

        let mut d = Decomposer::new(&g);
        d.with_full_simp()
            .with_driver(Driver::BssTOnly(false))
            .with_save(true)
            .decompose();
        assert_eq!(d.done.len(), 7 * 2 * 2);
    }

    #[test]
    fn cat4() {
        let mut g = Graph::new();

        let mut outputs = vec![];
        let z = g.add_vertex(VType::Z);
        for _ in 0..4 {
            let t = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
            g.add_edge_with_type(z, t, EType::H);

            let out = g.add_vertex(VType::B);
            g.add_edge(t, out);
            outputs.push(out);
        }
        g.set_outputs(outputs);

        let mut d = Decomposer::new(&g);
        d.with_full_simp()
            .with_driver(Driver::BssWithCats(false))
            .with_save(true)
            .decompose();
        assert_eq!(d.done.len(), 2);
    }

    #[test]
    fn cat6() {
        let mut g = Graph::new();

        let mut outputs = vec![];
        let z = g.add_vertex(VType::Z);
        for _ in 0..6 {
            let t = g.add_vertex_with_phase(VType::Z, Rational64::new(1, 4));
            g.add_edge_with_type(z, t, EType::H);

            let out = g.add_vertex(VType::B);
            g.add_edge(t, out);
            outputs.push(out);
        }
        g.set_outputs(outputs);

        let mut d = Decomposer::new(&g);
        d.with_full_simp()
            .with_driver(Driver::BssWithCats(false))
            .with_save(true)
            .decompose();

        assert_eq!(d.done.len(), 3);
    }
}
