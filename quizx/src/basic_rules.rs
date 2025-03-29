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

//! # Basic ZX-calculus rules
//!
//! These rules always come in triples of functions. For a rule X,
//! there is a function `check_X(&g, ...) -> bool` which checks
//! whether a rule is applicable at the given vertex or vertices,
//! `X_unchecked(&mut g, ...)` applies the rule without doing any
//! checking, and `X(&mut g, ...) -> bool` is the composition of the
//! first two.
//!
//! Note calling `X_unchecked` is allowed to make unsound ZX-diagram
//! transformations, or even panic, if `check_X` doesn't return true.

use crate::graph::*;
use crate::phase::Phase;
use crate::scalar::*;
use num::traits::Zero;
use num::Rational64;
use rustc_hash::FxHashSet;
use std::iter::FromIterator;

/// Define a checked rule that takes 1 vertex
macro_rules! checked_rule1 {
    ( $check:ident, $unchecked:ident, $name:ident ) => {
        /// A checked implementation of the rule
        ///
        /// See e.g. [spider_fusion] for an example.
        pub fn $name(g: &mut impl GraphLike, v: V) -> bool {
            if $check(g, v) {
                $unchecked(g, v);
                true
            } else {
                false
            }
        }
    };
}

/// Define a checked rule that takes 2 vertices
macro_rules! checked_rule2 {
    ( $check:ident, $unchecked:ident, $name:ident ) => {
        /// A checked implementation of the rule
        ///
        /// See e.g. [spider_fusion] for an example.
        pub fn $name(g: &mut impl GraphLike, v0: V, v1: V) -> bool {
            if $check(g, v0, v1) {
                $unchecked(g, v0, v1);
                true
            } else {
                false
            }
        }
    };
}

/// Check [spider_fusion_unchecked] applies
///
/// Both vertices must be Z or X, have the same type, and be connected
/// by a normal (i.e. non-Hadamard) edge.
///
/// ```
/// # use quizx::graph::*;
/// # use quizx::vec_graph::Graph;
/// # use quizx::basic_rules::check_spider_fusion;
/// let mut g = Graph::new();
/// let v0 = g.add_vertex(VType::Z);
/// let v1 = g.add_vertex(VType::Z);
/// let v2 = g.add_vertex(VType::X);
/// g.add_edge(v0, v1);
/// g.add_edge(v1, v2);
///
/// assert!(check_spider_fusion(&g, v0, v1));
/// assert!(!check_spider_fusion(&g, v1, v2));
/// ```
pub fn check_spider_fusion(g: &impl GraphLike, v0: V, v1: V) -> bool {
    v0 != v1
        && g.edge_type_opt(v0, v1) == Some(EType::N)
        && ((g.vertex_type(v0) == VType::Z && g.vertex_type(v1) == VType::Z)
            || (g.vertex_type(v0) == VType::X && g.vertex_type(v1) == VType::X))
}

/// Apply spider fusion
///
/// Note the first vertex is preserved by the fusion, and the second
/// is deleted.
///
/// ```
/// # use quizx::graph::*;
/// # use quizx::tensor::ToTensor;
/// # use quizx::vec_graph::Graph;
/// # use quizx::basic_rules::spider_fusion_unchecked;
/// let mut g = Graph::new();
/// let v0 = g.add_vertex(VType::Z);
/// let v1 = g.add_vertex(VType::Z);
/// let v2 = g.add_vertex(VType::X);
/// g.add_edge(v0, v1);
/// g.add_edge(v1, v2);
///
/// let h = g.clone();
/// spider_fusion_unchecked(&mut g, v0, v1);
/// assert_eq!(g.to_tensor4(), h.to_tensor4());
///
/// let h = g.clone();
/// spider_fusion_unchecked(&mut g, v0, v2); // oops! Not a valid spider fusion
/// assert_ne!(g.to_tensor4(), h.to_tensor4());
/// ```
pub fn spider_fusion_unchecked(g: &mut impl GraphLike, v0: V, v1: V) {
    for (v, et) in Vec::from_iter(g.incident_edges(v1)) {
        if v != v0 {
            g.add_edge_smart(v0, v, et);
        }
    }

    g.add_to_phase(v0, g.phase(v1));
    g.remove_vertex(v1);
}

/// A checked implementation of the rule
///
/// Note the first vertex is preserved by the fusion, and the second
/// is deleted.
///
/// ```
/// # use quizx::graph::*;
/// # use quizx::tensor::ToTensor;
/// # use quizx::vec_graph::Graph;
/// # use quizx::basic_rules::spider_fusion;
/// let mut g = Graph::new();
/// let v0 = g.add_vertex(VType::Z);
/// let v1 = g.add_vertex(VType::Z);
/// let v2 = g.add_vertex(VType::X);
/// g.add_edge(v0, v1);
/// g.add_edge(v1, v2);
///
/// let h = g.clone();
/// let success = spider_fusion(&mut g, v0, v1);
/// assert!(success);
/// assert_eq!(g.to_tensor4(), h.to_tensor4());
///
/// let h = g.clone();
/// let success = spider_fusion(&mut g, v0, v2); // should fail
/// assert!(!success);
/// assert_eq!(g, h); // g is unchanged
/// ```
pub fn spider_fusion(g: &mut impl GraphLike, v0: V, v1: V) -> bool {
    if check_spider_fusion(g, v0, v1) {
        spider_fusion_unchecked(g, v0, v1);
        true
    } else {
        false
    }
}

/// Check [pi_copy_unchecked] applies
pub fn check_pi_copy(g: &impl GraphLike, v: V) -> bool {
    let vt = g.vertex_type(v);
    // Find the opposite color of this spider, or fail
    // if it is not a colored spider.
    let ovt = match vt {
        VType::Z => VType::X,
        VType::X => VType::Z,
        _ => return false,
    };

    // No pi-copy on empty spiders.
    if g.degree(v) == 0 {
        return false;
    }

    for neighbor in g.neighbors(v) {
        // Every neighbor must be the same color and connected by
        // a hadamard edge or the opposite color and connected by a normal edge.
        match g.edge_type(v, neighbor) {
            EType::N if g.vertex_type(neighbor) != ovt => return false,
            EType::H if g.vertex_type(neighbor) != vt => return false,
            _ => (),
        }
    }

    true
}

/// Apply a pi-copy
///
/// This will flip the phase of a spider by adding a pi phase
/// to all its neighbours, assuming that they are either the
/// same color and connected by a Hadamard edge, or the opposite
/// color and connected by a normal edge. In particular, this means
/// that this rule can always be applied when the graph is in gh form.
pub fn pi_copy_unchecked(g: &mut impl GraphLike, v: V) {
    // Flip the phase of this node
    let phase = g.phase(v);
    g.scalar_mut().mul_phase(phase);
    g.set_phase(v, -phase);

    // Push a pi to all the surrounding nodes
    for neighbor in g.neighbor_vec(v) {
        g.add_to_phase(neighbor, 1);
    }
}

checked_rule1!(check_pi_copy, pi_copy_unchecked, pi_copy);

/// Check [remove_id_unchecked] applies
pub fn check_remove_id(g: &impl GraphLike, v: V) -> bool {
    let vt = g.vertex_type(v);

    (vt == VType::Z || vt == VType::X) && g.phase(v).is_zero() && g.degree(v) == 2
}

/// Remove an arity-2 spider with phase 0
///
/// Removes the spider and connects its two neighbors. The type
/// of the resulting edge is the parity of the types of
/// original 2 edges, namely: {N,N} -> N, {N,H} -> H, and
/// {H, H} -> N.
pub fn remove_id_unchecked(g: &mut impl GraphLike, v: V) {
    let nhd: Vec<(V, EType)> = g.incident_edges(v).collect();
    let new_et = match (nhd[0].1, nhd[1].1) {
        (EType::N, EType::N) => EType::N,
        (EType::N, EType::H) => EType::H,
        (EType::H, EType::N) => EType::H,
        (EType::H, EType::H) => EType::N,
        (EType::Wio, _) | (_, EType::Wio) => unimplemented!("W nodes not supported"),
    };
    g.add_edge_smart(nhd[0].0, nhd[1].0, new_et);
    g.remove_vertex(v);
}

checked_rule1!(check_remove_id, remove_id_unchecked, remove_id);

/// Check [color_change_unchecked] applies
pub fn check_color_change(g: &impl GraphLike, v: V) -> bool {
    let vt = g.vertex_type(v);
    vt == VType::X || vt == VType::Z
}

/// Change the color of a Z or X spider
///
/// All of the neighboring edge types are toggled, i.e. N -> H,
/// H -> N.
pub fn color_change_unchecked(g: &mut impl GraphLike, v: V) {
    let vt = g.vertex_type(v);
    g.set_vertex_type(v, if vt == VType::X { VType::Z } else { VType::X });
    for w in Vec::from_iter(g.neighbors(v)) {
        g.toggle_edge_type(v, w);
    }
}

checked_rule1!(check_color_change, color_change_unchecked, color_change);

/// Check [local_comp_unchecked] applies
///
/// The vertex must be Z, have a phase pi/2 or -pi/2, and be
/// surrounded by H-edges connected to other Z spiders.
pub fn check_local_comp(g: &impl GraphLike, v: V) -> bool {
    g.vertex_type(v) == VType::Z
        && g.phase(v).is_proper_clifford()
        && g.incident_edges(v)
            .all(|(v0, et)| g.vertex_type(v0) == VType::Z && et == EType::H)
}

/// Apply a local complementation
///
/// This is the version that deletes the targeted vertex. In
/// other words, it is an N-ary generalization of the Euler
/// decomposition rule.
pub fn local_comp_unchecked(g: &mut impl GraphLike, v: V) {
    let p = g.phase(v);

    // add a totally connected graph of the nhd of v
    let ns: Vec<V> = g.neighbors(v).collect();
    for i in 0..ns.len() {
        g.add_to_phase(ns[i], -p);
        for j in (i + 1)..ns.len() {
            g.add_edge_smart(ns[i], ns[j], EType::H);
        }
    }
    g.remove_vertex(v);

    let x = ns.len() as i32;
    g.scalar_mut().mul_sqrt2_pow(((x - 1) * (x - 2)) / 2);
    g.scalar_mut()
        .mul_phase(Rational64::new(*p.to_rational().numer(), 4));
}

checked_rule1!(check_local_comp, local_comp_unchecked, local_comp);

/// Check [pivot_unchecked] applies
///
/// Both vertices must be Z, have a phase 0 or pi, and be
/// surrounded by H-edges connected to other Z spiders.
pub fn check_pivot(g: &impl GraphLike, v0: V, v1: V) -> bool {
    g.vertex_type(v0) == VType::Z
        && g.vertex_type(v1) == VType::Z
        && g.edge_type_opt(v0, v1) == Some(EType::H)
        && g.phase(v0).is_pauli()
        && g.phase(v1).is_pauli()
        && g.incident_edges(v0)
            .all(|(w, et)| g.vertex_type(w) == VType::Z && et == EType::H)
        && g.incident_edges(v1)
            .all(|(w, et)| g.vertex_type(w) == VType::Z && et == EType::H)
}

/// Apply pivoting to a pair of vertices
///
/// This is the version that deletes both vertices, so it is
/// effectively a generalised version of the strong complementarity
/// rule.
pub fn pivot_unchecked(g: &mut impl GraphLike, v0: V, v1: V) {
    let p0 = g.phase(v0);
    let p1 = g.phase(v1);

    // add a complete bipartite graph between the neighbors of v0
    // and the neighbors of v1
    let ns0: Vec<V> = g.neighbors(v0).collect();
    let ns1: Vec<V> = g.neighbors(v1).collect();
    // let mut z: i32 = 0; // the number of neighbors of v0 and v1
    for &n0 in &ns0 {
        g.add_to_phase(n0, p1);
        for &n1 in &ns1 {
            if n0 != v1 && n1 != v0 {
                // unlike PyZX, add_edge_smart handles self-loops
                g.add_edge_smart(n0, n1, EType::H);
            }
        }
    }

    for &n1 in &ns1 {
        g.add_to_phase(n1, p0);
    }

    g.remove_vertex(v0);
    g.remove_vertex(v1);

    let x = ns0.len() as i32; // the number of neighbors of v0
    let y = ns1.len() as i32; // the number of neighbors of v1
    g.scalar_mut().mul_sqrt2_pow((x - 2) * (y - 2));

    if !p0.is_zero() && !p1.is_zero() {
        g.scalar_mut().mul_phase(Rational64::new(1, 1));
    }
}

checked_rule2!(check_pivot, pivot_unchecked, pivot);

/// Insert an identity spider so v is no longer adjancent to boundary b
///
/// If b is not a boundary, this is a noop. The new vertex will be connected
/// to v by a Hadamard edge.
fn unfuse_boundary(g: &mut impl GraphLike, v: V, b: V) {
    if g.vertex_type(b) != VType::B {
        return;
    }
    let vd = VData {
        ty: VType::Z,
        phase: Phase::zero(),
        row: g.row(v),
        qubit: g.qubit(v),
    };
    let v1 = g.add_vertex_with_data(vd);
    g.add_edge_with_type(v, v1, EType::H);
    g.add_edge_with_type(v1, b, g.edge_type(v, b).opposite());
    g.remove_edge(v, b);
}

/// Unfuse a non-Pauli phase as a degree-1 phase gadget
///
/// If the vertex already has a Pauli phase, this is a noop.
fn unfuse_gadget(g: &mut impl GraphLike, v: V) {
    if g.phase(v).is_pauli() {
        return;
    }
    let vd1 = VData {
        ty: VType::Z,
        phase: Phase::zero(),
        row: g.row(v),
        qubit: -1.0,
    };

    let vd2 = VData {
        ty: VType::Z,
        phase: g.phase(v),
        row: g.row(v),
        qubit: -2.0,
    };
    let v1 = g.add_vertex_with_data(vd1);
    let v2 = g.add_vertex_with_data(vd2);
    g.set_phase(v, Phase::zero());
    g.add_edge_with_type(v, v1, EType::H);
    g.add_edge_with_type(v1, v2, EType::H);
}

/// Check gen_pivot applies
///
/// This checks that we can create a valid matching for the pivot rule
/// by gadgetizing any non-Pauli phases and inserting identities to make
/// both vertices interior. Note that repeatedly applying `gen_pivot` with
/// this checker will not always terminate.
pub fn check_gen_pivot(g: &impl GraphLike, v0: V, v1: V) -> bool {
    if v0 == v1 {
        return false;
    }
    if g.edge_type_opt(v0, v1) != Some(EType::H) {
        return false;
    }

    for &v in &[v0, v1] {
        if g.vertex_type(v) != VType::Z {
            return false;
        }
        for (w, et) in g.incident_edges(v) {
            let t = g.vertex_type(w);
            if !((t == VType::Z && et == EType::H) || t == VType::B) {
                return false;
            }
        }
    }

    true
}

// check that a vertex is interior, has phase 0 or pi, and is not
// a phase gadget
fn is_interior_pauli(g: &impl GraphLike, v: V) -> bool {
    g.phase(v).is_pauli()
        && g.neighbors(v)
            .all(|n| g.vertex_type(n) == VType::Z && g.degree(n) > 1)
}

// check that a vertex is interior, has phase 0 or pi, and is not
// a phase gadget
fn is_boundary_pauli(g: &impl GraphLike, v: V) -> bool {
    g.phase(v).is_pauli() && g.neighbors(v).any(|n| g.vertex_type(n) == VType::B)
}

/// Check gen_pivot applies and at least one vertex is interior Pauli
///
/// Unlike `check_gen_pivot`, this guarantees applying `gen_pivot` will
/// strictly decrease the number of interior Pauli vertices, hence it
/// will terminate.
pub fn check_gen_pivot_reduce(g: &impl GraphLike, v0: V, v1: V) -> bool {
    check_gen_pivot(g, v0, v1) && (is_interior_pauli(g, v0) || is_interior_pauli(g, v1))
}

pub fn check_boundary_pivot(g: &impl GraphLike, v0: V, v1: V) -> bool {
    check_gen_pivot(g, v0, v1) && is_boundary_pauli(g, v0)
}

/// Generic version of the pivot rule
///
/// This version of the pivoting rule allows either of the vertices
/// to have non-Pauli phases and/or be connected to boundaries. To handle
/// these situations, some spiders are first unfused, such that interior
/// non-Pauli spiders produce phase gadgets and boundary non-Pauli spiders
/// produce phase gates on inputs/outputs.
pub fn gen_pivot_unchecked(g: &mut impl GraphLike, v0: V, v1: V) {
    let nhd0 = g.neighbor_vec(v0);
    unfuse_gadget(g, v0);
    for &n in &nhd0 {
        unfuse_boundary(g, v0, n);
    }

    let nhd1 = g.neighbor_vec(v1);
    unfuse_gadget(g, v1);
    for &n in &nhd1 {
        unfuse_boundary(g, v1, n);
    }

    pivot_unchecked(g, v0, v1);

    // for &n in nhd0.iter().chain(nhd1.iter()) {
    //      if g.contains_vertex(n) {
    //          if remove_id(g, n) { println!("REMOVED EXTRA: {}", n); }
    //      }
    //  }
}

checked_rule2!(check_gen_pivot, gen_pivot_unchecked, gen_pivot);
checked_rule2!(check_boundary_pivot, gen_pivot_unchecked, boundary_pivot);

pub fn check_gadget_fusion(g: &impl GraphLike, v0: V, v1: V) -> bool {
    if v0 == v1 {
        return false;
    }

    let vs = [v0, v1];
    let mut nhd = [FxHashSet::default(), FxHashSet::default()];

    for i in 0..2 {
        let mut found_gphase = false;
        for (n, et) in g.incident_edges(vs[i]) {
            if et != EType::H {
                return false;
            }
            if g.vertex_type(n) != VType::Z {
                return false;
            }
            if g.degree(n) == 1 {
                if found_gphase {
                    return false;
                } else {
                    found_gphase = true;
                }
            } else {
                nhd[i].insert(n);
            }
        }

        if !found_gphase {
            return false;
        }
    }

    nhd[0] == nhd[1]
}

pub fn gadget_fusion_unchecked(g: &mut impl GraphLike, v0: V, v1: V) {
    let gphase0 = g
        .neighbors(v0)
        .find(|&n| g.degree(n) == 1)
        .expect("v0 isn't a gadget");
    let gphase1 = g
        .neighbors(v1)
        .find(|&n| g.degree(n) == 1)
        .expect("v1 isn't a gadget");
    g.add_to_phase(gphase0, g.phase(gphase1));
    g.remove_vertex(v1);
    g.remove_vertex(gphase1);

    let d = g.degree(v0) as i32;
    g.scalar_mut().mul_sqrt2_pow(2 - d);
}

checked_rule2!(check_gadget_fusion, gadget_fusion_unchecked, gadget_fusion);

pub fn check_remove_single(g: &impl GraphLike, v: V) -> bool {
    let t = g.vertex_type(v);
    g.neighbors(v).len() == 0 && (t == VType::Z || t == VType::X)
}

/// Remove an isolated Z or X vertex and add it as a global scalar
pub fn remove_single_unchecked(g: &mut impl GraphLike, v: V) {
    let p = g.phase(v);
    *g.scalar_mut() *= ScalarN::one_plus_phase(p);
    g.remove_vertex(v);
}

checked_rule1!(check_remove_single, remove_single_unchecked, remove_single);

pub fn check_remove_pair(g: &impl GraphLike, v0: V, v1: V) -> bool {
    let t0 = g.vertex_type(v0);
    let t1 = g.vertex_type(v1);

    g.neighbors(v0).len() == 1
        && g.neighbors(v1).len() == 1
        && (t0 == VType::Z || t0 == VType::X)
        && (t1 == VType::Z || t1 == VType::X)
        && g.connected(v0, v1)
}

/// Remove an isolated Z or X vertex and add it as a global scalar
pub fn remove_pair_unchecked(g: &mut impl GraphLike, v0: V, v1: V) {
    let t0 = g.vertex_type(v0);
    let t1 = g.vertex_type(v1);
    let et = g.edge_type(v0, v1);
    let p0 = g.phase(v0);
    let p1 = g.phase(v1);

    // same color
    if (t0 == t1 && et == EType::N) || (t0 != t1 && et == EType::H) {
        *g.scalar_mut() *= ScalarN::one_plus_phase(p0 + p1);
    // different colors
    } else {
        let p2 = Phase::one() + p0 + p1;
        *g.scalar_mut() *= ScalarN::one()
            + ScalarN::from_phase(p0)
            + ScalarN::from_phase(p1)
            + ScalarN::from_phase(p2);
        g.scalar_mut().mul_sqrt2_pow(-1);
    }

    g.remove_vertex(v0);
    g.remove_vertex(v1);
}

checked_rule2!(check_remove_pair, remove_pair_unchecked, remove_pair);

// Tests {{{

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::*;
    use crate::vec_graph::Graph;
    use num::Rational64;

    #[test]
    fn spider_fusion_simple() {
        let mut g = Graph::new();
        let vs = [
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];

        g.set_phase(vs[2], Rational64::new(1, 2));
        g.set_phase(vs[3], Rational64::new(1, 4));

        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[2]);
        g.add_edge(vs[2], vs[3]);
        g.add_edge(vs[3], vs[4]);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[3], vs[6]);

        g.set_inputs(vec![vs[0], vs[1]]);
        g.set_outputs(vec![vs[4], vs[5], vs[6]]);

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 6);

        let h = g.clone();
        spider_fusion(&mut g, vs[2], vs[3]);

        assert_eq!(g.num_vertices(), 6);
        assert_eq!(g.num_edges(), 5);
        assert_eq!(g.degree(vs[2]), 5);

        assert_eq!(g.to_tensor4(), h.to_tensor4());

        assert_eq!(g.phase(vs[2]), Rational64::new(3, 4).into());
    }

    #[test]
    fn spider_fusion_smart() {
        // a spider fusion that creates parallel edges, which get
        // removed by complementarity
        let mut g = Graph::new();
        let vs = [
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::X),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];

        g.set_phase(vs[2], Rational64::new(1, 2));
        g.set_phase(vs[3], Rational64::new(1, 4));

        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[2]);
        g.add_edge(vs[2], vs[3]);
        g.add_edge(vs[2], vs[4]);
        g.add_edge(vs[3], vs[4]);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[4], vs[6]);

        g.set_inputs(vec![vs[0], vs[1]]);
        g.set_outputs(vec![vs[5], vs[6]]);

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 7);
        assert_eq!(g.degree(vs[2]), 4);
        assert_eq!(g.degree(vs[4]), 3);
        assert_eq!(*g.scalar(), Scalar::one());

        let h = g.clone();
        let success = spider_fusion(&mut g, vs[2], vs[3]);
        assert!(success, "Spider fusion should match");

        assert_eq!(g.num_vertices(), 6);
        assert_eq!(g.num_edges(), 4);
        assert_eq!(g.degree(vs[2]), 3);
        assert_eq!(g.degree(vs[4]), 1);
        assert_eq!(*g.scalar(), Scalar::sqrt2_pow(-2));

        let tg = g.to_tensor4();
        let th = h.to_tensor4();
        println!("\n\ntg =\n{}", tg);
        println!("\n\nth =\n{}", th);
        assert_eq!(tg, th);

        assert_eq!(g.phase(vs[2]), Rational64::new(3, 4).into());
    }

    #[test]
    fn local_comp_1() {
        let mut g = Graph::new();
        g.add_vertex(VType::Z);
        g.set_phase(0, Rational64::new(1, 2));
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_edge_with_type(0, 1, EType::H);
        g.add_edge_with_type(0, 2, EType::H);
        g.add_edge_with_type(0, 3, EType::H);
        g.add_edge_with_type(0, 4, EType::H);

        let b1 = g.add_vertex(VType::B);
        let b2 = g.add_vertex(VType::B);
        let b3 = g.add_vertex(VType::B);
        let b4 = g.add_vertex(VType::B);
        g.add_edge_with_type(b1, 1, EType::H);
        g.add_edge_with_type(b2, 2, EType::H);
        g.add_edge_with_type(b3, 3, EType::H);
        g.add_edge_with_type(b4, 4, EType::H);
        g.set_inputs(vec![b1, b2]);
        g.set_outputs(vec![b3, b4]);

        assert_eq!(g.num_vertices(), 9);
        assert_eq!(g.num_edges(), 8);

        let h = g.clone();
        let success = local_comp(&mut g, 0);
        assert!(success, "Local comp should match");

        assert_eq!(g.num_vertices(), 8);
        assert_eq!(g.num_edges(), 10);

        let tg = g.to_tensor4();
        let th = h.to_tensor4();
        println!("\n\ntg =\n{}", tg);
        println!("\n\nth =\n{}", th);
        assert_eq!(tg, th);

        for i in 1..5 {
            assert_eq!(g.phase(i), Rational64::new(-1, 2).into());
        }

        assert_eq!(
            *g.scalar(),
            Scalar::sqrt2_pow((4 - 1) * (4 - 2) / 2) * Scalar::from_phase(Rational64::new(1, 4))
        );

        let h = g.clone();
        let fail = local_comp(&mut g, 1);
        assert!(!fail, "Local comp should not match");
        assert_eq!(g, h);
    }

    #[test]
    fn pivot_1() {
        let mut g = Graph::new();

        for _ in 0..7 {
            g.add_vertex(VType::Z);
        }
        g.set_phase(3, Rational64::new(1, 1));
        for i in 0..3 {
            g.add_edge_with_type(i, 3, EType::H);
        }
        g.add_edge_with_type(3, 4, EType::H);
        for i in 5..7 {
            g.add_edge_with_type(4, i, EType::H);
        }

        // g.set_inputs(vec![0,1,2]);
        // g.set_outputs(vec![5,6]);

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 6);

        let mut h = g.clone();
        let success = pivot(&mut h, 3, 4);
        assert!(success, "Pivot should match");

        assert_eq!(g.to_tensor4(), h.to_tensor4());

        assert_eq!(h.num_vertices(), 5);
        assert_eq!(h.num_edges(), 6);

        assert_eq!(h.phase(0), Rational64::new(0, 1).into());
        assert_eq!(h.phase(6), Rational64::new(1, 1).into());

        let mut inputs: Vec<usize> = Vec::new();
        let mut outputs: Vec<usize> = Vec::new();
        for i in 0..3 {
            let inp = g.add_vertex(VType::B);
            inputs.push(inp);
            g.add_edge(i, inp);
        }

        for i in 5..7 {
            let outp = g.add_vertex(VType::B);
            outputs.push(outp);
            g.add_edge(i, outp);
        }

        g.set_inputs(inputs);
        g.set_outputs(outputs);

        let mut h = g.clone();
        let success = pivot(&mut h, 3, 4);
        assert!(success, "Second pivot should match");
        assert_eq!(g.to_tensor4(), h.to_tensor4());
    }

    #[test]
    fn pivot_2() {
        let mut g = Graph::new();

        for _ in 0..7 {
            g.add_vertex(VType::Z);
        }
        g.set_phase(3, Rational64::new(1, 1));
        g.set_phase(4, Rational64::new(1, 1));
        for i in 0..3 {
            g.add_edge_with_type(i, 3, EType::H);
        }
        g.add_edge_with_type(3, 4, EType::H);
        for i in 5..7 {
            g.add_edge_with_type(4, i, EType::H);
        }

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 6);

        let success = pivot(&mut g, 3, 4);
        assert!(success, "Pivot should match");

        assert_eq!(g.num_vertices(), 5);
        assert_eq!(g.num_edges(), 6);

        assert_eq!(g.phase(0), Rational64::new(1, 1).into());
        assert_eq!(g.phase(6), Rational64::new(1, 1).into());
    }

    #[test]
    fn get_pivot_1() {
        let mut g = Graph::new();

        for _ in 0..7 {
            g.add_vertex(VType::Z);
        }
        g.set_vertex_type(0, VType::B);
        // g.set_phase(3, Rational64::new(1,1));
        // g.set_phase(4, Rational64::new(1,1));
        for i in 0..3 {
            g.add_edge_with_type(i, 3, EType::H);
        }
        g.add_edge_with_type(3, 4, EType::H);
        for i in 5..7 {
            g.add_edge_with_type(4, i, EType::H);
        }
        g.set_inputs(vec![0]);

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 6);

        let mut h = g.clone();
        let success = pivot(&mut h, 3, 4);
        assert!(!success, "Pivot should not match");

        let mut h = g.clone();
        let success = gen_pivot(&mut h, 3, 4);
        assert!(success, "gen_pivot should match");

        println!("g=\n{}\n\nh=\n{}", g.to_dot(), h.to_dot());
        println!("gt=\n{}\n\nht=\n{}", g.to_tensor4(), h.to_tensor4());
        assert_eq!(g.to_tensor4(), h.to_tensor4());
    }

    #[test]
    fn get_pivot_2() {
        let mut g = Graph::new();

        for _ in 0..7 {
            g.add_vertex(VType::Z);
        }
        g.set_phase(3, Rational64::new(1, 1));
        g.set_phase(4, Rational64::new(1, 4));
        for i in 0..3 {
            g.add_edge_with_type(i, 3, EType::H);
        }
        g.add_edge_with_type(3, 4, EType::H);
        for i in 5..7 {
            g.add_edge_with_type(4, i, EType::H);
        }

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 6);

        let mut h = g.clone();
        let success = pivot(&mut h, 3, 4);
        assert!(!success, "Pivot should not match");

        let mut h = g.clone();
        let success = gen_pivot(&mut h, 3, 4);
        assert!(success, "gen_pivot should match");

        println!("g=\n{}\n\nh=\n{}", g.to_dot(), h.to_dot());
        println!("gt=\n{}\n\nht=\n{}", g.to_tensor4(), h.to_tensor4());
        assert_eq!(g.to_tensor4(), h.to_tensor4());
    }

    #[test]
    fn gadget_fusion_1() {
        // fuse gadgets of various sizes
        for n in 1..5 {
            let mut graph = Graph::new();
            let bs: Vec<_> = (0..n).map(|_| graph.add_vertex(VType::B)).collect();
            let vs: Vec<_> = (0..n).map(|_| graph.add_vertex(VType::Z)).collect();
            let gs: Vec<_> = (0..2).map(|_| graph.add_vertex(VType::Z)).collect();
            let ps: Vec<_> = (0..2).map(|_| graph.add_vertex(VType::Z)).collect();
            graph.set_inputs(bs.clone());

            for i in 0..n {
                graph.add_edge(bs[i], vs[i]);
            }
            for (&g, &p) in gs.iter().zip(ps.iter()) {
                graph.add_edge_with_type(g, p, EType::H);
                for &v in &vs {
                    graph.add_edge_with_type(v, g, EType::H);
                }
            }

            graph.set_phase(ps[0], Rational64::new(1, 4));
            graph.set_phase(ps[1], Rational64::new(1, 2));

            let h = graph.clone();

            assert!(gadget_fusion(&mut graph, gs[0], gs[1]));
            assert!(graph
                .find_vertex(|v| graph.phase(v) == Rational64::new(3, 4).into())
                .is_some());
            assert!(graph
                .find_vertex(|v| graph.phase(v) == Rational64::new(1, 4).into())
                .is_none());
            assert!(graph
                .find_vertex(|v| graph.phase(v) == Rational64::new(1, 2).into())
                .is_none());
            // println!("{}", g.to_tensor4());
            // println!("{}", h.to_tensor4());
            // println!("g = {} * \n {} \n\n", g.scalar(), g.to_dot());
            // println!("h = {} * \n {} \n\n", h.scalar(), h.to_dot());
            assert_eq!(graph.to_tensor4(), h.to_tensor4());
        }
    }

    #[test]
    fn scalar_rules() {
        for &t in &[VType::Z, VType::X] {
            let mut g = Graph::new();
            g.add_vertex_with_phase(t, Rational64::new(1, 4));
            let mut h = g.clone();
            assert!(remove_single(&mut h, 0));
            assert_eq!(h.num_vertices(), 0, "h still has vertices");
            assert_eq!(g.to_tensor4(), h.to_tensor4());
        }

        for &t0 in &[VType::Z, VType::X] {
            for &t1 in &[VType::Z, VType::X] {
                for &et in &[EType::N, EType::H] {
                    let mut g = Graph::new();
                    g.add_vertex_with_phase(t0, Rational64::new(1, 4));
                    g.add_vertex_with_phase(t1, Rational64::new(-1, 2));
                    g.add_edge_with_type(0, 1, et);
                    let mut h = g.clone();
                    assert!(remove_pair(&mut h, 0, 1));
                    assert_eq!(h.num_vertices(), 0, "h still has vertices");
                    assert_eq!(
                        g.to_tensor4(),
                        h.to_tensor4(),
                        "Eq failed on case: {:?}, {:?}, {:?}",
                        t0,
                        t1,
                        et
                    );
                }
            }
        }
    }

    #[test]
    fn pi_copy_1() {
        let mut g = Graph::new();
        let vs = [
            g.add_vertex_with_phase(VType::Z, 1),
            g.add_vertex_with_phase(VType::Z, (1, 2)),
            g.add_vertex_with_phase(VType::X, (1, 4)),
            g.add_vertex_with_phase(VType::Z, 0),
            g.add_vertex_with_phase(VType::X, (3, 2)),
            g.add_vertex_with_phase(VType::Z, (3, 4)),
            g.add_vertex_with_phase(VType::B, 0),
        ];

        g.set_inputs(vec![vs[6]]);

        g.add_edge_smart(vs[0], vs[1], EType::H);
        g.add_edge_smart(vs[1], vs[2], EType::N);
        g.add_edge_smart(vs[2], vs[3], EType::H);
        g.add_edge_smart(vs[3], vs[4], EType::N);
        g.add_edge_smart(vs[4], vs[5], EType::N);
        g.add_edge_smart(vs[5], vs[6], EType::N);

        assert!(check_pi_copy(&g, vs[0]));
        assert!(check_pi_copy(&g, vs[1]));
        assert!(!check_pi_copy(&g, vs[2]));
        assert!(!check_pi_copy(&g, vs[3]));
        assert!(check_pi_copy(&g, vs[4]));
        assert!(!check_pi_copy(&g, vs[5]));
        assert!(!check_pi_copy(&g, vs[6]));

        let h = g.clone();
        pi_copy(&mut g, vs[0]);
        assert_eq!(g.to_tensor4(), h.to_tensor4());
        pi_copy(&mut g, vs[1]);
        assert_eq!(g.to_tensor4(), h.to_tensor4());
        pi_copy(&mut g, vs[4]);
        assert_eq!(g.to_tensor4(), h.to_tensor4());
    }
}

// }}}
// vim:foldlevel=0:
