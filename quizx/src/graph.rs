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

use crate::params::Expr;
use crate::phase::Phase;
use crate::util::*;
use crate::{fscalar::*, params::Parity};
use derive_more::{Display, From};
use num::rational::Rational64;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::iter::FromIterator;

pub type V = usize;

/// The type of a vertex in a graph.
///
/// The serialized names may differ.
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum VType {
    B, // Boundary
    #[default]
    Z, // Z-spider
    X, // X-spider
    #[serde(rename = "hadamard")]
    H, // H-box
    #[serde(rename = "W_input")]
    WInput,
    #[serde(rename = "W_output")]
    WOutput,
    #[serde(rename = "Z_box")]
    ZBox,
}

#[derive(Debug, Clone, PartialEq)]
pub struct VData {
    pub ty: VType,
    pub phase: Phase,
    pub vars: Parity,
    pub qubit: f64,
    pub row: f64,
}

impl Default for VData {
    fn default() -> Self {
        VData {
            ty: VType::B,
            phase: Phase::zero(),
            vars: Parity::zero(),
            qubit: 0.0,
            row: 0.0,
        }
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum EType {
    /// Normal edge.
    #[default]
    #[serde(rename = "simple")]
    N,
    /// Hadamard edge.
    #[serde(rename = "hadamard")]
    H,
    /// W input/output
    #[serde(rename = "w_io")]
    Wio,
}

impl EType {
    pub fn opposite(&self) -> EType {
        match self {
            EType::N => EType::H,
            EType::H => EType::N,
            EType::Wio => EType::Wio,
        }
    }

    pub fn merge(et0: EType, et1: EType) -> EType {
        if et0 == EType::N {
            et1
        } else {
            et1.opposite()
        }
    }
}

/// An enum specifying an X or Z basis element
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum BasisElem {
    Z0, // |0>
    Z1, // |1>
    X0, // |+>
    X1, // |->
    SKIP,
}

impl BasisElem {
    pub fn phase(&self) -> Rational64 {
        if *self == BasisElem::Z1 || *self == BasisElem::X1 {
            Rational64::one()
        } else {
            Rational64::zero()
        }
    }

    pub fn is_z(&self) -> bool {
        *self == BasisElem::Z0 || *self == BasisElem::Z1
    }

    pub fn is_x(&self) -> bool {
        *self == BasisElem::X0 || *self == BasisElem::X1
    }

    pub fn flipped(&self) -> BasisElem {
        match self {
            BasisElem::Z0 => BasisElem::Z1,
            BasisElem::Z1 => BasisElem::Z0,
            BasisElem::X0 => BasisElem::X1,
            BasisElem::X1 => BasisElem::X0,
            BasisElem::SKIP => BasisElem::SKIP,
        }
    }
}

/// Coordinates for rendering a node.
#[derive(Display, Debug, Default, Clone, Copy, PartialEq, PartialOrd, From)]
#[display("({},{})", x, y)]
pub struct Coord {
    pub x: f64,
    pub y: f64,
}

impl Coord {
    /// Create a new coordinate.
    pub fn new(x: f64, y: f64) -> Self {
        Coord { x, y }
    }

    /// Infer the qubit index from the y-coordinate.
    pub fn qubit(&self) -> f64 {
        self.y
    }

    /// Infer the row index from the x-coordinate.
    pub fn row(&self) -> f64 {
        self.x
    }
}

pub trait GraphLike: Clone + Sized + Send + Sync + std::fmt::Debug {
    /// Initialise a new empty graph
    fn new() -> Self;

    /// Next fresh vertex index
    fn vindex(&self) -> V;

    /// Number of vertices
    fn num_vertices(&self) -> usize;

    /// Number of edges
    fn num_edges(&self) -> usize;

    /// Get iterator over all vertices
    fn vertices(&self) -> impl Iterator<Item = V>;

    /// Get iterator over all edges
    ///
    /// An "edge" is a triple (s, t, edge_type), where s <= t.
    fn edges(&self) -> impl Iterator<Item = (V, V, EType)>;

    /// List of boundary vertices which serve as inputs
    fn inputs(&self) -> &Vec<V>;

    /// Mutable list of boundary vertices which serve as inputs
    fn inputs_mut(&mut self) -> &mut Vec<V>;

    /// Set inputs for the graph
    fn set_inputs(&mut self, inputs: Vec<V>);

    /// List of boundary vertices which serve as outputs
    fn outputs(&self) -> &Vec<V>;

    /// Mutable list of boundary vertices which serve as outputs
    fn outputs_mut(&mut self) -> &mut Vec<V>;

    /// Set outputs for the graph
    fn set_outputs(&mut self, outputs: Vec<V>);

    /// Add a vertex with the given type
    fn add_vertex(&mut self, ty: VType) -> V;

    /// Add a vertex with the given VData struct
    fn add_vertex_with_data(&mut self, d: VData) -> V;

    /// Add a vertex with the given name and VData struct
    ///
    /// Returns an error if vertex already exists
    fn add_named_vertex_with_data(&mut self, v: V, d: VData) -> Result<(), &str>;

    /// Remove a vertex from a graph
    ///
    /// Behavior is undefined if the vertex is not in the graph.
    fn remove_vertex(&mut self, v: V);

    /// Add an edge with the given type
    ///
    /// Behaviour is undefined if an edge already exists between s and t.
    fn add_edge_with_type(&mut self, s: V, t: V, ety: EType);

    /// Remove an edge from a graph
    ///
    /// Behaviour is undefined if there is no edge between s and t.
    fn remove_edge(&mut self, s: V, t: V);

    /// Get the data associated to the given vertex
    fn vertex_data(&self, v: V) -> &VData;

    /// Get a mutable ref to the data associated to the given vertex
    fn vertex_data_mut(&mut self, v: V) -> &mut VData;
    fn set_edge_type(&mut self, s: V, t: V, ety: EType);
    fn edge_type_opt(&self, s: V, t: V) -> Option<EType>;
    fn neighbors(&self, v: V) -> impl Iterator<Item = V>;
    fn incident_edges(&self, v: V) -> impl Iterator<Item = (V, EType)>;
    fn degree(&self, v: V) -> usize;
    fn scalar(&self) -> &FScalar;
    fn scalar_mut(&mut self) -> &mut FScalar;
    fn find_edge<F>(&self, f: F) -> Option<(V, V, EType)>
    where
        F: Fn(V, V, EType) -> bool;
    fn find_vertex<F>(&self, f: F) -> Option<V>
    where
        F: Fn(V) -> bool;
    fn contains_vertex(&self, v: V) -> bool;

    fn scalar_factors(&self) -> impl Iterator<Item = (&Expr, &FScalar)>;
    fn get_scalar_factor(&self, e: &Expr) -> Option<FScalar>;
    fn mul_scalar_factor(&mut self, e: Expr, s: FScalar);

    /// Returns the phase and any boolean variables at a vertex
    fn phase_and_vars(&self, v: V) -> (Phase, Parity) {
        let vd = self.vertex_data(v);
        (vd.phase, vd.vars.clone())
    }

    /// Set the phase of a vertex
    fn set_phase(&mut self, v: V, phase: impl Into<Phase>) {
        self.vertex_data_mut(v).phase = phase.into();
    }

    /// Returns the phase of vertex `v`
    fn phase(&self, v: V) -> Phase {
        self.vertex_data(v).phase
    }

    /// Adds a value to the phase of a vertex
    fn add_to_phase(&mut self, v: V, phase: impl Into<Phase>) {
        let vd = self.vertex_data_mut(v);
        vd.phase = (vd.phase + phase.into()).normalize();
    }

    /// Sets the type of a vertex
    fn set_vertex_type(&mut self, v: V, ty: VType) {
        self.vertex_data_mut(v).ty = ty;
    }

    /// Returns the type of a vertex
    fn vertex_type(&self, v: V) -> VType {
        self.vertex_data(v).ty
    }

    /// Sets the coordinate of a vertex
    ///
    /// This method takes a Coord as an argument and sets the qubit of the vertex to `coord.y` and
    /// the row of the vertex to `coord.x`.
    fn set_coord(&mut self, v: V, coord: impl Into<Coord>) {
        let coord = coord.into();
        let d = self.vertex_data_mut(v);
        d.qubit = coord.y;
        d.row = coord.x;
    }

    /// Returns the coordinate of a vertex
    ///
    /// This method returns a Coord for a given vertex, where the qubit of the vertex is `coord.y` and
    /// the row of the vertex is `coord.x`.
    fn coord(&self, v: V) -> Coord {
        let d = self.vertex_data(v);
        Coord::new(d.row, d.qubit)
    }

    /// Sets the qubit index of the given vertex
    ///
    /// This is primarily used for visual layout of a vertex, hence `qubit` and `row` are allowed to take
    /// fractional values to allow arbitrary placements in 2D space.
    fn set_qubit(&mut self, v: V, qubit: f64) {
        self.vertex_data_mut(v).qubit = qubit;
    }

    /// Returns the qubit index of the given vertex
    fn qubit(&self, v: V) -> f64 {
        self.vertex_data(v).qubit
    }

    /// Sets the row of the given vertex
    fn set_row(&mut self, v: V, row: f64) {
        self.vertex_data_mut(v).row = row;
    }

    /// Returns the row of the given vertex
    fn row(&self, v: V) -> f64 {
        self.vertex_data(v).row
    }

    /// Sets the boolean variables that affect the phase of this vertex
    ///
    /// This allows the phase to depend on an XOR of boolean variables, represented by a
    /// `Parity`. If the parity is even, the phase is considered to be `self.vertex_data(v).phase`,
    /// whereas if the parity is odd, the phase should be `self.vertex_data(v).phase + Phase::one()`,
    /// i.e. it gains a pi term.
    fn set_vars(&mut self, v: V, vars: Parity) {
        self.vertex_data_mut(v).vars = vars;
    }

    /// Returns the boolean variables that affect the phase of this vertex
    fn vars(&self, v: V) -> Parity {
        self.vertex_data(v).vars.clone()
    }

    /// Adds the given variables to the parity expression of the vertex
    fn add_to_vars(&mut self, v: V, vars: &Parity) {
        let vars1 = &self.vertex_data(v).vars + vars;
        self.vertex_data_mut(v).vars = vars1;
    }

    fn add_edge(&mut self, s: V, t: V) {
        self.add_edge_with_type(s, t, EType::N);
    }

    fn edge_type(&self, s: V, t: V) -> EType {
        self.edge_type_opt(s, t).expect("Edge not found")
    }

    fn connected(&self, v0: V, v1: V) -> bool {
        self.edge_type_opt(v0, v1).is_some()
    }

    fn toggle_edge_type(&mut self, v0: V, v1: V) {
        self.set_edge_type(v0, v1, self.edge_type(v0, v1).opposite());
    }

    fn vertex_vec(&self) -> Vec<V> {
        self.vertices().collect()
    }
    fn edge_vec(&self) -> Vec<(V, V, EType)> {
        self.edges().collect()
    }
    fn neighbor_vec(&self, v: V) -> Vec<V> {
        self.neighbors(v).collect()
    }
    fn incident_edge_vec(&self, v: V) -> Vec<(V, EType)> {
        self.incident_edges(v).collect()
    }

    /// Convert all X spiders to Z with the colour-change rule
    fn x_to_z(&mut self) {
        for v in Vec::from_iter(self.vertices()) {
            if self.vertex_type(v) == VType::X {
                self.set_vertex_type(v, VType::Z);
                for w in Vec::from_iter(self.neighbors(v)) {
                    self.toggle_edge_type(v, w);
                }
            }
        }
    }

    /// Add a vertex to the graph with the given type and phase
    fn add_vertex_with_phase(&mut self, ty: VType, phase: impl Into<Phase>) -> V {
        let v = self.add_vertex(ty);
        self.set_phase(v, phase.into());
        v
    }

    /// Add an edge and simplify if necessary to remove parallel edges
    ///
    /// The behaviour of this function depends on the type of source/target
    /// vertex as well as the type of the existing edge (if there is one).
    fn add_edge_smart(&mut self, s: V, t: V, ety: EType) {
        let st = self.vertex_type(s);
        if s == t {
            if st == VType::Z || st == VType::X {
                if ety == EType::H {
                    self.add_to_phase(s, Phase::one());
                    self.scalar_mut().mul_sqrt2_pow(-1);
                }
            } else {
                panic!("Self-loops only supported on Z and X nodes");
            }
        } else if let Some(ety0) = self.edge_type_opt(s, t) {
            let tt = self.vertex_type(t);
            match (st, tt) {
                (VType::Z, VType::Z) | (VType::X, VType::X) => {
                    match (ety0, ety) {
                        (EType::N, EType::N) => {} // ignore new edge
                        (EType::H, EType::H) => {
                            self.remove_edge(s, t);
                            self.scalar_mut().mul_sqrt2_pow(-2);
                        }
                        (EType::H, EType::N) => {
                            self.set_edge_type(s, t, EType::N);
                            self.add_to_phase(s, Rational64::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::N, EType::H) => {
                            self.add_to_phase(s, Rational64::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::Wio, _) | (_, EType::Wio) => {
                            unimplemented!("W nodes not supported")
                        }
                    }
                }
                (VType::Z, VType::X) | (VType::X, VType::Z) => {
                    match (ety0, ety) {
                        (EType::N, EType::N) => {
                            self.remove_edge(s, t);
                            self.scalar_mut().mul_sqrt2_pow(-2);
                        }
                        (EType::N, EType::H) => {
                            self.set_edge_type(s, t, EType::H);
                            self.add_to_phase(s, Rational64::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::H, EType::N) => {
                            self.add_to_phase(s, Rational64::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::H, EType::H) => {} // ignore new edge
                        (EType::Wio, _) | (_, EType::Wio) => {
                            unimplemented!("W nodes not supported")
                        }
                    }
                }
                _ => panic!(
                    "Parallel edges only supported between Z and X vertices ({st:?} --> {tt:?})"
                ),
            }
        } else {
            self.add_edge_with_type(s, t, ety);
        }
    }

    /// Replace a boundary vertex with the given basis element
    ///
    /// Note this does not replace the vertex from the input/output list or do
    /// normalisation.
    fn plug_vertex(&mut self, v: V, b: BasisElem) {
        if b != BasisElem::SKIP {
            self.set_vertex_type(v, VType::Z);
            self.set_phase(v, b.phase());

            if b.is_z() {
                let n = self
                    .neighbors(v)
                    .next()
                    .expect("Boundary should have 1 neighbor.");
                self.toggle_edge_type(v, n);
            }
        }
    }

    /// Plug the given basis vertex into the i-th output.
    fn plug_output(&mut self, i: usize, b: BasisElem) {
        self.plug_vertex(self.outputs()[i], b);
        self.outputs_mut().remove(i);
        self.scalar_mut().mul_sqrt2_pow(-1);
    }

    /// Plug the given basis vertex into the i-th input.
    fn plug_input(&mut self, i: usize, b: BasisElem) {
        self.plug_vertex(self.inputs()[i], b);
        self.inputs_mut().remove(i);
        self.scalar_mut().mul_sqrt2_pow(-1);
    }

    /// Plug the given list of normalised basis elements in as inputs, starting from the left
    ///
    /// The list `plug` should be of length <= the number of inputs.
    fn plug_inputs(&mut self, plug: &[BasisElem]) {
        let mut inp: Vec<V> = vec![];
        let mut num_plugged = 0;
        let inputs = self.inputs().clone();
        for (i, &v) in inputs.iter().enumerate() {
            if plug[i] != BasisElem::SKIP && i < plug.len() {
                self.plug_vertex(v, plug[i]);
                num_plugged += 1;
            } else {
                inp.push(v);
            }
        }
        self.set_inputs(inp);
        self.scalar_mut().mul_sqrt2_pow(-num_plugged);
    }

    /// Plug the given list of normalised basis elements in as outputs, starting from the left
    ///
    /// The list `plug` should of length <= the number of outputs.
    fn plug_outputs(&mut self, plug: &[BasisElem]) {
        let mut outp: Vec<V> = vec![];
        let mut num_plugged = 0;
        let outputs = self.outputs().clone();
        for (i, &v) in outputs.iter().enumerate() {
            if plug[i] != BasisElem::SKIP && i < plug.len() {
                self.plug_vertex(v, plug[i]);
                num_plugged += 1;
            } else {
                outp.push(v);
            }
        }
        self.set_outputs(outp);
        self.scalar_mut().mul_sqrt2_pow(-num_plugged);
    }

    /// Appends the given graph to the current one, with fresh names.
    ///
    /// The renaming map is returned. The scalars are multiplied, but the inputs/outputs
    /// of `self` are NOT updated.
    fn append_graph(&mut self, other: &impl GraphLike) -> FxHashMap<V, V> {
        let mut vmap = FxHashMap::default();

        for v in other.vertices() {
            let v1 = self.add_vertex_with_data(other.vertex_data(v).clone());
            vmap.insert(v, v1);
        }

        for (v0, v1, et) in other.edges() {
            self.add_edge_with_type(vmap[&v0], vmap[&v1], et);
        }

        *self.scalar_mut() *= other.scalar();

        vmap
    }

    /// Plug the given graph into the outputs and multiply scalars
    ///
    /// Panics if the outputs of `self` are not the same length as the inputs of `other`.
    fn plug(&mut self, other: &impl GraphLike) {
        if other.inputs().len() != self.outputs().len() {
            panic!("Outputs and inputs must match");
        }

        let vmap = self.append_graph(other);

        for k in 0..self.outputs().len() {
            let o = self.outputs()[k];
            let i = other.inputs()[k];
            let (no, et0) = self
                .incident_edges(o)
                .next()
                .unwrap_or_else(|| panic!("Bad output: {o}"));
            let (ni, et1) = other
                .incident_edges(i)
                .next()
                .unwrap_or_else(|| panic!("Bad input: {i}"));
            let et = EType::merge(et0, et1);

            self.add_edge_smart(no, vmap[&ni], et);
            self.remove_vertex(o);
            self.remove_vertex(vmap[&i]);
        }

        let outp = other.outputs().iter().map(|o| vmap[o]).collect();
        self.set_outputs(outp);
    }

    /// Checks if the given graph only consists of wires from the inputs to outputs (in order)
    fn is_identity(&self) -> bool {
        let n = self.inputs().len();
        self.inputs().len() == n
            && self.outputs().len() == n
            && self.num_vertices() == 2 * n
            && (0..n).all(|i| self.connected(self.inputs()[i], self.outputs()[i]))
    }

    /// Return number of Z or X spiders with non-Clifford phase
    fn tcount(&self) -> usize {
        let mut n = 0;
        for v in self.vertices() {
            let t = self.vertex_type(v);
            if (t == VType::Z || t == VType::X) && !self.phase(v).is_clifford() {
                n += 1;
            }
        }
        n
    }

    /// Return a graphviz-friendly string representation of the graph
    fn to_dot(&self) -> String {
        let mut dot = String::from("graph {\n");
        for v in self.vertices() {
            let t = self.vertex_type(v);
            let p = self.phase(v);
            dot += &format!(
                "  {} [color={}, label=\"{}\"",
                v,
                match t {
                    VType::B => "black",
                    VType::Z => "green",
                    VType::X => "red",
                    VType::H => "yellow",
                    VType::WInput => "blue",
                    VType::WOutput => "blue",
                    VType::ZBox => "purple",
                },
                if self.inputs().contains(&v) {
                    format!("{v}:i")
                } else if self.outputs().contains(&v) {
                    format!("{v}:o")
                } else if !p.is_zero() {
                    format!("{v}:{p}")
                } else {
                    format!("{v}")
                }
            );
            let q = self.qubit(v);
            let r = self.row(v);
            if q != 0.0 || r != 0.0 {
                dot += &format!(", pos=\"{q},{r}!\"");
            }
            dot += "]\n";
        }

        dot += "\n";

        for (s, t, ty) in self.edges() {
            dot += &format!("  {s} -- {t}");
            if ty == EType::H {
                dot += " [color=blue]";
            }
            dot += "\n";
        }

        dot += "}\n";

        dot
    }

    /// Exchange inputs and outputs and reverse all phases
    fn adjoint(&mut self) {
        for v in self.vertex_vec() {
            let p = self.phase(v);
            self.set_phase(v, -p);
        }

        let inp = self.inputs().clone();
        self.set_inputs(self.outputs().clone());
        self.set_outputs(inp);
        let s = self.scalar().conj();
        *(self.scalar_mut()) = s;
    }

    /// Same as GraphLike::adjoint(), but return as a copy
    fn to_adjoint(&self) -> Self {
        let mut g = self.clone();
        g.adjoint();
        g
    }

    /// Returns vertices in the components of g
    fn component_vertices(&self) -> Vec<FxHashSet<V>> {
        // vec of vecs storing components
        let mut comps = vec![];

        // set of vertices left to visit
        let mut vset: FxHashSet<V> = self.vertices().collect();

        // stack used in the DFS
        let mut stack = vec![];

        while let Some(&v) = vset.iter().next() {
            // start a new component
            comps.push(FxHashSet::default());
            let i = comps.len() - 1;

            // fill last vec in comps by DFS
            stack.push(v);
            while let Some(v) = stack.pop() {
                comps[i].insert(v);
                vset.remove(&v);

                for w in self.neighbors(v) {
                    if vset.contains(&w) {
                        stack.push(w);
                    }
                }
            }
        }

        comps
    }

    /// Returns the full subgraph containing the given vertices
    fn subgraph_from_vertices(&self, verts: Vec<V>) -> Self {
        let mut g = Self::new();
        let mut vert_map: FxHashMap<V, V> = FxHashMap::default();
        for v in verts {
            let w = g.add_vertex_with_data(self.vertex_data(v).clone());
            vert_map.insert(v, w);
        }

        for (s, t, ety) in self.edges() {
            if vert_map.contains_key(&s) && vert_map.contains_key(&t) {
                g.add_edge_with_type(vert_map[&s], vert_map[&t], ety);
            }
        }

        g
    }

    /// Returns max row of any vertex
    fn depth(&self) -> f64 {
        pmax(self.vertices().map(|v| self.row(v))).unwrap_or(-1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::ToTensor;
    use crate::vec_graph::Graph;
    #[test]
    fn smart_edges() {
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::X);
        g.add_vertex(VType::B);
        g.add_edge_smart(0, 1, EType::N);
        g.add_edge_smart(1, 2, EType::N);
        g.add_edge_smart(2, 3, EType::N);
        g.add_edge_smart(1, 3, EType::N);
        g.add_edge_smart(3, 4, EType::N);
        g.set_inputs(vec![0]);
        g.set_outputs(vec![4]);

        let mut h = Graph::new();
        h.add_vertex(VType::B);
        h.add_vertex(VType::Z);
        h.add_vertex(VType::X);
        h.add_vertex(VType::B);
        h.add_edge_smart(0, 1, EType::N);
        h.add_edge_smart(1, 2, EType::N);
        h.add_edge_smart(1, 2, EType::N);
        h.add_edge_smart(2, 3, EType::N);
        h.set_inputs(vec![0]);
        h.set_outputs(vec![3]);

        let tg = g.to_tensorf();
        let th = h.to_tensorf();
        println!("\n\ntg =\n{tg}");
        println!("\n\nth =\n{th}");
        assert_eq!(tg, th);
    }

    #[test]
    fn plugs() {
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(1, 3);
        g.set_inputs(vec![0]);
        g.set_outputs(vec![2, 3]);

        let mut h = Graph::new();
        h.add_vertex(VType::B);
        h.add_vertex(VType::B);
        h.add_vertex(VType::Z);
        h.add_vertex(VType::B);
        h.add_edge(0, 2);
        h.add_edge(1, 2);
        h.add_edge(2, 3);
        h.set_inputs(vec![0, 1]);
        h.set_outputs(vec![3]);

        g.plug(&h);
        assert_eq!(g.num_vertices(), 4);
        assert_eq!(g.num_edges(), 3);
        let zs: Vec<_> = g
            .vertices()
            .filter(|&v| g.vertex_type(v) == VType::Z)
            .collect();
        assert_eq!(zs.len(), 2);
        assert!(g.connected(zs[0], zs[1]));
    }

    #[test]
    fn dedupe() {
        let mut g: Graph = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::B);

        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(1, 3);
        g.add_edge(0, 3);

        assert_eq!(g.component_vertices().first().unwrap().len(), 4)
    }
}
