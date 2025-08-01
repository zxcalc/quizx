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

pub use crate::graph::*;
use crate::params::Expr;
use crate::scalar::*;
use num::rational::Rational64;
use rustc_hash::FxHashMap;
use std::mem;

pub type VTab<T> = Vec<Option<T>>;

/// When VecGraph::pack(false) is called, only do the packing if the number of holes * PACK_RATIO
/// is greater than the number of allocated vertex indices.
const PACK_RATIO: usize = 10;

#[derive(Debug, Clone, PartialEq)]
pub struct Graph {
    vdata: VTab<VData>,
    edata: VTab<Vec<(V, EType)>>,
    holes: Vec<V>, // places where a vertex has been deleted
    inputs: Vec<V>,
    outputs: Vec<V>,
    numv: usize,
    nume: usize,
    scalar: Scalar4,
    scalar_factors: FxHashMap<Expr, Scalar4>,
}

impl Default for Graph {
    fn default() -> Self {
        Graph {
            vdata: Vec::new(),
            edata: Vec::new(),
            holes: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            numv: 0,
            nume: 0,
            scalar: Scalar4::one(),
            scalar_factors: FxHashMap::default(),
        }
    }
}

impl Graph {
    /// Explicitly index neighbors of a vertex. Used for iteration.
    pub fn neighbor_at(&self, v: V, n: usize) -> V {
        if let Some(d) = &self.edata[v] {
            d[n].0
        } else {
            0
        }
    }

    fn index<U>(nhd: &[(V, U)], v: V) -> Option<usize> {
        nhd.iter().position(|&(v0, _)| v == v0)
    }

    fn value<U: Copy>(nhd: &[(V, U)], v: V) -> Option<U> {
        for (v0, u) in nhd.iter() {
            if v == *v0 {
                return Some(*u);
            }
        }
        None
    }

    /// Removes vertex 't' from the adjacency map of 's'. This private method
    /// is used by remove_edge and remove_vertex to make the latter slightly
    /// more efficient.
    fn remove_half_edge(&mut self, s: V, t: V) {
        if let Some(Some(nhd)) = self.edata.get_mut(s) {
            Graph::index(nhd, t).map(|i| nhd.swap_remove(i));
        }
    }

    // Here are some simpler implementations of the vertices and edges functions,
    // but they can't be moved into the trait because they return "impl" types.
    // pub fn vertices2(&self) -> impl Iterator<Item=V> + '_ {
    //     self.vdata.iter().enumerate().filter_map(|(v,d)| d.map(|_| v))
    // }

    // pub fn edges2(&self) -> impl Iterator<Item=(V,V,EType)> + '_ {
    //     self.edata
    //         .iter()
    //         .enumerate()
    //         .filter_map(|(v,tab)| match tab {
    //             Some(x) => Some(x.iter().filter_map(move |(v1,et)|
    //                 if v <= *v1 { Some((v, *v1, *et)) } else { None }
    //                 )),
    //             None => None
    //         })
    //         .flatten()
    // }
}

impl GraphLike for Graph {
    fn new() -> Graph {
        Graph {
            vdata: Vec::new(),
            edata: Vec::new(),
            holes: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            numv: 0,
            nume: 0,
            scalar: 1.into(),
            scalar_factors: FxHashMap::default(),
        }
    }

    fn vindex(&self) -> V {
        self.vdata.len()
    }

    fn num_vertices(&self) -> usize {
        self.numv
    }

    fn num_edges(&self) -> usize {
        self.nume
    }

    fn vertices(&self) -> impl Iterator<Item = V> {
        (0..self.vdata.len()).filter(|&i| self.vdata[i].is_some())
    }

    fn edges(&self) -> impl Iterator<Item = (V, V, EType)> {
        self.edata
            .iter()
            .enumerate()
            .filter_map(|(u, tab)| tab.as_ref().map(|tab1| (u, tab1)))
            .flat_map(|(u, tab)| {
                tab.iter()
                    .filter_map(move |&(v, t)| if u <= v { Some((u, v, t)) } else { None })
            })
    }

    fn inputs(&self) -> &Vec<V> {
        &self.inputs
    }
    fn inputs_mut(&mut self) -> &mut Vec<V> {
        &mut self.inputs
    }
    fn set_inputs(&mut self, inputs: Vec<V>) {
        self.inputs = inputs;
    }
    fn outputs(&self) -> &Vec<V> {
        &self.outputs
    }
    fn set_outputs(&mut self, outputs: Vec<V>) {
        self.outputs = outputs;
    }
    fn outputs_mut(&mut self) -> &mut Vec<V> {
        &mut self.outputs
    }

    fn add_vertex(&mut self, ty: VType) -> V {
        self.add_vertex_with_data(VData {
            ty,
            phase: Rational64::new(0, 1).into(),
            ..Default::default()
        })
    }

    fn add_vertex_with_data(&mut self, d: VData) -> V {
        self.numv += 1;
        if let Some(v) = self.holes.pop() {
            self.vdata[v] = Some(d);
            self.edata[v] = Some(Vec::new());
            v
        } else {
            self.vdata.push(Some(d));
            self.edata.push(Some(Vec::new()));
            self.vdata.len() - 1
        }
    }

    fn add_named_vertex_with_data(&mut self, v: V, d: VData) -> Result<(), &str> {
        if v < self.vdata.len() {
            let h = self.holes.iter().position(|&h| h == v);
            if h.is_none() {
                return Err("Vertex already in graph");
            }
            self.holes.remove(h.unwrap());
        } else {
            for i in self.vdata.len()..(v) {
                self.vdata.push(None);
                self.edata.push(None);
                self.holes.push(i);
            }
        }

        self.numv += 1;
        self.vdata[v] = Some(d);
        self.edata[v] = Some(Vec::new());
        Ok(())
    }

    fn remove_vertex(&mut self, v: V) {
        self.numv -= 1;
        self.holes.push(v);

        self.vdata[v] = None;
        let adj = mem::take(&mut self.edata[v]).expect("No such vertex.");

        for (v1, _) in adj {
            self.nume -= 1;
            self.remove_half_edge(v1, v);
        }
    }

    fn add_edge_with_type(&mut self, s: V, t: V, ety: EType) {
        self.nume += 1;
        // if self.connected(s,t) { panic!("introducing parallel edge!"); }

        if let Some(Some(nhd)) = self.edata.get_mut(s) {
            nhd.push((t, ety));
        } else {
            panic!("Source vertex not found");
        }

        if let Some(Some(nhd)) = self.edata.get_mut(t) {
            nhd.push((s, ety));
        } else {
            panic!("Target vertex not found");
        }
    }

    fn remove_edge(&mut self, s: V, t: V) {
        self.nume -= 1;
        self.remove_half_edge(s, t);
        self.remove_half_edge(t, s);
    }

    fn vertex_data(&self, v: V) -> &VData {
        self.vdata[v].as_ref().expect("Vertex not found")
    }

    fn vertex_data_opt(&self, v: V) -> Option<&VData> {
        self.vdata.get(v)?.as_ref()
    }

    fn vertex_data_mut(&mut self, v: V) -> &mut VData {
        let vd = self.vdata.get_mut(v).expect("Vertex not found");
        vd.as_mut().expect("Vertex not found")
    }

    fn set_edge_type(&mut self, s: V, t: V, ety: EType) {
        if let Some(Some(nhd)) = self.edata.get_mut(s) {
            let i = Graph::index(nhd, t).expect("Edge not found");
            nhd[i] = (t, ety);
        } else {
            panic!("Source vertex not found");
        }

        if let Some(Some(nhd)) = self.edata.get_mut(t) {
            let i = Graph::index(nhd, s).expect("Edge not found");
            nhd[i] = (s, ety);
        } else {
            panic!("Target vertex not found");
        }
    }

    fn edge_type_opt(&self, s: V, t: V) -> Option<EType> {
        if let Some(Some(nhd)) = self.edata.get(s) {
            Graph::value(nhd, t)
        } else {
            None
        }
    }

    fn set_coord(&mut self, v: V, coord: impl Into<Coord>) {
        let coord = coord.into();
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.qubit = coord.y;
            d.row = coord.x;
        } else {
            panic!("Vertex not found")
        }
    }

    fn coord(&self, v: V) -> Coord {
        let d = self.vdata[v].as_ref().expect("Vertex not found");
        Coord::new(d.row, d.qubit)
    }

    fn set_qubit(&mut self, v: V, qubit: f64) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.qubit = qubit;
        } else {
            panic!("Vertex not found")
        }
    }

    fn qubit(&self, v: V) -> f64 {
        self.vdata[v].as_ref().expect("Vertex not found").qubit
    }

    fn set_row(&mut self, v: V, row: f64) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.row = row;
        } else {
            panic!("Vertex not found")
        }
    }

    fn row(&self, v: V) -> f64 {
        self.vdata[v].as_ref().expect("Vertex not found").row
    }

    fn neighbors(&self, v: V) -> impl Iterator<Item = usize> {
        if let Some(Some(nhd)) = self.edata.get(v) {
            nhd.iter().map(|(u, _)| *u)
        } else {
            panic!("Vertex not found")
        }
    }

    fn incident_edges(&self, v: V) -> impl Iterator<Item = (usize, EType)> {
        if let Some(Some(nhd)) = self.edata.get(v) {
            nhd.iter().cloned()
        } else {
            panic!("Vertex not found")
        }
    }

    fn degree(&self, v: V) -> usize {
        if let Some(Some(nhd)) = self.edata.get(v) {
            nhd.len()
        } else {
            panic!("Vertex not found")
        }
    }

    fn scalar(&self) -> &Scalar4 {
        &self.scalar
    }
    fn scalar_mut(&mut self) -> &mut Scalar4 {
        &mut self.scalar
    }

    fn find_edge<F>(&self, f: F) -> Option<(V, V, EType)>
    where
        F: Fn(V, V, EType) -> bool,
    {
        for (v0, nhd0) in self.edata.iter().enumerate() {
            if let Some(nhd) = nhd0 {
                for &(v1, et) in nhd.iter() {
                    if v0 <= v1 && f(v0, v1, et) {
                        return Some((v0, v1, et));
                    }
                }
            };
        }

        None
    }

    fn find_vertex<F>(&self, f: F) -> Option<V>
    where
        F: Fn(V) -> bool,
    {
        for (v, d) in self.vdata.iter().enumerate() {
            if d.is_some() && f(v) {
                return Some(v);
            }
        }

        None
    }

    fn contains_vertex(&self, v: V) -> bool {
        v < self.vdata.len() && self.vdata[v].is_some()
    }

    fn scalar_factors(&self) -> impl Iterator<Item = (&Expr, &Scalar4)> {
        self.scalar_factors.iter()
    }

    fn get_scalar_factor(&self, e: &Expr) -> Option<Scalar4> {
        self.scalar_factors.get(e).copied()
    }

    fn mul_scalar_factor(&mut self, e: Expr, s: Scalar4) {
        if let Some(t) = self.scalar_factors.get_mut(&e) {
            *t *= s;
        } else {
            self.scalar_factors.insert(e, s);
        }
    }

    #[allow(clippy::needless_range_loop)]
    fn pack(&mut self, force: bool) {
        if force || self.holes.len() * PACK_RATIO > self.vdata.len() {
            let new_size = self.num_vertices();
            let mut vtab = vec![0; self.vdata.len()];
            let mut j = 0;
            for i in 0..self.vdata.len() {
                if self.vdata[i].is_some() {
                    self.vdata[j] = self.vdata[i].take();
                    self.edata[j] = self.edata[i].take();
                    vtab[i] = j;
                    j += 1;
                }
            }

            self.vdata.truncate(new_size);
            self.edata.truncate(new_size);
            self.holes = vec![];

            for et in self.edata.iter_mut() {
                if let Some(tab) = et.as_mut() {
                    tab.iter_mut()
                        .for_each(|pair| *pair = (vtab[pair.0], pair.1));
                }
            }

            self.inputs = self.inputs.iter().map(|v| vtab[*v]).collect();
            self.outputs = self.outputs.iter().map(|v| vtab[*v]).collect();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_empty_graph() {
        let g = Graph::new();
        assert_eq!(g.num_vertices(), 0);
        assert_eq!(g.num_edges(), 0);
    }

    fn simple_graph() -> (Graph, Vec<V>) {
        let mut g = Graph::new();
        let vs = vec![
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::X),
            g.add_vertex(VType::X),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];
        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[3]);
        g.add_edge(vs[2], vs[4]);
        g.add_edge(vs[2], vs[5]);
        g.add_edge(vs[3], vs[4]);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[4], vs[6]);
        g.add_edge(vs[5], vs[7]);
        (g, vs)
    }

    #[test]
    fn create_simple_graph() {
        let (g, _) = simple_graph();
        assert_eq!(g.num_vertices(), 8);
        assert_eq!(g.num_edges(), 8);
    }

    #[test]
    fn clone_graph() {
        let (g, _) = simple_graph();
        let h = g.clone();
        assert!(g.num_vertices() == h.num_vertices());
        assert!(g.num_edges() == h.num_edges());
        // assert!(g == h);
    }

    #[test]
    fn vertex_iterator() {
        let (g, mut expected_vs) = simple_graph();
        let mut vs: Vec<_> = g.vertices().collect();
        vs.sort();
        expected_vs.sort();
        assert_eq!(expected_vs, vs);
    }

    #[test]
    fn edge_iterator() {
        let (mut g, vs) = simple_graph();
        g.set_edge_type(vs[1], vs[3], EType::H);

        let mut edges: Vec<_> = g.edges().collect();
        let mut expected_edges = vec![
            (vs[0], vs[2], EType::N),
            (vs[1], vs[3], EType::H),
            (vs[2], vs[4], EType::N),
            (vs[2], vs[5], EType::N),
            (vs[3], vs[4], EType::N),
            (vs[3], vs[5], EType::N),
            (vs[4], vs[6], EType::N),
            (vs[5], vs[7], EType::N),
        ];

        edges.sort();
        expected_edges.sort();
        assert_eq!(expected_edges, edges);
    }

    #[test]
    fn smart_edges_zx() {
        let mut g = Graph::new();
        let vs = [
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::X),
            g.add_vertex(VType::B),
        ];
        g.add_edge(vs[0], vs[1]);
        g.add_edge(vs[2], vs[3]);

        let mut h = g.clone();
        h.add_edge_smart(vs[1], vs[2], EType::N);
        h.add_edge_smart(vs[1], vs[2], EType::N);
        assert_eq!(h.num_vertices(), 4);
        // assert_eq!(h.num_edges(), 2,
        //     "Wrong edges in NN test: {:?}",
        //     Vec::from_iter(h.edges()));

        let mut h = g.clone();
        h.add_edge_smart(vs[1], vs[2], EType::H);
        h.add_edge_smart(vs[1], vs[2], EType::H);
        assert_eq!(h.num_vertices(), 4);
        // assert_eq!(h.num_edges(), 3,
        //     "Wrong edges in HH test: {:?}",
        //     Vec::from_iter(h.edges()));
        assert_eq!(h.edge_type(vs[1], vs[2]), EType::H);
    }

    #[test]
    fn pack() {
        let mut g = Graph::new();
        g.add_vertex(VType::X);
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::X);
        g.add_vertex(VType::B);
        g.add_edge(0, 1);
        g.add_edge(0, 2);
        g.add_edge(0, 4);

        g.remove_vertex(1);
        g.remove_vertex(3);
        assert!(g.contains_vertex(0));
        assert!(!g.contains_vertex(1));
        assert!(g.contains_vertex(2));
        assert!(!g.contains_vertex(3));
        assert!(g.contains_vertex(4));
        assert!(g.connected(0, 2));
        assert!(g.connected(0, 4));
        assert!(g.vertex_type(0) == VType::X);
        assert!(g.vertex_type(2) == VType::Z);
        assert!(g.vertex_type(4) == VType::B);

        g.pack(true);
        assert!(g.contains_vertex(0));
        assert!(g.contains_vertex(1));
        assert!(g.contains_vertex(2));
        assert!(!g.contains_vertex(3));
        assert!(!g.contains_vertex(4));
        assert!(g.connected(0, 1));
        assert!(g.connected(0, 2));
        assert!(g.vertex_type(0) == VType::X);
        assert!(g.vertex_type(1) == VType::Z);
        assert!(g.vertex_type(2) == VType::B);
    }
}
