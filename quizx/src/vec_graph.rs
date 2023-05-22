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
use crate::scalar::*;
use num::rational::Rational;
use std::mem;

pub type VTab<T> = Vec<Option<T>>;

#[derive(Debug, Clone, PartialEq)]
pub struct Graph {
    vdata: VTab<VData>,
    edata: VTab<Vec<(V, EType)>>,
    holes: Vec<V>, // places where a vertex has been deleted
    inputs: Vec<V>,
    outputs: Vec<V>,
    numv: usize,
    nume: usize,
    scalar: ScalarN,
}

impl Graph {
    /// Explicitly index neighbors of a vertex. Used for iteration.
    pub fn neighbor_at(&self, v: V, n: usize) -> V {
        let Some(d) = &self.edata[v] else { return 0 };
        d[n].0
    }

    fn index<U>(nhd: &Vec<(V, U)>, v: V) -> Option<usize> {
        nhd.iter().position(|&(v0, _)| v == v0)
    }

    fn value<U: Copy>(nhd: &[(V, U)], v: V) -> Option<U> {
        nhd.iter().find(|t| v == t.0).map(|t| t.1)
    }

    /// Removes vertex 't' from the adjacency map of 's'. This private method
    /// is used by remove_edge and remove_vertex to make the latter slightly
    /// more efficient.
    fn remove_half_edge(&mut self, s: V, t: V) {
        let Some(Some(nhd)) = self.edata.get_mut(s) else { return };
        Self::index(&nhd, t).map(|i| nhd.swap_remove(i));
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
    fn new() -> Self {
        Self {
            vdata: vec![],
            edata: vec![],
            holes: vec![],
            inputs: vec![],
            outputs: vec![],
            numv: 0,
            nume: 0,
            scalar: Scalar::one(),
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

    fn vertices(&self) -> VIter {
        VIter::Vec(self.numv, self.vdata.iter().enumerate())
    }

    fn edges(&self) -> EIter {
        EIter::Vec(self.nume, self.edata.iter().enumerate(), None)
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
            phase: Rational::new(0, 1),
            qubit: 0,
            row: 0,
        })
    }

    fn add_vertex_with_data(&mut self, d: VData) -> V {
        self.numv += 1;
        if let Some(v) = self.holes.pop() {
            self.vdata[v] = Some(d);
            self.edata[v] = Some(vec![]);
            v
        } else {
            self.vdata.push(Some(d));
            self.edata.push(Some(vec![]));
            self.vdata.len() - 1
        }
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

        let Some(Some(nhd)) = self.edata.get_mut(s) else { panic!("Source vertex not found") };
        nhd.push((t, ety));

        let Some(Some(nhd)) = self.edata.get_mut(t) else { panic!("Target vertex not found") };
        nhd.push((s, ety));
    }

    fn remove_edge(&mut self, s: V, t: V) {
        self.nume -= 1;
        self.remove_half_edge(s, t);
        self.remove_half_edge(t, s);
    }

    fn set_phase(&mut self, v: V, phase: Rational) {
        let Some(Some(d)) = self.vdata.get_mut(v) else {  panic!("Vertex not found") };
        d.phase = phase.mod2();
    }

    fn phase(&self, v: V) -> Rational {
        self.vdata[v].expect("Vertex not found").phase
    }

    fn add_to_phase(&mut self, v: V, phase: Rational) {
        let Some(Some(d)) = self.vdata.get_mut(v) else { panic!("Vertex not found") };
        d.phase = (d.phase + phase).mod2();
    }

    fn set_vertex_type(&mut self, v: V, ty: VType) {
        let Some(Some(d)) = self.vdata.get_mut(v) else { panic!("Vertex not found") };
        d.ty = ty;
    }

    fn vertex_data(&self, v: V) -> VData {
        self.vdata[v].expect("Vertex not found")
    }

    fn vertex_type(&self, v: V) -> VType {
        self.vertex_data(v).ty
    }

    fn set_edge_type(&mut self, s: V, t: V, ety: EType) {
        let Some(Some(nhd)) = self.edata.get_mut(s) else { panic!("Source vertex not found") };
        let i = Self::index(&nhd, t).expect("Edge not found");
        nhd[i] = (t, ety);

        let Some(Some(nhd)) = self.edata.get_mut(t) else { panic!("Target vertex not found") };
        let i = Self::index(&nhd, s).expect("Edge not found");
        nhd[i] = (s, ety);
    }

    fn edge_type_opt(&self, s: V, t: V) -> Option<EType> {
        let Some(Some(nhd)) = self.edata.get(s) else { return None };
        Self::value(&nhd, t)
    }

    fn set_coord(&mut self, v: V, coord: (i32, i32)) {
        let Some(Some(d)) = self.vdata.get_mut(v) else { panic!("Vertex not found") };
        d.qubit = coord.0;
        d.row = coord.1;
    }

    fn coord(&mut self, v: V) -> (i32, i32) {
        let d = self.vdata[v].expect("Vertex not found");
        (d.qubit, d.row)
    }

    fn set_qubit(&mut self, v: V, qubit: i32) {
        let Some(Some(d)) = self.vdata.get_mut(v) else { panic!("Vertex not found") };
        d.qubit = qubit;
    }

    fn qubit(&self, v: V) -> i32 {
        self.vdata[v].expect("Vertex not found").qubit
    }

    fn set_row(&mut self, v: V, row: i32) {
        let Some(Some(d)) = self.vdata.get_mut(v) else { panic!("Vertex not found") };
        d.row = row;
    }

    fn row(&self, v: V) -> i32 {
        self.vdata[v].expect("Vertex not found").row
    }

    fn neighbors(&self, v: V) -> NeighborIter {
        let Some(Some(nhd)) = self.edata.get(v) else { panic!("Vertex not found") };
        NeighborIter::Vec(nhd.iter())
    }

    fn incident_edges(&self, v: V) -> IncidentEdgeIter {
        let Some(Some(nhd)) = self.edata.get(v) else { panic!("Vertex not found") };
        IncidentEdgeIter::Vec(nhd.iter())
    }

    fn degree(&self, v: V) -> usize {
        let Some(Some(nhd)) = self.edata.get(v) else { panic!("Vertex not found") };
        nhd.len()
    }

    fn scalar(&self) -> &ScalarN {
        &self.scalar
    }
    fn scalar_mut(&mut self) -> &mut ScalarN {
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
}
