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
use crate::json::JsonGraph;
use crate::phase::Phase;
use crate::scalar::*;
use num::rational::Rational64;
use rustc_hash::FxHashMap;
use serde::de::Error as _;
use serde::ser::Error as _;
use std::iter::FromIterator;

pub type VTab<T> = FxHashMap<V, T>;

#[derive(Debug, Clone, PartialEq)]
pub struct Graph {
    vdata: VTab<VData>,
    edata: VTab<VTab<EType>>,
    inputs: Vec<V>,
    outputs: Vec<V>,
    numv: usize,
    nume: usize,
    freshv: V,
    scalar: ScalarN,
}

impl Graph {
    /// Removes vertex 't' from the adjacency map of 's'. This private method
    /// is used by remove_edge and remove_vertex to make the latter slightly
    /// more efficient.
    fn remove_half_edge(&mut self, s: V, t: V) {
        self.edata.get_mut(&s).map(|nhd| nhd.remove(&t));
    }
}

impl GraphLike for Graph {
    fn new() -> Graph {
        Graph {
            vdata: FxHashMap::default(),
            edata: FxHashMap::default(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            numv: 0,
            nume: 0,
            freshv: 0,
            scalar: Scalar::one(),
        }
    }

    fn vindex(&self) -> V {
        self.freshv
    }

    fn num_vertices(&self) -> usize {
        self.numv
    }

    fn num_edges(&self) -> usize {
        self.nume
    }

    fn vertices(&self) -> VIter {
        VIter::Hash(self.vdata.keys())
    }

    fn edges(&self) -> EIter {
        EIter::Hash(self.nume, self.edata.iter(), None)
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
            qubit: 0,
            row: 0,
        })
    }

    fn add_vertex_with_data(&mut self, d: VData) -> V {
        let v = self.freshv;
        self.freshv += 1;
        self.numv += 1;
        self.vdata.insert(v, d);
        self.edata.insert(v, FxHashMap::default());
        v
    }

    fn remove_vertex(&mut self, v: V) {
        self.numv -= 1;

        for v1 in Vec::from_iter(self.neighbors(v)) {
            self.nume -= 1;
            self.remove_half_edge(v1, v);
        }

        self.vdata.remove(&v);
        self.edata.remove(&v);
    }

    fn add_edge_with_type(&mut self, s: V, t: V, ety: EType) {
        self.nume += 1;

        self.edata
            .get_mut(&s)
            .expect("Source vertex not found")
            .insert(t, ety);
        self.edata
            .get_mut(&t)
            .expect("Target vertex not found")
            .insert(s, ety);
    }

    fn remove_edge(&mut self, s: V, t: V) {
        self.nume -= 1;
        self.remove_half_edge(s, t);
        self.remove_half_edge(t, s);
    }

    fn set_phase(&mut self, v: V, phase: impl Into<Phase>) {
        self.vdata.get_mut(&v).expect("Vertex not found").phase = phase.into();
    }

    fn phase(&self, v: V) -> Phase {
        self.vdata.get(&v).expect("Vertex not found").phase
    }

    fn add_to_phase(&mut self, v: V, phase: impl Into<Phase>) {
        if let Some(d) = self.vdata.get_mut(&v) {
            d.phase = (d.phase + phase.into()).normalize();
        } else {
            panic!("Vertex not found");
        }
    }

    fn set_vertex_type(&mut self, v: V, ty: VType) {
        self.vdata.get_mut(&v).expect("Vertex not found").ty = ty;
    }

    fn vertex_data(&self, v: V) -> VData {
        *self.vdata.get(&v).expect("Vertex not found")
    }

    fn vertex_type(&self, v: V) -> VType {
        self.vdata.get(&v).expect("Vertex not found").ty
    }

    fn set_edge_type(&mut self, s: V, t: V, ety: EType) {
        *self
            .edata
            .get_mut(&s)
            .expect("Source vertex not found")
            .get_mut(&t)
            .expect("Edge not found") = ety;
        *self
            .edata
            .get_mut(&t)
            .expect("Target vertex not found")
            .get_mut(&s)
            .expect("Edge not found") = ety;
    }

    fn edge_type_opt(&self, s: V, t: V) -> Option<EType> {
        self.edata
            .get(&s)
            .expect("Source vertex not found")
            .get(&t)
            .copied()
    }

    fn set_coord(&mut self, v: V, coord: impl Into<Coord>) {
        let coord = coord.into();
        let d = self.vdata.get_mut(&v).expect("Vertex not found");
        d.qubit = coord.x;
        d.row = coord.y;
    }

    fn coord(&self, v: V) -> Coord {
        let d = self.vdata.get(&v).expect("Vertex not found");
        Coord::new(d.qubit, d.row)
    }

    fn set_qubit(&mut self, v: V, qubit: i32) {
        self.vdata.get_mut(&v).expect("Vertex not found").qubit = qubit;
    }

    fn qubit(&self, v: V) -> i32 {
        self.vdata.get(&v).expect("Vertex not found").qubit
    }

    fn set_row(&mut self, v: V, row: i32) {
        self.vdata.get_mut(&v).expect("Vertex not found").row = row;
    }

    fn row(&self, v: V) -> i32 {
        self.vdata.get(&v).expect("Vertex not found").row
    }

    fn neighbors(&self, v: V) -> NeighborIter {
        NeighborIter::Hash(self.edata.get(&v).expect("Vertex not found").keys())
    }

    fn incident_edges(&self, v: V) -> IncidentEdgeIter {
        IncidentEdgeIter::Hash(self.edata.get(&v).expect("Vertex not found").iter())
    }

    fn degree(&self, v: V) -> usize {
        self.edata.get(&v).expect("Vertex not found").len()
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
        for (&v0, tab) in self.edata.iter() {
            for (&v1, &et) in tab.iter() {
                if f(v0, v1, et) {
                    return Some((v0, v1, et));
                }
            }
        }

        None
    }

    fn find_vertex<F>(&self, f: F) -> Option<V>
    where
        F: Fn(V) -> bool,
    {
        self.vdata.keys().find(|&&v| f(v)).copied()
    }

    fn contains_vertex(&self, v: V) -> bool {
        self.vdata.contains_key(&v)
    }
}

impl serde::Serialize for Graph {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let jg = JsonGraph::from_graph(self)
            .map_err(|e| S::Error::custom(format!("Failed to convert graph to JSON graph: {e}")))?;
        jg.serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for Graph {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let jg = JsonGraph::deserialize(deserializer)?;
        let g = jg
            .to_graph()
            .map_err(|e| D::Error::custom(format!("Failed to convert JSON graph to graph: {e}")))?;
        Ok(g)
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
        let mut vs = Vec::from_iter(g.vertices());
        vs.sort();
        expected_vs.sort();
        assert_eq!(expected_vs, vs);
    }

    #[test]
    fn edge_iterator() {
        let (mut g, vs) = simple_graph();
        g.set_edge_type(vs[1], vs[3], EType::H);

        let mut edges = Vec::from_iter(g.edges());
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
