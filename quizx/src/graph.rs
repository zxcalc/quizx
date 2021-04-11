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

use crate::scalar::*;
use num::rational::Rational;
use std::iter::FromIterator;

pub type V = usize;

#[derive(Debug,Copy,Clone,PartialEq,Eq,PartialOrd,Ord)]
pub enum VType {
    B, // Boundary
    Z, // Z-spider
    X, // X-spider
    H, // H-box
}

#[derive(Debug,Copy,Clone,PartialEq,Eq)]
pub struct VData {
    pub ty: VType,
    pub phase: Rational,
    pub qubit: i32,
    pub row: i32,
}

#[derive(Debug,Copy,Clone,PartialEq,Eq,PartialOrd,Ord)]
pub enum EType {
    N, // normal edge
    H, // hadamard edge
}

impl EType {
    pub fn opposite(&self) -> EType {
        match self {
            EType::N => EType::H,
            EType::H => EType::N,
        }
    }
}

pub enum VIter<'a> {
    Vec(usize,std::iter::Enumerate<std::slice::Iter<'a,Option<VData>>>),
    Hash(std::collections::hash_map::Keys<'a,V,VData>)
}

impl<'a> Iterator for VIter<'a> {
    type Item = V;
    fn next(&mut self) -> Option<V> {
        match self {
            VIter::Vec(_,inner)  => {
                match inner.next() {
                    Some((v, Some(_))) => Some(v),
                    Some((_, None)) => self.next(),
                    None => None
                }
            },
            VIter::Hash(inner) => inner.next().map(|&v| v)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            VIter::Vec(sz,_)  => *sz,
            VIter::Hash(inner) => inner.len(),
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for VIter<'a> {}

pub enum EIter<'a> {
    Vec(usize,
        std::iter::Enumerate<std::slice::Iter<'a,Option<Vec<(V,EType)>>>>,
        Option<(V,std::slice::Iter<'a,(V,EType)>)>),
    Hash(usize,
         std::collections::hash_map::Iter<'a,V,rustc_hash::FxHashMap<V,EType>>,
         Option<(V,std::collections::hash_map::Iter<'a,V,EType>)>)
}

impl<'a> Iterator for EIter<'a> {
    type Item = (V,V,EType);
    fn next(&mut self) -> Option<(V,V,EType)> {
        match self {
            EIter::Vec(_,outer,inner) => match inner {
                Some((v, inner1)) => match inner1.next() {
                    Some((v1,et)) =>
                        if *v <= *v1 { Some((*v,*v1,*et)) }
                        else { self.next() },
                    None => { *inner = None; self.next() }
                },
                None => match outer.next() {
                    Some((v, Some(tab))) => {
                        *inner = Some((v, tab.iter()));
                        self.next()
                    },
                    Some((_, None)) => self.next(),
                    None => None,
                }
            },

            EIter::Hash(_, outer, inner) => match inner {
                Some((v, inner1)) => match inner1.next() {
                    Some((v1,et)) =>
                        if *v <= *v1 { Some((*v,*v1,*et)) }
                        else { self.next() },
                    None => { *inner = None; self.next() }
                },
                None => match outer.next() {
                    Some((v, tab)) => {
                        *inner = Some((*v, tab.iter()));
                        self.next()
                    },
                    None => None
                }
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            EIter::Vec(sz, ..)  => *sz,
            EIter::Hash(sz, ..) => *sz,
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for EIter<'a> {}

pub enum NeighborIter<'a> {
    Vec(std::slice::Iter<'a,(V,EType)>),
    Hash(std::collections::hash_map::Keys<'a,V,EType>)
}

impl<'a> Iterator for NeighborIter<'a> {
    type Item = V;
    fn next(&mut self) -> Option<V> {
        match self {
            NeighborIter::Vec(inner)  => inner.next().map(|&(v,_)| v),
            NeighborIter::Hash(inner) => inner.next().map(|&v| v)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            NeighborIter::Vec(inner)  => inner.len(),
            NeighborIter::Hash(inner) => inner.len(),
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for NeighborIter<'a> {}


pub enum IncidentEdgeIter<'a> {
    Vec(std::slice::Iter<'a,(V,EType)>),
    Hash(std::collections::hash_map::Iter<'a,V,EType>)
}

impl<'a> Iterator for IncidentEdgeIter<'a> {
    type Item = (V,EType);
    fn next(&mut self) -> Option<(V,EType)> {
        match self {
            IncidentEdgeIter::Vec(inner)  => inner.next().map(|&x| x),
            IncidentEdgeIter::Hash(inner) => inner.next().map(|(&v,&et)| (v,et))
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            IncidentEdgeIter::Vec(inner)  => inner.len(),
            IncidentEdgeIter::Hash(inner) => inner.len(),
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for IncidentEdgeIter<'a> {}

pub trait GraphLike: Clone + std::fmt::Debug {
    /// Initialise a new empty graph
    fn new() -> Self;

    /// Next fresh vertex index
    fn vindex(&self) -> V;

    /// Number of vertices
    fn num_vertices(&self) -> usize;

    /// Number of edges
    fn num_edges(&self) -> usize;

    /// Get iterator over all vertices
    fn vertices(&self) -> VIter;

    /// Get iterator over all edges
    ///
    /// An "edge" is a triple (s, t, edge_type), where s <= t.
    fn edges(&self) -> EIter;

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

    fn set_phase(&mut self, v: V, phase: Rational);
    fn phase(&self, v: V) -> Rational;
    fn add_to_phase(&mut self, v: V, phase: Rational);
    fn set_vertex_type(&mut self, v: V, ty: VType);
    fn vertex_type(&self, v: V) -> VType;
    fn set_edge_type(&mut self, s: V, t: V, ety: EType);
    fn edge_type_opt(&self, s: V, t: V) -> Option<EType>;
    fn set_coord(&mut self, v: V, coord: (i32,i32));
    fn coord(&mut self, v: V) -> (i32,i32);
    fn set_qubit(&mut self, v: V, qubit: i32);
    fn qubit(&self, v: V) -> i32;
    fn set_row(&mut self, v: V, row: i32);
    fn row(&self, v: V) -> i32;
    fn neighbors(&self, v: V) -> NeighborIter;
    fn incident_edges(&self, v: V) -> IncidentEdgeIter;
    fn degree(&self, v: V) -> usize;
    fn scalar(&self) -> &ScalarN;
    fn scalar_mut(&mut self) -> &mut ScalarN;
    fn find_edge<F>(&self, f: F) -> Option<(V,V,EType)>
        where F : Fn(V,V,EType) -> bool;
    fn find_vertex<F>(&self, f: F) -> Option<V>
        where F : Fn(V) -> bool;
    fn contains_vertex(&self, v: V) -> bool;

    fn add_edge(&mut self, s: V, t: V) {
        self.add_edge_with_type(s, t, EType::N);
    }

    fn edge_type(&self, s: V, t: V) -> EType {
        self.edge_type_opt(s,t).expect("Edge not found")
    }

    fn connected(&self, v0: V, v1: V) -> bool {
        self.edge_type_opt(v0, v1).is_some()
    }

    fn toggle_edge_type(&mut self, v0: V, v1: V) {
        self.set_edge_type(v0, v1, self.edge_type(v0, v1).opposite());
    }

    fn vertex_vec(&self) -> Vec<V> { self.vertices().collect() }
    fn edge_vec(&self) -> Vec<(V,V,EType)> { self.edges().collect() }
    fn neighbor_vec(&self, v: V) -> Vec<V> { self.neighbors(v).collect() }
    fn incident_edge_vec(&self, v: V) -> Vec<(V,EType)> { self.incident_edges(v).collect() }

    /// Convert all X spiders to Z with the colour-change rule
    fn x_to_z(&mut self) {
        for v in Vec::from_iter(self.vertices()) {
            if self.vertex_type(v) == VType::X {
                self.set_vertex_type(v, VType::Z);
                for w in Vec::from_iter(self.neighbors(v)) {
                    self.toggle_edge_type(v,w);
                }
            }
        }
    }

    fn add_vertex_with_phase(&mut self, ty: VType, phase: Rational) -> V {
        let v = self.add_vertex(ty);
        self.set_phase(v, phase);
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
                   self.add_to_phase(s, Rational::new(1,1));
                   self.scalar_mut().mul_sqrt2_pow(-1);
               }
           } else {
               panic!("Self-loops only supported on Z and X nodes");
           }
        } else if let Some(ety0) = self.edge_type_opt(s,t) {
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
                            self.add_to_phase(s, Rational::new(1,1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::N, EType::H) => {
                            self.add_to_phase(s, Rational::new(1,1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
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
                            self.add_to_phase(s, Rational::new(1,1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::H, EType::N) => {
                            self.add_to_phase(s, Rational::new(1,1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::H, EType::H) => {} // ignore new edge
                    }
                }
                _ => panic!("Parallel edges only supported between Z and X vertices")
            }
        } else {
            self.add_edge_with_type(s, t, ety);
        }
    }

    /// Return a graphviz-friendly string representation of the graph
    fn to_dot(&self) -> String {
        let mut dot = String::from("graph {\n");
        for v in self.vertices() {
            let t = self.vertex_type(v);
            let p = self.phase(v);
            dot += &format!("  {} [color={}, label=\"{}\"", v,
                            match t {
                                VType::B => "black",
                                VType::Z => "green",
                                VType::X => "red",
                                VType::H => "yellow",
                            },
                            if self.inputs().contains(&v) { format!("{}:i", v) }
                            else if self.outputs().contains(&v) { format!("{}:o", v) }
                            else if !p.is_zero() { format!("{}:{}", v, p) }
                            else { format!("{}", v) }
                           );
            let q = self.qubit(v);
            let r = self.row(v);
            if q != 0 || r != 0 {
                dot += &format!(", pos=\"{},{}!\"", q, r);
            }
            dot += "]\n";
        }

        dot += "\n";

        for (s, t, ty) in self.edges() {
            dot += &format!("  {} -- {}", s, t);
            if ty == EType::H { dot += " [color=blue]"; }
            dot += "\n";
        }

        dot += "}\n";

        dot
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec_graph::Graph;
    use crate::tensor::ToTensor;
    #[test]
    fn smart_edges() {
       let mut g = Graph::new();
       g.add_vertex(VType::B);
       g.add_vertex(VType::Z);
       g.add_vertex(VType::Z);
       g.add_vertex(VType::X);
       g.add_vertex(VType::B);
       g.add_edge_smart(0,1,EType::N);
       g.add_edge_smart(1,2,EType::N);
       g.add_edge_smart(2,3,EType::N);
       g.add_edge_smart(1,3,EType::N);
       g.add_edge_smart(3,4,EType::N);
       g.set_inputs(vec![0]);
       g.set_outputs(vec![4]);

       let mut h = Graph::new();
       h.add_vertex(VType::B);
       h.add_vertex(VType::Z);
       h.add_vertex(VType::X);
       h.add_vertex(VType::B);
       h.add_edge_smart(0,1,EType::N);
       h.add_edge_smart(1,2,EType::N);
       h.add_edge_smart(1,2,EType::N);
       h.add_edge_smart(2,3,EType::N);
       h.set_inputs(vec![0]);
       h.set_outputs(vec![3]);

       let tg = g.to_tensor4();
       let th = h.to_tensor4();
       println!("\n\ntg =\n{}", tg);
       println!("\n\nth =\n{}", th);
       assert_eq!(tg, th);
    }
}
