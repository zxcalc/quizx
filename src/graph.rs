use crate::scalar::Scalar;
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
            EIter::Vec(_,outer,inner)  => {
                match inner {
                    Some((v, inner1)) => match inner1.next() {
                        Some((v1,et)) =>
                            if *v <= *v1 {
                                Some((*v,*v1,*et))
                            } else {
                                self.next()
                            },
                        None => {
                            *inner = None;
                            self.next()
                        }
                    },
                    None => match outer.next() {
                        Some((v, Some(tab))) => {
                            *inner = Some((v, tab.iter()));
                            self.next()
                        },
                        Some((_, None)) => self.next(),
                        None => None,
                    }
                }
            },

            EIter::Hash(_, outer, inner) => {
                match inner {
                    Some((v, inner1)) => match inner1.next() {
                        Some((v1,et)) =>
                            if *v <= *v1 {
                                Some((*v,*v1,*et))
                            } else {
                                self.next()
                            },
                        None => {
                            *inner = None;
                            self.next()
                        }
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

pub trait IsGraph {
    fn num_vertices(&self) -> usize;
    fn num_edges(&self) -> usize;
    fn vertices(&self) -> VIter;
    fn edges(&self) -> EIter;
    fn inputs(&self) -> &Vec<V>;
    fn set_inputs(&mut self, inputs: Vec<V>);
    fn outputs(&self) -> &Vec<V>;
    fn set_outputs(&mut self, outputs: Vec<V>);
    fn add_vertex(&mut self, ty: VType) -> V;
    fn add_vertex_with_data(&mut self, d: VData) -> V;
    fn remove_vertex(&mut self, v: V);
    fn add_edge_with_type(&mut self, s: V, t: V, ety: EType);
    fn remove_edge(&mut self, s: V, t: V);
    fn add_edge_smart(&mut self, s: V, t: V, ety: EType);
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
    fn qubit(&mut self, v: V) -> i32;
    fn set_row(&mut self, v: V, row: i32);
    fn row(&mut self, v: V) -> i32;
    fn neighbors(&self, v: V) -> NeighborIter;
    fn incident_edges(&self, v: V) -> IncidentEdgeIter;
    fn degree(&self, v: V) -> usize;
    fn scalar(&mut self) -> &mut Scalar;
    fn find_edge<F>(&self, f: F) -> Option<(V,V,EType)>
        where F : Fn(V,V,EType) -> bool;
    fn find_vertex<F>(&self, f: F) -> Option<V>
        where F : Fn(V) -> bool;

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
        match self.edge_type(v0, v1) {
            EType::N => self.set_edge_type(v0,v1,EType::H),
            EType::H => self.set_edge_type(v0,v1,EType::N),
        }
    }

    fn vertex_vec(&self) -> Vec<V> { self.vertices().collect() }
    fn edge_vec(&self) -> Vec<(V,V,EType)> { self.edges().collect() }
    fn neighbor_vec(&self, v: V) -> Vec<V> { self.neighbors(v).collect() }
    fn incident_edge_vec(&self, v: V) -> Vec<(V,EType)> { self.incident_edges(v).collect() }

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
}

