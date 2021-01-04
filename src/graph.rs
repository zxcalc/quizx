use num::rational::Rational;
use crate::scalar::Scalar;

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

pub trait IsGraph {
    fn num_vertices(&self) -> usize;
    fn num_edges(&self) -> usize;
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
    fn neighbors(&self, v: V) -> Vec<V>;
    fn incident_edges(&self, v: V) -> Vec<(V,EType)>;
    fn degree(&self, v: V) -> usize;
    fn scalar(&self) -> &Scalar;
    fn set_scalar(&mut self, s: Scalar);

    fn add_edge(&mut self, s: V, t: V) {
        self.add_edge_with_type(s, t, EType::N);
    }

    fn edge_type(&self, s: V, t: V) -> EType {
        self.edge_type_opt(s,t).expect("Edge not found")
    }

    fn connected(&self, v0: V, v1: V) -> bool {
        self.edge_type_opt(v0, v1).is_some()
    }

}

