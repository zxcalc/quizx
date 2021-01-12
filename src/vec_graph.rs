use crate::graph::*;
use crate::scalar::*;
use num::rational::Rational;
use std::iter::FromIterator;

pub type VTab<T> = Vec<Option<T>>;

#[derive(Debug,Clone,PartialEq)]
pub struct Graph {
    vdata: VTab<VData>,
    edata: VTab<Vec<(V,EType)>>,
    inputs: Vec<V>,
    outputs: Vec<V>,
    numv: usize,
    nume: usize,
    pub scalar: ScalarN,
}

// pub struct VertexIter<'a> {
//     inner: std::collections::hash_map::Keys<'a,V,VData>,
// }

// impl<'a> Iterator for VertexIter<'a> {
//     type Item = V;
//     fn next(&mut self) -> Option<Self::Item> {
//         self.inner.next().map(|x| *x)
//     }
// }

// pub struct EdgeIter<'a> {
//     outer: std::collections::hash_map::Iter<'a,V,VTab<EType>>,
//     inner: Option<(V, std::collections::hash_map::Iter<'a,V,EType>)>,
// }

// impl<'a> Iterator for EdgeIter<'a> {
//     /// Iterate over the edges in a graph. An edge is returned as a triple
//     /// (s: V, t: V, ety: EType), where we enforce s <= t to avoid double-
//     /// counting edges.
//     type Item = (V,V,EType);

//     fn next(&mut self) -> Option<Self::Item> {
//        match &mut self.inner {
//            Some((s, iter)) =>
//                match iter.next() {
//                    Some((t,ety)) => if *s <= *t { Some((*s,*t,*ety)) } else { self.next() }
//                    None => match self.outer.next() {
//                        Some((k,v)) => { self.inner = Some((*k,v.iter())); self.next() }
//                        None => None
//                    }
//                }
//            None => None
//        }
//     }
// }

impl Graph {
    pub fn new() -> Graph {
        Graph {
            vdata: Vec::new(),
            edata: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            numv: 0,
            nume: 0,
            scalar: Scalar::one(),
        }
    }

    fn index<U>(nhd: &Vec<(V,U)>, v: V) -> Option<usize> {
        nhd.iter().position(|&(v0,_)| v == v0)
    }

    fn value<U: Copy>(nhd: &Vec<(V,U)>, v: V) -> Option<U> {
        for (v0,u) in nhd.iter() {
            if v == *v0 { return Some(*u); }
        }
        None
    }

    /// Removes vertex 't' from the adjacency map of 's'. This private method
    /// is used by remove_edge and remove_vertex to make the latter slightly
    /// more efficient.
    fn remove_half_edge(&mut self, s: V, t: V) {
        if let Some(Some(nhd)) = self.edata.get_mut(s) {
            let i = Graph::index(&nhd,t).expect("Target vertex not found");
            nhd.swap_remove(i);
        } else {
            panic!("Source vertex not found");
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

impl IsGraph for Graph {
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

    fn inputs(&self) -> &Vec<V> { &self.inputs }
    fn set_inputs(&mut self, inputs: Vec<V>) { self.inputs = inputs; }
    fn outputs(&self) -> &Vec<V> { &self.outputs }
    fn set_outputs(&mut self, outputs: Vec<V>) { self.outputs = outputs; }

    // fn vertices(&self) -> VertexIter {
    //     VertexIter { inner: self.vdata.keys() }
    // }

    // fn edges(&self) -> EdgeIter {
    //     let mut outer = self.edata.iter();
    //     let inner = match outer.next() {
    //         Some((s, h)) => Some((*s, h.iter())),
    //         None => None,
    //     };
    //     EdgeIter { outer, inner }
    // }

    fn add_vertex(&mut self, ty: VType) -> V {
        self.add_vertex_with_data(VData { ty, phase: Rational::new(0,1), qubit: 0, row: 0 })
    }

    fn add_vertex_with_data(&mut self, d: VData) -> V {
        self.numv += 1;
        self.vdata.push(Some(d));
        self.edata.push(Some(Vec::new()));
        self.vdata.len() - 1
    }

    fn remove_vertex(&mut self, v: V) {
        self.numv -= 1;

        for v1 in Vec::from_iter(self.neighbors(v)) {
            self.nume -= 1;
            self.remove_half_edge(v1,v);
        }

        self.vdata[v] = None;
        self.edata[v] = None;
    }

    fn add_edge_with_type(&mut self, s: V, t: V, ety: EType) {
        self.nume += 1;

        if let Some(Some(nhd)) = self.edata.get_mut(s) {
            nhd.push((t,ety));
        } else {
            panic!("Source vertex not found");
        }

        if let Some(Some(nhd)) = self.edata.get_mut(t) {
            nhd.push((s,ety));
        } else {
            panic!("Target vertex not found");
        }
    }


    fn remove_edge(&mut self, s: V, t: V) {
        self.nume -= 1;
        self.remove_half_edge(s,t);
        self.remove_half_edge(t,s);
    }

    fn set_phase(&mut self, v: V, phase: Rational) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.phase = phase.mod2();
        } else {
            panic!("Vertex not found");
        }
    }

    fn phase(&self, v: V) -> Rational {
        self.vdata[v]
            .expect("Vertex not found")
            .phase
    }

    fn add_to_phase(&mut self, v: V, phase: Rational) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.phase = (d.phase + phase).mod2();
        } else {
            panic!("Vertex not found");
        }
    }

    fn set_vertex_type(&mut self, v: V, ty: VType) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.ty = ty;
        } else {
            panic!("Vertex not found");
        }
    }

    fn vertex_type(&self, v: V) -> VType {
        self.vdata[v]
            .expect("Vertex not found")
            .ty
    }

    fn set_edge_type(&mut self, s: V, t: V, ety: EType) {
        if let Some(Some(nhd)) = self.edata.get_mut(s) {
            let i = Graph::index(&nhd, t).expect("Edge not found");
            nhd[i] = (t, ety);
        } else {
            panic!("Source vertex not found");
        }

        if let Some(Some(nhd)) = self.edata.get_mut(t) {
            let i = Graph::index(&nhd, s).expect("Edge not found");
            nhd[i] = (s, ety);
        } else {
            panic!("Target vertex not found");
        }
    }

    fn edge_type_opt(&self, s: V, t: V) -> Option<EType> {
        if let Some(Some(nhd)) = self.edata.get(s) {
            Graph::value(&nhd, t)
        } else {
            None
        }
    }


    fn set_coord(&mut self, v: V, coord: (i32,i32)) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.qubit = coord.0;
            d.row = coord.1;
        } else {
            panic!("Vertex not found")
        }
    }

    fn coord(&mut self, v: V) -> (i32,i32) {
        let d = self.vdata[v].expect("Vertex not found");
        (d.qubit, d.row)
    }

    fn set_qubit(&mut self, v: V, qubit: i32) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.qubit = qubit;
        } else {
            panic!("Vertex not found")
        }
    }

    fn qubit(&mut self, v: V) -> i32 {
        self.vdata[v]
            .expect("Vertex not found").qubit
    }

    fn set_row(&mut self, v: V, row: i32) {
        if let Some(Some(d)) = self.vdata.get_mut(v) {
            d.row = row;
        } else {
            panic!("Vertex not found")
        }
    }

    fn row(&mut self, v: V) -> i32 {
        self.vdata[v]
            .expect("Vertex not found").row
    }

    fn neighbors(&self, v: V) -> NeighborIter {
        if let Some(Some(nhd)) = self.edata.get(v) {
            NeighborIter::Vec(nhd.iter())
        } else {
            panic!("Vertex not found")
        }
    }

    fn incident_edges(&self, v: V) -> IncidentEdgeIter {
        if let Some(Some(nhd)) = self.edata.get(v) {
            IncidentEdgeIter::Vec(nhd.iter())
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

    fn scalar(&mut self) -> &mut ScalarN { &mut self.scalar }

    fn find_edge<F>(&self, f: F) -> Option<(V,V,EType)>
        where F : Fn(V,V,EType) -> bool
    {
        for (v0, nhd0) in self.edata.iter().enumerate() {
            if let Some(nhd) = nhd0 {
                for &(v1, et) in nhd.iter() {
                    if v0 <= v1 && f(v0,v1,et) { return Some((v0,v1,et)); }
                }
            };
        }

        None
    }

    fn find_vertex<F>(&self, f: F) -> Option<V>
        where F : Fn(V) -> bool
    {
        for (v, d) in self.vdata.iter().enumerate() {
            if d.is_some() && f(v) { return Some(v); }
        }

        None
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

    fn simple_graph() -> (Graph,Vec<V>) {
        let mut g = Graph::new();
        let vs = vec![
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::X),
            g.add_vertex(VType::X),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B)];
        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[3]);
        g.add_edge(vs[2], vs[4]);
        g.add_edge(vs[2], vs[5]);
        g.add_edge(vs[3], vs[4]);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[4], vs[6]);
        g.add_edge(vs[5], vs[7]);
        (g,vs)
    }

    #[test]
    fn create_simple_graph() {
        let (g,_) = simple_graph();
        assert_eq!(g.num_vertices(), 8);
        assert_eq!(g.num_edges(), 8);
    }

    #[test]
    fn clone_graph() {
       let (g,_) = simple_graph();
       let h = g.clone();
       assert!(g.num_vertices() == h.num_vertices());
       assert!(g.num_edges() == h.num_edges());
       assert!(g == h);
    }

    // #[test]
    // fn nhd() {
    //     let (g,vs) = simple_graph();
    //     let mut ws = Vec::from_iter(g.nhd(vs[2]));
    //     ws.sort();
    //     let mut expected_ws = vec![
    //         (&vs[0], &EType::N),
    //         (&vs[4], &EType::N),
    //         (&vs[5], &EType::N)];
    //     expected_ws.sort();
    //     assert_eq!(expected_ws, ws);
    // }

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
            g.add_vertex(VType::B)];
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

