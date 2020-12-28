use std::collections::HashMap;
use num_rational::Rational;

pub type V = u32;
pub type E = u32;

#[derive(Debug,Copy,Clone,PartialEq,Eq,PartialOrd,Ord)]
pub enum VType {
    B, // Boundary
    Z, // Z-spider
    X, // X-spider
    H, // H-box
}

#[derive(Debug,Copy,Clone,PartialEq,Eq)]
pub struct VData {
    ty: VType,
    phase: Rational,
    qubit: i32,
    row: i32,
}

#[derive(Debug,Copy,Clone,PartialEq,Eq,PartialOrd,Ord)]
pub enum EType {
    N, // normal edge
    H, // hadamard edge
}

pub type AdjIter<'a> = std::collections::hash_map::Keys<'a,V,EType>;

#[derive(Debug,Clone,PartialEq,Eq)]
pub struct Graph {
    vdata: HashMap<V,VData>,
    edata: HashMap<V,HashMap<V,EType>>,
    inputs: Vec<V>,
    outputs: Vec<V>,
    numv: usize,
    nume: usize,
    freshv: V,
}

pub struct EdgeIter<'a> {
    outer: std::collections::hash_map::Iter<'a,V,HashMap<V,EType>>,
    inner: Option<(V, std::collections::hash_map::Iter<'a,V,EType>)>,
}

impl Iterator for EdgeIter<'_> {
    /// Iterate over the edges in a graph. An edge is returned as a triple
    /// (s: V, t: V, ety: EType), where we enforce s <= t to avoid double-
    /// counting edges.
    type Item = (V,V,EType);

    fn next(&mut self) -> Option<Self::Item> {
       match &mut self.inner {
           Some((s, iter)) =>
               match iter.next() {
                   Some((t,ety)) => if *s <= *t { Some((*s,*t,*ety)) } else { self.next() }
                   None => match self.outer.next() {
                       Some((k,v)) => { self.inner = Some((*k,v.iter())); self.next() }
                       None => None
                   }
               }
           None => None
       }
    }
}

impl Graph {
    pub fn new() -> Graph {
        Graph {
            vdata: HashMap::new(),
            edata: HashMap::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            numv: 0,
            nume: 0,
            freshv: 0,
        }
    }

    pub fn num_vertices(&self) -> usize {
        self.numv
    }

    pub fn num_edges(&self) -> usize {
        self.nume
    }

    pub fn edges(&self) -> EdgeIter {
        let mut outer = self.edata.iter();
        let inner = match outer.next() {
            Some((s, h)) => Some((*s, h.iter())),
            None => None,
        };
        EdgeIter { outer, inner }
    }

    pub fn add_vertex(&mut self, ty: VType) -> V {
        self.add_vertex_with_data(VData { ty, phase: Rational::new(0,1), qubit: 0, row: 0 })
    }

    pub fn add_vertex_with_data(&mut self, d: VData) -> V {
        let v = self.freshv;
        self.freshv += 1;
        self.numv += 1;
        self.vdata.insert(v, d);
        self.edata.insert(v, HashMap::new());
        v
    }

    pub fn add_edge(&mut self, s: V, t: V) {
        self.add_edge_with_type(s, t, EType::N);
    }

    pub fn add_edge_with_type(&mut self, s: V, t: V, ety: EType) {
        self.nume += 1;

        self.edata.get_mut(&s)
            .expect("Source vertex not found")
            .insert(t, ety);
        self.edata.get_mut(&t)
            .expect("Target vertex not found")
            .insert(s, ety);
    }

    pub fn remove_edge(&mut self, s: V, t: V) {
        self.nume -= 1;

        self.edata.get_mut(&s)
            .expect("Source vertex not found")
            .remove(&t);
        self.edata.get_mut(&t)
            .expect("Target vertex not found")
            .remove(&s);
    }

    pub fn add_edge_smart(&mut self, s: V, t: V, ety: EType) {
        // TODO: scalars
        if let Some(ety0) = self.edata.get(&s).and_then(|x| x.get(&t)) {
            let st = self.vdata.get(&s).expect("Source vertex not found").ty;
            let tt = self.vdata.get(&t).expect("Target vertex not found").ty;
            match (st, tt) {
                (VType::Z, VType::Z) | (VType::X, VType::X) => {
                    match (ety0, ety) {
                        (EType::N, EType::N) => {} // ignore new edge
                        (EType::H, EType::H) => {
                            self.remove_edge(s, t);
                        }
                        (EType::H, EType::N) => {
                            self.set_edge_type(s, t, EType::N);
                            self.add_to_phase(s, Rational::new(1,1));
                        }
                        (EType::N, EType::H) => {
                            self.add_to_phase(s, Rational::new(1,1));
                        }
                    }
                }
                (VType::Z, VType::X) | (VType::X, VType::Z) => {
                    match (ety0, ety) {
                        (EType::N, EType::N) => {
                            self.remove_edge(s, t);
                        }
                        (EType::N, EType::H) => {
                            self.set_edge_type(s, t, EType::H);
                            self.add_to_phase(s, Rational::new(1,1));
                        }
                        (EType::H, EType::N) => {
                            self.add_to_phase(s, Rational::new(1,1));
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

    pub fn set_phase(&mut self, v: V, phase: Rational) {
        self.vdata.get_mut(&v)
            .expect("Vertex not found")
            .phase = phase;
    }

    pub fn phase(&self, v: V) -> Rational {
        self.vdata.get(&v)
            .expect("Vertex not found")
            .phase
    }

    pub fn add_to_phase(&mut self, v: V, phase: Rational) {
        self.vdata.get_mut(&v)
            .expect("Vertex not found")
            .phase += phase;
    }

    pub fn set_vertex_type(&mut self, v: V, ty: VType) {
        self.vdata.get_mut(&v)
            .expect("Vertex not found")
            .ty = ty;
    }

    pub fn vertex_type(&self, v: V) -> VType {
        self.vdata.get(&v)
            .expect("Vertex not found")
            .ty
    }

    pub fn set_edge_type(&mut self, s: V, t: V, ety: EType) {
        *self.edata.get_mut(&s)
            .expect("Source vertex not found")
            .get_mut(&t)
            .expect("Edge not found") = ety;
        *self.edata.get_mut(&t)
            .expect("Target vertex not found")
            .get_mut(&s)
            .expect("Edge not found") = ety;
    }

    pub fn edge_type(&self, s: V, t: V) -> Option<EType> {
        self.edata.get(&s).and_then(|x| x.get(&t)).copied()
    }

    pub fn set_coord(&mut self, v: V, coord: (i32,i32)) {
        let d = self.vdata.get_mut(&v).expect("Vertex not found");
        d.qubit = coord.0;
        d.row = coord.1;
    }

    pub fn coord(&mut self, v: V) -> (i32,i32) {
        let d = self.vdata.get(&v).expect("Vertex not found");
        (d.qubit, d.row)
    }

    pub fn set_qubit(&mut self, v: V, qubit: i32) {
        self.vdata.get_mut(&v)
            .expect("Vertex not found").qubit = qubit;
    }

    pub fn qubit(&mut self, v: V) -> i32 {
        self.vdata.get(&v)
            .expect("Vertex not found").qubit
    }

    pub fn set_row(&mut self, v: V, row: i32) {
        self.vdata.get_mut(&v)
            .expect("Vertex not found").row = row;
    }

    pub fn row(&mut self, v: V) -> i32 {
        self.vdata.get(&v)
            .expect("Vertex not found").row
    }

    pub fn connected(&self, v0: V, v1: V) -> bool {
        self.edata.get(&v0)
            .expect("Vertex not found")
            .contains_key(&v1)
    }

    pub fn neighbors(&self, v: V) -> AdjIter {
        self.edata.get(&v)
            .expect("Vertex not found")
            .keys()
    }

    pub fn nhd(&self, v: V) -> &HashMap<V,EType> {
        self.edata.get(&v)
            .expect("Vertex not found")
    }
}

#[cfg(test)]
mod tests {
    use std::iter::FromIterator;
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
        assert_eq!(g.num_vertices(), 4);
        assert_eq!(g.num_edges(), 2);

        let mut h = g.clone();
        h.add_edge_smart(vs[1], vs[2], EType::H);
        h.add_edge_smart(vs[1], vs[2], EType::H);
        assert_eq!(g.num_vertices(), 4);
        assert_eq!(g.num_edges(), 3);
    }
}

