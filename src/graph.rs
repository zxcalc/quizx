use std::collections::HashMap;
use num_rational::Rational;

pub type V = u32;
pub type E = u32;

#[derive(Debug,Copy,Clone)]
pub enum VType {
    B, // Boundary
    Z, // Z-spider
    X, // X-spider
    H, // H-box
}

#[derive(Debug,Copy,Clone)]
pub struct VData {
    ty: VType,
    phase: Rational,
    qubit: i32,
    row: i32,
}

#[derive(Debug,Copy,Clone)]
pub enum EType {
    N, // normal edge
    H, // hadamard edge
}

pub type AdjIter<'a> = std::collections::hash_map::Keys<'a,V,EType>;

#[derive(Debug)]
pub struct Graph {
    vdata: HashMap<V,VData>,
    edata: HashMap<V,HashMap<V,EType>>,
    inputs: Vec<V>,
    outputs: Vec<V>,
    numv: usize,
    nume: usize,
    freshv: V,
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

    pub fn add_edge(&mut self, s: V, t: V, et: EType) {
        self.nume += 1;

        self.edata.get_mut(&s)
            .expect("Source vertex not found")
            .insert(t, et);
        self.edata.get_mut(&t)
            .expect("Target vertex not found")
            .insert(s, et);
    }

    pub fn add_edge_smart(&mut self, s: V, t: V, et: EType) {
        // TODO: scalars
        if let Some(et0) = self.edata.get(&s).and_then(|x| x.get(&t)) {
            let st = self.vdata.get(&s).expect("Source vertex not found").ty;
            let tt = self.vdata.get(&t).expect("Target vertex not found").ty;
            match (st, tt) {
                (VType::Z, VType::Z) | (VType::X, VType::X) => {
                    match (et0, et) {
                        (EType::N, EType::N) => {}
                        (EType::H, EType::H) => {}
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
                }
                _ => panic!("Parallel edges only supported between Z and X vertices")
            }
        } else {
            self.add_edge(s, t, et);
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

    pub fn set_edge_type(&mut self, s: V, t: V, et: EType) {
        *(self.edata.get_mut(&s)
            .expect("Source vertex not found")
            .get_mut(&t)
            .expect("Edge not found")) = et;
        *(self.edata.get_mut(&t)
            .expect("Target vertex not found")
            .get_mut(&s)
            .expect("Edge not found")) = et;
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
    use super::*;

    #[test]
    fn create_empty_graph() {
        let g = Graph::new();
        assert_eq!(g.num_vertices(), 0);
        assert_eq!(g.num_edges(), 0);
    }
}

