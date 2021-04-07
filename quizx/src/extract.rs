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

use crate::circuit::*;
use crate::gate::*;
use crate::graph::*;
// use crate::tensor::*;
use crate::linalg::*;
use std::fmt;
use num::{Rational, Zero};
use rustc_hash::FxHashSet;
use crate::basic_rules::{boundary_pivot, remove_id};

/// Extraction couldn't finish. Returns a message, a
/// partially-extracted circuit, and the remainder of
/// the graph.
pub struct ExtractError<G: GraphLike>(pub String, pub Circuit, pub G);

impl<G: GraphLike> fmt::Display for ExtractError<G> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<G: GraphLike> fmt::Debug for ExtractError<G> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<G: GraphLike> std::error::Error for ExtractError<G> {}

pub trait ToCircuit: GraphLike {
    fn into_circuit(&mut self) -> Result<Circuit, ExtractError<Self>>;
    fn to_circuit(&self) -> Result<Circuit, ExtractError<Self>> {
        self.clone().into_circuit()
    }

    fn extractor(&mut self) -> Extractor<Self> {
        Extractor::new(self)
    }
}

pub struct Extractor<'a, G: GraphLike> {
    g: &'a mut G,
    frontier: Vec<(usize,V)>,
    up_to_perm: bool,
    gaussf: fn(&mut Extractor<'a,G>, &mut Circuit),
}

impl<'a, G: GraphLike> Extractor<'a, G> {
    pub fn new(g: &'a mut G) -> Extractor<G> {
        Extractor {
            g,
            frontier: Vec::new(),
            up_to_perm: false,
            gaussf: Extractor::single_sln_set,
        }
    }

    pub fn with_gaussf(&mut self, f: fn(&mut Extractor<'a,G>, &mut Circuit)) -> &mut Self {
        self.gaussf = f;
        self
    }

    pub fn up_to_perm(&mut self) -> &mut Self {
        self.up_to_perm = true;
        self
    }

    pub fn flow(&mut self) -> &mut Self {
        self.with_gaussf(Extractor::no_gauss)
    }

    pub fn gflow(&mut self) -> &mut Self {
        self.with_gaussf(Extractor::single_sln_set)
    }

    pub fn gflow_simple_gauss(&mut self) -> &mut Self {
        self.with_gaussf(Extractor::simple_gauss)
    }

    /// Build a biadjacency matrix of frontier with its neighbors
    ///
    /// Frontier elements are rows and neighbors are columns. The computed
    /// vec of neighbors and the matrix are returned.
    fn frontier_biadj(&self) -> (Vec<V>, Mat2) {
        let mut neighbor_set = FxHashSet::default();
        for &(_,v) in &self.frontier {
            for n in self.g.neighbors(v) {
                if self.g.vertex_type(n) == VType::Z { neighbor_set.insert(n); }
            }
        }

        let neighbors: Vec<_> = neighbor_set.iter().copied().collect();

        // Build an adjacency matrix between the frontier and its neighbors
        let m = Mat2::build(self.frontier.len(), neighbors.len(), |i,j| {
            self.g.connected(self.frontier[i].1, neighbors[j])
        });

        (neighbors, m)
    }

    /// Set edges between frontier and given neighbors to match biadj. matrix
    fn update_frontier_biadj(&mut self, neighbors: &Vec<V>, m: Mat2) {
        for (i, &(_,v)) in self.frontier.iter().enumerate() {
            for (j, &w) in neighbors.iter().enumerate() {
                if m[(i,j)] == 1 {
                    if !self.g.connected(v, w) { self.g.add_edge_with_type(v, w, EType::H); }
                } else {
                    if self.g.connected(v, w) { self.g.remove_edge(v, w); }
                }
            }
        }
    }


    /// Push the gates in `c1` on to the front of `c`
    ///
    /// Note the order of gates will get reversed when doing this. Since c1 only refers to
    /// gates between frontier qubits, qubit indexes may need to be translated.
    fn update_frontier_circuit(&mut self, c1: &Circuit, c: &mut Circuit) {
        for gate in &c1.gates {
            // note the frontier might only be a subset of the qubits, so we should
            // lift to the global qubit index before adding a CNOT to the circuit
            let mut gate = gate.clone();
            gate.qs[0] = self.frontier[gate.qs[0]].0;
            gate.qs[1] = self.frontier[gate.qs[1]].0;
            c.push_front(gate);
        }
    }

    /// Don't do gaussian elimination on frontier
    ///
    /// This function is intended for extracting graphs that aleady have
    /// a causal flow.
    pub fn no_gauss(_: &mut Extractor<G>, _: &mut Circuit) {}

    /// Perform gaussian elimination on the frontier as CNOT gates
    ///
    /// At this point, we assume the frontier is phase-free and there are no edges
    /// between frontier vertices.
    pub fn simple_gauss(e: &mut Extractor<G>, c: &mut Circuit) {
        let (neighbors, mut m) = e.frontier_biadj();
        // Extract CNOTs until adj. matrix is in reduced echelon form. N.b. the
        // generated CNOTs correspond exactly to the row operations of the gaussian
        // elimination, but they are pushed on to the front of the circuit, so they
        // should end up in reverse order.
        let mut c1 = Circuit::new(c.num_qubits());
        m.gauss_x(true, 3, &mut c1);

        e.update_frontier_circuit(&c1, c);
        e.update_frontier_biadj(&neighbors, m);
    }

    /// Perform row operations to free a single vertex with the smallest solution set
    pub fn single_sln_set(e: &mut Extractor<G>, c: &mut Circuit) {
        let (neighbors, mut m) = e.frontier_biadj();
        let mut row_ops = Mat2::id(m.num_rows());
        let mut m1 = m.clone();
        m1.gauss_x(true, 1, &mut row_ops);
        let mut min_weight = row_ops.num_cols() as u8;
        let mut extr_rows = Vec::new();
        let mut min_weight_row = 0;

        // find the vertex with the smallest solution set
        for i in 0..m1.num_rows() {
            if m1.row_weight(i) == 1 {
                extr_rows.push(i);
                let weight = row_ops.row_weight(i);
                if weight <= min_weight {
                    min_weight_row = i;
                    min_weight = weight;
                }
            }
        }

        // compute the solution set
        let sln_set: Vec<_> = (0..row_ops.num_cols())
            .filter(|&i| row_ops[min_weight_row][i] == 1)
            .collect();

        if sln_set.len() < 2 { return; }


        let target = sln_set[0];

        // We can choose any qubit in the solution set as the target for the CNOTs. There are
        // lots of ways to choose this. Here, we choose the one that minimises the size of
        // the solution set of the next extractable vertex.

        // if extr_rows.len() > 1 {
        //     let cost_m = Mat2::build(extr_rows.len(), row_ops.num_cols(), |i,j| {
        //         row_ops[extr_rows[i]][j] == 1
        //     });
        //     let mut min_cost = (row_ops.num_cols() * sln_set.len()) as u8;

        //     for &t in &sln_set {
        //         let mut cm = cost_m.clone();
        //         for &i in &sln_set {
        //             cm.col_add(t, i);
        //         }

        //         let cost = (0..cm.num_rows()).map(|i| cm.row_weight(i)).min().unwrap();
        //         if cost < min_cost {
        //             target = t;
        //             min_cost = cost;
        //         }
        //     }
        // }


        let mut c1 = Circuit::new(c.num_qubits());
        for &i in &sln_set {
            if i != target {
                m.row_add(i, target);
                c1.row_add(i, target);
            }
        }

        // println!("Got {} extractable verts.", num_extr);
        // println!("New adj:\n{}", m);
        // let mut m2 = m.clone();
        // m2.gauss(true);
        // println!("Reduced:\n{}", m2);

        e.update_frontier_circuit(&c1, c);
        e.update_frontier_biadj(&neighbors, m);
    }

    /// Converts a permutation graph to a circuit of CNOT gates
    ///
    /// A permutation graph contains only inputs, outputs, and normal edges
    /// connecting inputs to outputs.
    fn perm_to_cnots(&mut self, c: &mut Circuit, blocksize: usize) {
        let mut m = Mat2::build(self.g.outputs().len(), self.g.inputs().len(), |i,j| {
            self.g.connected(self.g.outputs()[i], self.g.inputs()[j])
        });

        // Extract CNOTs until adj. matrix is in reduced echelon form
        let mut c1 = Circuit::new(c.num_qubits());
        m.gauss_x(true, blocksize, &mut c1);
        for g in c1.gates { c.push_front(g); }
    }

    /// Prepare the frontier for circuit extraction
    ///
    /// Identifies the frontier, and pulls Hadamards, phases, and CZ
    /// gates into the circuit. Returns the frontier as a Vec of pairs
    /// (qubit, frontier_vertex).
    fn prepare_frontier(&mut self, c: &mut Circuit) -> Result<(), ExtractError<G>> {
        self.frontier = Vec::new();

        for q in 0..self.g.outputs().len() {
            let o = self.g.outputs()[q];
            if let Some((v,et)) = self.g.incident_edges(o).next() {
                // replace a Hadamard edge from the output with a Hadamard gate
                if et == EType::H {
                    c.push_front(Gate::new(HAD, vec![q]));
                    self.g.set_edge_type(v, o, EType::N);
                }

                // output connects to an input, so skip. When the MAIN PHASE of
                // extraction is done, all vertices will be skipped this way.
                if self.g.vertex_type(v) == VType::B { continue; }

                self.frontier.push((q,v));

                // replace a non-zero phase on the frontier with a phase gate
                let p = self.g.phase(v);
                if !p.is_zero() {
                    c.push_front(Gate::new_with_phase(ZPhase, vec![q], p));
                    self.g.set_phase(v, Rational::zero());
                }

                // inspect neighbors of the frontier vertex
                for n in self.g.neighbor_vec(v) {
                    // ignore the output neighbour
                    if n == o {
                        continue;
                        // if another boundary is encountered...
                    } else if self.g.vertex_type(n) == VType::B {
                        // for unitary circuits, an additional boundary must be an input
                        if !self.g.inputs().contains(&n) {
                            return Err(ExtractError(format!("Two outputs connected to a single vertex {}.", v),
                            c.clone(), self.g.clone()));
                        }

                        // if a vertex is connected to more than just an input, pad the input
                        // with a dummy identity spider
                        if self.g.degree(v) > 2 {
                            let vd = VData {
                                ty: VType::Z,
                                phase: Rational::zero(),
                                qubit: self.g.qubit(n),
                                row: self.g.row(n)+1 };
                            let n1 = self.g.add_vertex_with_data(vd);
                            self.g.add_edge_with_type(n, n1, self.g.edge_type(n, v).opposite());
                            self.g.add_edge_with_type(n1, v, EType::H);
                            self.g.remove_edge(n, v);
                        }

                        // if the frontier vertex is connected to another frontier vertex, replace
                        // the edge with a CZ gate
                    } else if let Some(&(r,_)) = self.frontier.iter().find(|&&(_,n1)| n == n1) {
                        // TODO: CZ optimisation (maybe)
                        self.g.remove_edge(v, n);
                        c.push_front(Gate::new(CZ, vec![q,r]));

                        // we should not encounter any non-Z vertices at this point
                    } else if self.g.vertex_type(n) != VType::Z {
                        return Err(ExtractError(format!("Bad neighbour: {}", n), c.clone(), self.g.clone()));
                    }
                }
            } else {
                // this will happen if there is an output vertex not connected to anything, which
                // is a mal-formed graph
                return Err(ExtractError(format!("Bad output vertex {}", o), c.clone(), self.g.clone()));
            }
        }

        Ok(())
    }

    /// Pivot to remove gadgets adjacent to the frontier
    fn fix_gadgets(&mut self,
                   c: &Circuit,
                   gadgets: &mut FxHashSet<V>)
        -> Result<bool, ExtractError<G>>
        {

            for &(_,v) in &self.frontier {
                for n in self.g.neighbor_vec(v) {
                    if gadgets.contains(&n) {
                        // TODO: this can be probably be done with
                        // gen_pivot_unsafe
                        // let t = g.to_tensor4();
                        if boundary_pivot(self.g, v, n) {
                            // println!("FIXED GADGET: ({}, gad = {})", v, n);
                            // assert_eq!(t, g.to_tensor4());
                            // println!("{}", g.to_dot());
                            gadgets.remove(&n);
                            return Ok(true);
                        } else {
                            return Err(ExtractError(format!("Could not remove gadget by pivoting: ({}, {})", v, n),
                            c.clone(), self.g.clone()));
                        }
                    }
                }
            }
            Ok(false)
        }

    /// Extract vertices from the frontier
    ///
    /// Look for frontier elements that are phase-free and degree 2, and replace them
    /// with identity. Returns true if we got any.
    fn extract_from_frontier(&mut self) -> bool {
        let mut found = false;
        for &(_,v) in &self.frontier {
            if remove_id(self.g, v) {
                found = true;
                // println!("EXTRACTED: {}", v);
            }
        }
        found
    }

    pub fn extract(&mut self) -> Result<Circuit, ExtractError<G>> {
        // let t = self.to_tensor4(); // DEBUG
        let mut c = Circuit::new(self.g.outputs().len());

        // Pre-generate a set of all the phase gadgets. The extraction should
        // only ever eliminate phase gadgets, never create new ones.
        let mut gadgets = FxHashSet::default();
        for v in self.g.vertices() {
            if self.g.degree(v) == 1 && self.g.vertex_type(v) == VType::Z
            {
                let n = self.g.neighbors(v).next().unwrap();
                if self.g.vertex_type(n) == VType::Z {
                    gadgets.insert(n);
                }
            }
        }
        // println!("gadgets: {:?}", gadgets);

        loop {
            // PREPROCESSING PHASE
            //
            // Remove any phases, Hadamards, or CZs from the output and generate
            // a list of frontier vertices. If the frontier is empty after pre-processing,
            // we are done.
            self.prepare_frontier(&mut c)?;
            if self.frontier.is_empty() { break; }

            // Uncomment to debug {{{
            // println!("frontier: {:?}", frontier);
            // let t1 = self.g.to_tensor4().plug_n_qubits(c.num_qubits(), &c.to_tensor4());
            // assert!(Tensor4::scalar_eq(&t, &t1));
            // }}}

            // GADGET PHASE
            //
            // If any gadgets are adjacent to the frontier, do a generalised pivot to remove
            // them. In that case, some edges will change, so we need to re-generate the frontier.
            if self.fix_gadgets(&c, &mut gadgets)? { continue; }

            // MAIN PHASE
            //
            // Look for extractible vertices. If we found some, loop. If not, try gaussian
            // elimination via CNOTs and look again.
            if self.extract_from_frontier() { continue; }

            let gaussf = self.gaussf;
            gaussf(self, &mut c);

            if self.extract_from_frontier() { continue; }

            // If we didn't make progress, terminate with an error. This prevents infinite loops
            // in the case where a graph is not extractible.
            return Err(ExtractError("No extractible vertex found.".into(), c, self.g.clone()));
        }

        // FINAL PERMUTATION PHASE
        //
        // Generate CNOTs to turn the final permutation into the identity
        if !self.up_to_perm {
            self.perm_to_cnots(&mut c, 3);
        }

        Ok(c)
    }
}

impl<G: GraphLike + Clone> ToCircuit for G {
    fn into_circuit(&mut self) -> Result<Circuit, ExtractError<G>> {
        Extractor::new(self).extract()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec_graph::Graph;
    use crate::tensor::*;
    use crate::simplify::*;

    #[test]
    fn id_test() {
        let mut g = Graph::new();

        let is: Vec<_> = (0..4).map(|_| g.add_vertex(VType::B)).collect();
        let os: Vec<_> = (0..4).map(|_| g.add_vertex(VType::B)).collect();
        g.add_edge(is[0], os[0]);
        g.add_edge(is[1], os[1]);
        g.add_edge(is[2], os[2]);
        g.add_edge(is[3], os[3]);
        g.set_inputs(is);
        g.set_outputs(os);

        let mut c = Circuit::new(4);
        let mut g1 = g.clone();
        let mut e = Extractor::new(&mut g1);
        e.perm_to_cnots(&mut c, 3);
        // c.adjoint();
        println!("{}", c);
        // panic!("foo");

        assert_eq!(g.to_tensor4(), c.to_tensor4());
    }

    #[test]
    fn perm_test() {
        let mut g = Graph::new();

        let is: Vec<_> = (0..4).map(|_| g.add_vertex(VType::B)).collect();
        let os: Vec<_> = (0..4).map(|_| g.add_vertex(VType::B)).collect();
        g.add_edge(is[0], os[1]);
        g.add_edge(is[1], os[2]);
        g.add_edge(is[2], os[0]);
        g.add_edge(is[3], os[3]);
        g.set_inputs(is);
        g.set_outputs(os);

        let mut c = Circuit::new(4);
        let mut g1 = g.clone();
        let mut e = Extractor::new(&mut g1);
        e.perm_to_cnots(&mut c, 3);
        // c.adjoint();
        println!("{}", c);
        // panic!("foo");

        assert_eq!(g.to_tensor4(), c.to_tensor4());
    }

    #[test]
    fn extract_h() {
        let c = Circuit::from_qasm(r#"
            qreg q[1];
            h q[0];
        "#).unwrap();
        let mut g: Graph = c.to_graph();
        clifford_simp(&mut g);
        println!("GRAPH BEFORE: {}\n\n", g.to_dot());

        match g.to_circuit() {
            Ok(c1) => {
                println!("CIRCUIT: {}\n", c1);
                assert_eq!(c.to_tensor4(), c1.to_tensor4());
            },
            Err(ExtractError(msg, c1, g)) => {
                println!("CIRCUIT: {}\n\nGRAPH: {}\n", c1, g.to_dot());
                panic!("Extraction failed: {}", msg);
            }
        }
    }

    #[test]
    fn extract_swap() {
        let c = Circuit::from_qasm(r#"
            qreg q[2];
            cx q[0], q[1];
            cx q[1], q[0];
            cx q[0], q[1];
        "#).unwrap();
        let mut g: Graph = c.to_graph();
        clifford_simp(&mut g);
        println!("GRAPH BEFORE: {}\n\n", g.to_dot());

        match g.to_circuit() {
            Ok(c1) => { assert_eq!(c.to_tensor4(), c1.to_tensor4()); },
            Err(ExtractError(msg, c1, g)) => {
                println!("CIRCUIT: {}\n\nGRAPH: {}\n", c1, g.to_dot());
                panic!("Extraction failed: {}", msg);
            }
        }
    }

    #[test]
    fn extract1() {
        let c = Circuit::from_qasm(r#"
            qreg q[4];
            cx q[0], q[1];
            cx q[0], q[2];
            cx q[0], q[3];
            cx q[1], q[2];
            cx q[2], q[1];
            cx q[1], q[2];
            cx q[1], q[3];
            cx q[1], q[0];
        "#).unwrap();
        let mut g: Graph = c.to_graph();
        clifford_simp(&mut g);
        println!("GRAPH BEFORE: {}\n\n", g.to_dot());

        match g.to_circuit() {
            Ok(c1) => {
                println!("CIRCUIT: {}\n", c1);
                assert_eq!(c.to_tensor4(), c1.to_tensor4());
            },
            Err(ExtractError(msg, c1, g)) => {
                println!("CIRCUIT: {}\n\nGRAPH: {}\n", c1, g.to_dot());
                panic!("Extraction failed: {}", msg);
            }
        }
    }

    #[test]
    fn random_flow_extract() {
        // this particular circuit never calls gauss_frontier
        let c = Circuit::random()
            .seed(1337)
            .qubits(5)
            .depth(20)
            .p_t(0.2)
            .with_cliffords()
            .build();
        let mut g: Graph = c.to_graph();
        clifford_simp(&mut g);

        assert_eq!(c.to_tensor4(), g.to_tensor4());
        let c1 = g.to_circuit().expect("Circuit should extract.");
        assert!(Tensor4::scalar_compare(&c, &c1));
    }

    #[test]
    fn random_gflow_extract() {
        // this particular circuit does call gauss_frontier
        let c = Circuit::random()
            .seed(1337)
            .qubits(5)
            .depth(30)
            .p_t(0.2)
            .with_cliffords()
            .build();

        let mut g: Graph = c.to_graph();
        clifford_simp(&mut g);

        assert_eq!(c.to_tensor4(), g.to_tensor4());
        let c1 = g.to_circuit().expect("Circuit should extract.");
        assert!(Tensor4::scalar_compare(&c, &c1));
    }

    #[test]
    fn random_extract() {
        let c = Circuit::random()
            .seed(1337)
            .qubits(10)
            .depth(40)
            .p_t(0.2)
            .with_cliffords()
            .build();
        let mut g: Graph = c.to_graph();
        clifford_simp(&mut g);
        println!("{}", g.to_dot());
        let _c1 = g.to_circuit().expect("Circuit should extract.");
    }

    #[test]
    fn regression_extract_1() {
        // caused bug toward the end of extraction, when frontier was only
        // a subset of qubits
        let c = Circuit::from_qasm(r#"
          qreg q[5];
          cx q[3], q[4];
          tdg q[4];
          cx q[0], q[3];
          tdg q[3];
          cx q[0], q[3];
          cx q[1], q[4];
          cx q[0], q[4];
          cx q[1], q[4];
          tdg q[4];
          t q[0];
          "#).unwrap();

          let mut g: Graph = c.to_graph();
          clifford_simp(&mut g);
          assert!(Tensor4::scalar_compare(&g, &c));
          let c1 = g.to_circuit().unwrap();
          assert!(Tensor4::scalar_compare(&c, &c1));
    }
}
