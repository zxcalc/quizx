use crate::circuit::*;
use crate::gate::*;
use crate::graph::*;
use crate::linalg::*;
use num::{Rational, Zero};
use rustc_hash::FxHashSet;
use crate::basic_rules::{gen_pivot, remove_id};

/// Extraction couldn't finish. Returns a message, a
/// partially-extracted circuit, and the remainder of
/// the graph.
pub type ExtractError<G> = (String, Circuit, G);

trait ToCircuit: Clone {
    fn into_circuit(self) -> Result<Circuit, ExtractError<Self>>;
    fn to_circuit(&self) -> Result<Circuit, ExtractError<Self>> {
        self.clone().into_circuit()
    }
}

/// Converts a permutation graph to a circuit of CNOT gates
///
/// A permutation graph contains only inputs, outputs, and normal edges
/// connecting inputs to outputs.
fn perm_to_cnots(g: &impl GraphLike, c: &mut Circuit, blocksize: usize) {
    let mut m = Mat2::build(g.outputs().len(), g.inputs().len(), |i,j| {
        g.connected(g.outputs()[i], g.inputs()[j])
    });

    // Extract CNOTs until adj. matrix is in reduced echelon form
    m.gauss_y(true, blocksize, c);
}

/// Prepare the frontier for circuit extraction
///
/// Identifies the frontier, and pulls Hadamards, phases, and CZ
/// gates into the circuit. Returns the frontier as a Vec of pairs
/// (qubit, frontier_vertex).
fn prepare_frontier<G: GraphLike>(g: &mut G, c: &mut Circuit) -> Result<Vec<(usize,V)>, ExtractError<G>> {
    let mut frontier = Vec::new();
    let outputs = g.outputs().clone();

    for (q,&o) in outputs.iter().enumerate() {
        if let Some((v,et)) = g.incident_edges(o).next() {
            if et == EType::H {
                c.push_front(Gate::new(HAD, vec![q]));
                g.set_edge_type(v, o, EType::N);
            }

            // output connects to an input, so skip
            if g.vertex_type(v) == VType::B { continue; }

            frontier.push((q,v));

            let p = g.phase(v);
            if !p.is_zero() {
                c.push_front(Gate::new_with_phase(ZPhase, vec![q], p));
                g.set_phase(v, Rational::zero());
            }

            for n in g.neighbor_vec(v) {
                if n == o {
                    continue;
                } else if g.vertex_type(n) == VType::B {
                    if g.degree(v) == 2 {
                        frontier.pop();

                        if g.edge_type(v, n) == EType::H {
                            c.push_front(Gate::new(HAD, vec![q]));
                        }

                        g.remove_vertex(v);
                        g.add_edge(o, n);

                        break;
                    } else {
                        let vd = VData {
                            ty: VType::Z,
                            phase: Rational::zero(),
                            qubit: g.qubit(n),
                            row: g.row(n)+1 };
                        let n1 = g.add_vertex_with_data(vd);
                        g.add_edge_with_type(n, n1,
                                                g.edge_type(n, v).opposite());
                        g.add_edge_with_type(n1, v, EType::H);
                        g.remove_edge(n, v);
                    }
                } else if let Some(&(r,_)) = frontier.iter().find(|&&(_,n1)| n == n1) {
                    // TODO: CZ optimisation (maybe)
                    g.remove_edge(v, n);
                    c.push_front(Gate::new(CZ, vec![q,r]));
                } else if g.vertex_type(n) != VType::Z {
                    return Err((format!("Bad neighbour: {}", n), c.clone(), g.clone()));
                }
            }
        } else {
            return Err((format!("Bad output vertex {}", o), c.clone(), g.clone()));
        }
    }

    Ok(frontier)
}


/// Extract vertices from the frontier
///
/// Look for frontier elements that are phase-free and degree 2, and replace them
/// with identity. Returns true if we got any.
fn extract_from_frontier<G: GraphLike>(g: &mut G, frontier: &Vec<(usize,V)>) -> bool {
    let mut found = false;
    for &(_,v) in frontier { found = found || remove_id(g, v); }
    found
}

/// Perform gaussian elimination on the frontier as CNOT gates
///
/// At this point, we assume the frontier is phase-free and there are no edges
/// between frontier vertices.
fn gauss_frontier<G: GraphLike>(g: &mut G, c: &mut Circuit, frontier: &Vec<(usize,V)>) {
    let mut neighbor_set = FxHashSet::default();
    for &(_,v) in frontier {
        for n in g.neighbors(v) {
            if g.vertex_type(v) == VType::Z { neighbor_set.insert(n); }
        }
    }

    let neighbors: Vec<_> = neighbor_set.iter().copied().collect();

    // Build an adjacency matrix between the frontier and its neighbors
    let mut m = Mat2::build(frontier.len(), neighbors.len(), |i,j| {
        g.connected(frontier[i].1, neighbors[j])
    });

    // Extract CNOTs until adj. matrix is in reduced echelon form
    m.gauss_y(true, 3, c);

    for (i, &(_,v)) in frontier.iter().enumerate() {
        for (j, &w) in neighbors.iter().enumerate() {
            if m[(i,j)] == 1 {
                if !g.connected(v, w) { g.add_edge_with_type(v, w, EType::H); }
            } else {
                if g.connected(v, w) { g.remove_edge(v, w); }
            }
        }
    }
}

impl<G: GraphLike + Clone> ToCircuit for G {
    fn into_circuit(mut self) -> Result<Circuit, ExtractError<G>> {
        // use GType::*;
        let mut c = Circuit::new(self.outputs().len());
        let mut gadgets = FxHashSet::default();

        for v in self.vertices() {
            if self.degree(v) == 1 &&
               self.vertex_type(v) == VType::Z
            {
                let n = self.neighbors(v).next().unwrap();
                if self.vertex_type(n) == VType::Z {
                    gadgets.insert(n);
                }
            }
        }

        'outer: loop {
            // PREPROCESSING PHASE
            let frontier = prepare_frontier(&mut self, &mut c)?;

            // if the frontier is empty after pre-processing, we are done
            if frontier.is_empty() { break; }

            //
            // GADGET PHASE
            //
            // If there are any phase gadgets adjacent to the frontier, we can
            // convert them to spiders with a pivot. This can change some
            // adjacency with the frontier, so we should restart the main
            // loop in this case.
            for &(_,v) in &frontier {
                for n in self.neighbor_vec(v) {
                    if gadgets.contains(&n) {
                        // TODO: this can be probably be done with
                        // gen_pivot_unsafe
                        if gen_pivot(&mut self, v, n) {
                            gadgets.remove(&n);
                            continue 'outer;
                        } else {
                            return Err((format!("Could not remove gadget by pivoting: ({}, {})", v, n),
                                        c, self));
                        }
                    }
                }
            }

            //
            // MAIN PHASE
            //
            if extract_from_frontier(&mut self, &frontier) { continue; }

            // gauss_frontier(&mut self, &mut c, &frontier);

            // if extract_from_frontier(&mut self, &frontier) { continue; }

            // If we didn't make progress, terminate with an error.
            return Err(("No extractible vertex found.".into(), c, self));
        }

        // FINAL PERMUATION PHASE
        //
        // Generate CNOTs to turn the final permutation into the identity
        perm_to_cnots(&self, &mut c, 3);

        Ok(c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec_graph::Graph;
    use crate::tensor::ToTensor;
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
        perm_to_cnots(&mut g.clone(), &mut c, 3);
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
        perm_to_cnots(&mut g.clone(), &mut c, 3);
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
            Err((msg, c1, g)) => {
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
            Err((msg, c1, g)) => {
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
            Err((msg, c1, g)) => {
                println!("CIRCUIT: {}\n\nGRAPH: {}\n", c1, g.to_dot());
                panic!("Extraction failed: {}", msg);
            }
        }
    }
}
