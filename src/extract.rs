use crate::circuit::*;
use crate::gate::*;
use crate::graph::*;
use crate::linalg::*;
use num::{Rational, Zero};
use rustc_hash::{FxHashMap,FxHashSet};
use crate::basic_rules::gen_pivot;

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

fn perm_to_cnots(g: &impl GraphLike, c: &mut Circuit, blocksize: usize) {
    let mut m = Mat2::build(g.inputs().len(), g.outputs().len(), |i,j| {
        g.connected(g.inputs()[i], g.outputs()[j])
    });

    // Extract CNOTs until adj. matrix is in reduced echelon form
    m.gauss_aux(true, blocksize, c);
}

impl<G: GraphLike + Clone> ToCircuit for G {
    fn into_circuit(mut self) -> Result<Circuit, ExtractError<G>> {
        use GType::*;
        let mut c = Circuit::new(self.outputs().len());
        let mut gadgets = FxHashSet::default();
        let mut qubit_map = FxHashMap::default();

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

        let outputs = self.outputs().clone();

        'outer: loop {
            let mut frontier = Vec::new();
            let mut neighbor_set = FxHashSet::default();

            //
            // PREPROCESSING PHASE
            //
            // build the frontier and extract any phases on the frontier
            // as phase gates and any edges between frontier nodes as
            // CZ gates.
            for (q,&o) in outputs.iter().enumerate() {
                if let Some(v) = self.neighbors(o).next() {
                    // output connects to an input, so skip
                    if self.vertex_type(v) == VType::B { continue; }

                    frontier.push(v);
                    qubit_map.insert(v, q);

                    let et = self.edge_type(v,o);
                    if et == EType::H {
                        c.push(Gate::new(HAD, vec![q]));
                        self.set_edge_type(v, o, EType::N);
                    }

                    let p = self.phase(v);
                    if !p.is_zero() {
                        c.push(Gate::new_with_phase(ZPhase, vec![q], p));
                        self.set_phase(v, Rational::zero());
                    }

                    for n in self.neighbor_vec(v) {
                        if n == o {
                            continue;
                        } else if self.vertex_type(n) == VType::B {
                            if self.degree(v) == 2 {
                                frontier.pop();

                                if self.edge_type(v, n) == EType::H {
                                    c.push(Gate::new(HAD, vec![q]));
                                }

                                self.remove_vertex(v);
                                self.add_edge(o, n);

                                break;
                            } else {
                                let vd = VData {
                                    ty: VType::Z,
                                    phase: Rational::zero(),
                                    qubit: self.qubit(n),
                                    row: self.row(n)+1 };
                                let n1 = self.add_vertex_with_data(vd);
                                self.add_edge_with_type(n, n1,
                                                        self.edge_type(n, v).opposite());
                                self.add_edge_with_type(n1, v, EType::H);
                                self.remove_edge(n, v);
                            }
                        } else if frontier.contains(&n) {
                            // TODO: CZ optimisation (maybe)
                            self.remove_edge(v, n);
                            c.push(Gate::new(CZ, vec![v,n]));
                        } else if self.vertex_type(n) == VType::Z {
                            neighbor_set.insert(n);
                        } else {
                            return Err((format!("Bad neighbour: {}", n), c, self));
                        }
                    }
                } else {
                    return Err((format!("Bad output vertex {}", o), c, self));
                }
            }

            // if the frontier is empty after pre-processing, we are done
            if frontier.is_empty() { break; }

            //
            // GADGET PHASE
            //
            // If there are any phase gadgets adjacent to the frontier, we can
            // convert them to spiders with a pivot. This can change some
            // adjacency with the frontier, so we should restart the main
            // loop in this case.
            for &v in &frontier {
                for n in self.neighbor_vec(v) {
                    if gadgets.contains(&n) {
                        // TODO: this can be probably be done with
                        // gen_pivot_unsafe
                        if gen_pivot(&mut self, v, n) {
                            gadgets.remove(&n);
                            continue 'outer;
                        } else {
                            return Err((format!("Could not remove gadget by pivoting: ({}, {})", v, n), c, self));
                        }
                    }
                }
            }

            //
            // MAIN PHASE
            //
            // Look for extractable spiders and extract them.

            // Build an adjacency matrix between the frontier and its neighbors
            let neighbors: Vec<_> = neighbor_set.iter().copied().collect();
            let mut m = Mat2::build(frontier.len(), neighbors.len(), |i,j| {
                self.connected(frontier[i], neighbors[j])
            });

            // Extract CNOTs until adj. matrix is in reduced echelon form
            m.gauss_aux(true, 3, &mut c);

            // Update graph with new adjacency, and look for elements of
            // the frontier with a unique neighbour. When this occurs,
            // save the frontier element and the neighbor in extract_pairs.
            //
            // n.b. add_edge_with_type should be a noop if an edge is already
            // there, and similarly remove_edge should be a noop if there is
            // no edge there.
            let mut extract_pairs = vec![];
            for (i, &v) in frontier.iter().enumerate() {
                let mut unique = true;
                let mut k = None;
                for (j, &w) in neighbors.iter().enumerate() {
                    if m[(i,j)] == 1 {
                        self.add_edge_with_type(v, w, EType::H);
                        if k.is_some() { unique = false; }
                        else { k = Some(j); }
                    } else {
                        self.remove_edge(v, w);
                    }
                }

                if unique {
                    k.map(|k1| extract_pairs.push((i,k1)));
                }
            }

            // If we can't make progress, return an error. This should prevent
            // infinite loops.
            if extract_pairs.is_empty() {
                return Err(("No extractible vertex found.".into(), c, self));
            }

            // Each extractible vertex takes the place of its neighbor in the
            // frontier. Since we have already removed the phase and the
            // hadamard edge from the frontier vertex, we can simply delete it
            // and attach its neighbor to the corresponding output.
            for (i,j) in extract_pairs {
                self.add_edge_with_type(neighbors[j], outputs[i], EType::H);
                self.remove_vertex(frontier[i]);
            }
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
}
