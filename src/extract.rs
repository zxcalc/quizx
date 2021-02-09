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
                    if self.inputs().contains(&v) { continue; }
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
            let neighbors: Vec<_> = neighbor_set.iter().copied().collect();
            let m = Mat2::build(frontier.len(), neighbors.len(), |i,j| {
                self.connected(frontier[i], neighbors[j])
            });

            break; // TODO: finish!
        }

        Ok(c)
    }
}
