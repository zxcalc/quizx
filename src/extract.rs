use crate::circuit::*;
use crate::graph::*;
use num::{Rational, Zero};
use rustc_hash::FxHashMap;

trait ToCircuit: Clone {
    fn into_circuit(self) -> Result<Circuit, String>;
    fn to_circuit(&self) -> Result<Circuit, String> {
        self.clone().into_circuit()
    }
}

impl<G: GraphLike + Clone> ToCircuit for G {
    fn into_circuit(mut self) -> Result<Circuit, String> {
        let mut c = Circuit::new(self.outputs().len());
        let mut gadgets = FxHashMap::default();
        let mut qubit_map = FxHashMap::default();

        for v in self.vertices() {
            if self.degree(v) == 1 &&
               !self.inputs().contains(&v) &&
               !self.outputs().contains(&v)
            {
                let n = self.neighbors(v).next().unwrap();
                gadgets.insert(n, v);
            }
        }

        let mut frontier = Vec::new();
        for (i,&o) in self.outputs().iter().enumerate() {
            let v = self.neighbors(o).next()
                .map_or(Err("Bad output vertex"), Ok)?;
            if self.inputs().contains(&v) { continue; }
            frontier.push(v);
            qubit_map.insert(v, i);
        }

        loop {
            for &v in &frontier {
                let q = qubit_map[&v];
                let b = self.neighbors(v)
                    .filter(|w| self.outputs().contains(w))
                    .next()
                    .expect("Frontier should be next to output.");
                let et = self.edge_type(v,b);
                if et == EType::H {
                    c.add_gate("had", vec![q]);
                    self.set_edge_type(v, b, EType::N);
                }

                let p = self.phase(v);
                if !p.is_zero() {
                    c.add_gate_with_phase("rz", vec![q], p);
                    self.set_phase(v, Rational::zero());
                }
            }
            break; // TODO: finish!
        }

        Ok(c)
    }
}
