use crate::graph::*;
use rand::{SeedableRng, Rng};
use rand::rngs::StdRng;
use num::Rational;

pub struct EquatorialStabilizerStateBuilder {
    pub rng: StdRng,
    pub qubits: usize,
}

impl EquatorialStabilizerStateBuilder {
    pub fn new() -> EquatorialStabilizerStateBuilder {
        EquatorialStabilizerStateBuilder {
            rng: StdRng::from_entropy(),
            qubits: 1,
        }
    }

    pub fn seed(&mut self, seed: u64) -> &mut Self { self.rng = StdRng::seed_from_u64(seed); self }
    pub fn qubits(&mut self, qubits: usize) -> &mut Self { self.qubits = qubits; self }
    pub fn build<G: GraphLike>(&mut self) -> G {
        let mut g = G::new();
        let outputs: Vec<_> = (0..self.qubits).map(|_| g.add_vertex(VType::B)).collect();
        let spiders: Vec<_> = (0..self.qubits).map(|_| g.add_vertex(VType::Z)).collect();

        let mut num_cz = 0;
        for i in 0..self.qubits {
            g.add_edge(spiders[i], outputs[i]);
            g.set_phase(spiders[i], Rational::new(self.rng.gen_range(0..3), 2));

            for j in 0..i {
                if self.rng.gen_bool(0.5) {
                    g.add_edge_with_type(spiders[i], spiders[j], EType::H);
                    num_cz += 1;
                }
            }
        }

        g.set_outputs(outputs);
        g.scalar_mut().mul_sqrt2_pow(num_cz - (self.qubits as i32));

        g
    }
}


