use crate::graph::*;
use num::Rational64;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

pub struct EquatorialStabilizerStateBuilder {
    pub rng: StdRng,
    pub qubits: usize,
}

impl Default for EquatorialStabilizerStateBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl EquatorialStabilizerStateBuilder {
    pub fn new() -> EquatorialStabilizerStateBuilder {
        EquatorialStabilizerStateBuilder {
            rng: StdRng::from_os_rng(),
            qubits: 1,
        }
    }

    pub fn seed(&mut self, seed: u64) -> &mut Self {
        self.rng = StdRng::seed_from_u64(seed);
        self
    }
    pub fn qubits(&mut self, qubits: usize) -> &mut Self {
        self.qubits = qubits;
        self
    }
    pub fn build<G: GraphLike>(&mut self) -> G {
        let mut g = G::new();
        let outputs: Vec<_> = (0..self.qubits).map(|_| g.add_vertex(VType::B)).collect();
        let spiders: Vec<_> = (0..self.qubits).map(|_| g.add_vertex(VType::Z)).collect();

        let mut num_cz = 0;
        for i in 0..self.qubits {
            g.add_edge(spiders[i], outputs[i]);
            g.set_phase(spiders[i], Rational64::new(self.rng.random_range(0..3), 2));

            for j in 0..i {
                if self.rng.random_bool(0.5) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec_graph::Graph;

    #[test]
    fn default_builder() {
        let builder = EquatorialStabilizerStateBuilder::default();
        assert_eq!(builder.qubits, 1);
    }

    #[test]
    fn seeded_deterministic() {
        let g1: Graph = EquatorialStabilizerStateBuilder::new()
            .seed(42)
            .qubits(5)
            .build();
        let g2: Graph = EquatorialStabilizerStateBuilder::new()
            .seed(42)
            .qubits(5)
            .build();

        assert_eq!(g1.num_vertices(), g2.num_vertices());
        assert_eq!(g1.num_edges(), g2.num_edges());

        // Same phases on corresponding vertices
        for v in g1.vertices() {
            assert_eq!(g1.phase(v), g2.phase(v));
        }
    }

    #[test]
    fn graph_structure() {
        let g: Graph = EquatorialStabilizerStateBuilder::new()
            .seed(123)
            .qubits(4)
            .build();

        assert_eq!(g.outputs().len(), 4);
        assert!(g.inputs().is_empty());

        // Each output is a boundary vertex connected to exactly one Z spider
        for &out in g.outputs() {
            assert_eq!(g.vertex_type(out), VType::B);
            assert_eq!(g.degree(out), 1);
            let neighbor = g.neighbors(out).next().unwrap();
            assert_eq!(g.vertex_type(neighbor), VType::Z);
        }
    }

    #[test]
    fn phases_are_clifford() {
        // Equatorial stabilizer states have phases in {0, 1/2, 1}.
        use crate::phase::Phase;
        use num::One;

        let g: Graph = EquatorialStabilizerStateBuilder::new()
            .seed(999)
            .qubits(10)
            .build();

        for v in g.vertices() {
            if g.vertex_type(v) == VType::Z {
                let phase = g.phase(v);
                let valid = phase == Phase::from((0, 1))
                    || phase == Phase::from((1, 2))
                    || phase == Phase::one();
                assert!(valid, "Phase {} is not in {{0, 1/2, 1}}", phase);
            }
        }
    }
}
