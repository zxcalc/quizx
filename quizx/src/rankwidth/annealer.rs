use rand::Rng;

use crate::graph::GraphLike;
use crate::rankwidth::decomp_tree::DecompTree;

pub struct RankwidthAnnealer<R: Rng, G: GraphLike> {
    rng: R,
    graph: G,
    init_decomp: DecompTree,
    init_temp: f64,
    min_temp: f64,
    cooling_rate: f64,
    adaptive_cooling: bool,
    iterations: usize,
}

impl<R: Rng, G: GraphLike> RankwidthAnnealer<R, G> {
    pub fn new(graph: G, mut rng: R) -> Self {
        let init_decomp = DecompTree::random_decomp(&graph, &mut rng);
        Self {
            rng,
            graph,
            init_decomp,
            init_temp: 5.0,
            min_temp: 0.01,
            cooling_rate: 0.95,
            adaptive_cooling: true,
            iterations: 1000,
        }
    }

    pub fn new_with_decomp(graph: G, init_decomp: DecompTree, rng: R) -> Self {
        Self {
            rng,
            graph,
            init_decomp,
            init_temp: 5.0,
            min_temp: 0.05,
            cooling_rate: 0.95,
            adaptive_cooling: true,
            iterations: 1000,
        }
    }

    pub fn set_init_decomp(&mut self, init_decomp: DecompTree) -> &mut Self {
        self.init_decomp = init_decomp;
        self
    }

    pub fn set_init_temp(&mut self, init_temp: f64) -> &mut Self {
        self.init_temp = init_temp;
        self
    }

    pub fn set_min_temp(&mut self, min_temp: f64) -> &mut Self {
        self.min_temp = min_temp;
        self
    }

    pub fn set_cooling_rate(&mut self, cooling_rate: f64) -> &mut Self {
        self.cooling_rate = cooling_rate;
        self
    }

    pub fn set_adaptive_cooling(&mut self, adaptive_cooling: bool) -> &mut Self {
        self.adaptive_cooling = adaptive_cooling;
        self
    }

    pub fn set_iterations(&mut self, iterations: usize) -> &mut Self {
        self.iterations = iterations;
        self
    }

    pub fn init_decomp(&self) -> &DecompTree {
        &self.init_decomp
    }

    pub fn init_temp(&self) -> f64 {
        self.init_temp
    }

    pub fn min_temp(&self) -> f64 {
        self.min_temp
    }

    pub fn cooling_rate(&self) -> f64 {
        self.cooling_rate
    }

    pub fn adaptive_cooling(&self) -> bool {
        self.adaptive_cooling
    }

    pub fn iterations(&self) -> usize {
        self.iterations
    }

    pub fn run(&mut self) -> DecompTree {
        let mut best_width = self.init_decomp.rankwidth(&self.graph);
        let mut best_score = self.init_decomp.rankwidth_score(&self.graph);
        let mut old_score = self.init_decomp.rankwidth_score(&self.graph);
        let mut best_decomp = self.init_decomp.clone();
        let mut old_decomp = self.init_decomp.clone();
        let mut temp = self.init_temp;

        // operators chosen with weights:
        //  - leaf swap: 1
        //  - local swap: 4
        //  - subtree move: 5
        for _ in 0..self.iterations {
            let mut decomp = old_decomp.clone();
            let op = self.rng.random_range(0..10);
            if op < 1 {
                decomp.swap_random_leaves(&mut self.rng);
            } else if op < 5 {
                decomp.random_local_swap(&mut self.rng);
            } else {
                decomp.move_random_subtree(&mut self.rng);
            }

            let score = decomp.rankwidth_score(&self.graph);

            let keep = if score < old_score {
                true
            } else {
                let delta = old_score as f64 - score as f64;
                let t = if self.adaptive_cooling {
                    temp * (1.0 + (score as f64 - best_score as f64) / best_score as f64)
                } else {
                    temp
                };
                let prob = (delta / t).exp();
                if prob > 1.0 {
                    true
                } else {
                    self.rng.random_bool(prob)
                }
            };

            if keep {
                let width = decomp.rankwidth(&self.graph);
                if width < best_width {
                    best_width = width;
                    best_decomp = decomp.clone();
                }

                if score < best_score {
                    best_score = score;
                }

                old_score = score;
                old_decomp = decomp;
            }

            temp *= self.cooling_rate;

            if temp < self.min_temp {
                break;
            }
        }
        best_decomp
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;

    use super::*;
    use crate::graph::{EType, VType};
    use crate::vec_graph::Graph;

    #[test]
    fn test_rankwidth_annealer() {
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);
        let mut graph = Graph::new();
        let num_v = 40;
        let num_e = 400;
        for _ in 0..num_v {
            graph.add_vertex(VType::Z);
        }

        let mut e = 0;
        while e < num_e {
            let v1 = rng.random_range(0..num_v);
            let v2 = rng.random_range(0..num_v);
            if v1 != v2 && !graph.connected(v1, v2) {
                graph.add_edge_with_type(v1, v2, EType::H);
                e += 1;
            }
        }

        println!("Finished constructing graph");

        let mut annealer = RankwidthAnnealer::new(graph.clone(), rng);
        assert!(
            annealer.init_decomp.is_valid_for_graph(&graph),
            "Initial decomposition tree is not valid"
        );
        annealer.set_iterations(1000);
        println!(
            "Initial rankwidth: {}",
            annealer.init_decomp.rankwidth(&graph)
        );
        println!(
            "Initial score: {}",
            annealer.init_decomp.rankwidth_score(&graph)
        );
        let mut decomp = annealer.run();
        assert!(
            decomp.is_valid_for_graph(&graph),
            "Final decomposition tree is not valid"
        );
        println!("Final rankwidth: {}", decomp.rankwidth(&graph));
        println!("Final score: {}", decomp.rankwidth_score(&graph));
    }
}
