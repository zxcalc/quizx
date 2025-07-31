use rand::Rng;

use crate::{graph::GraphLike, rankwidth::RankWidthDecomposer};

pub struct RankwidthAnnealer<'a, R: Rng, G: GraphLike> {
    rng: R,
    init_decomp: RankWidthDecomposer<'a, G>,
    init_temp: f64,
    min_temp: f64,
    cooling_rate: f64,
    adaptive_cooling: bool,
    iterations: usize,
}

impl<'a, R: Rng, G: GraphLike> RankwidthAnnealer<'a, R, G> {
    pub fn from_graph(graph: &'a G, mut rng: R) -> Self {
        let init_decomp = RankWidthDecomposer::random_decomp(graph, &mut rng);
        Self {
            rng,
            init_decomp,
            init_temp: 5.0,
            min_temp: 0.05,
            cooling_rate: 0.95,
            adaptive_cooling: true,
            iterations: 1000,
        }
    }

    pub fn from_decomp(init_decomp: RankWidthDecomposer<'a, G>, rng: R) -> Self {
        Self {
            rng,
            init_decomp,
            init_temp: 5.0,
            min_temp: 0.05,
            cooling_rate: 0.95,
            adaptive_cooling: true,
            iterations: 1000,
        }
    }

    pub fn set_init_temp(&mut self, init_temp: f64) {
        self.init_temp = init_temp;
    }

    pub fn set_min_temp(&mut self, min_temp: f64) {
        self.min_temp = min_temp;
    }

    pub fn set_cooling_rate(&mut self, cooling_rate: f64) {
        self.cooling_rate = cooling_rate;
    }

    pub fn set_adaptive_cooling(&mut self, adaptive_cooling: bool) {
        self.adaptive_cooling = adaptive_cooling;
    }

    pub fn set_iterations(&mut self, iterations: usize) {
        self.iterations = iterations;
    }

    pub fn run(&mut self) -> RankWidthDecomposer<'a, G> {
        let mut best_width = self.init_decomp.rankwidth();
        let mut best_score = self.init_decomp.rankwidth_score();
        let mut old_score = self.init_decomp.rankwidth_score();
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

            let score = decomp.rankwidth_score();

            let keep = if score < old_score {
                true
            } else {
                let delta = old_score as f64 - score as f64;
                let t = if self.adaptive_cooling {
                    temp * (1.0 + (score as f64 - best_score as f64) / best_score as f64)
                } else {
                    temp
                };
                let prob = (-delta / t).exp();
                self.rng.random_bool(prob)
            };

            if keep {
                let width = decomp.rankwidth();
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
