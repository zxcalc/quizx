// QuiZX - Rust library for quantum circuit rewriting and optimization
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

//! Causal-flow preserving superoptimizer.
//!
//! Based on the paper "Exhaustive Generation of Quantum Circuit Rewrite Rules using ZX" [arxiv][]
//!
//! The code in this module is adapted from the [Badger superoptimizer] from `tket2`.
//!
//! TODO: arxiv link
//! [arxiv]: TODO
//! [Badger superoptimizer]: https://github.com/CQCL/tket2

use crate::cost::{CostMetric, TwoQubitGateCount};
use crate::hash::weisfeiler_lehman_graph_hash;
use crate::rewriter::Rewriter;

use std::hash::{Hash, Hasher};
use std::time::{Duration, Instant};

use priority_queue::DoublePriorityQueue;
use quizx::graph::GraphLike;
use rustc_hash::FxHashSet;

/// Configuration options for the superoptimizer.
#[derive(Copy, Clone, Debug)]
pub struct SuperOptOptions {
    /// The maximum time (in seconds) to run the optimizer.
    ///
    /// Defaults to `None`, which means no timeout.
    pub timeout: Option<u64>,
    /// The maximum time (in seconds) to search for new improvements to the
    /// graph. If no progress is made in this time, the optimizer will stop.
    ///
    /// Defaults to `None`, which means no timeout.
    pub progress_timeout: Option<u64>,
    /// The maximum size of the graph candidates priority queue.
    ///
    /// Defaults to `200`.
    pub queue_size: usize,
}

impl Default for SuperOptOptions {
    fn default() -> Self {
        Self {
            timeout: None,
            progress_timeout: None,
            queue_size: 200,
        }
    }
}

/// Causal flow preserving superoptimizer for ZX diagrams.
///
/// Based on the paper "Exhaustive Generation of Quantum Circuit Rewrite Rules using ZX" [arxiv][].
/// Adapted from the [Badger superoptimizer] from `tket2`.
///
/// Using a rewriter and a rewrite strategy, the optimizer
/// will repeatedly rewrite the graph, optimizing the graph according to
/// the cost function provided.
///
/// Optimization is done by maintaining a priority queue of graphs and
/// always processing the graph with the lowest cost first. Rewrites are
/// computed for that graph and all new graph obtained are added to the queue.
///
/// TODO: arxiv link
/// [arxiv]: TODO
/// [Badger superoptimizer]: https://github.com/CQCL/tket2
#[derive(Clone, Debug)]
pub struct SuperOptimizer<R, C = TwoQubitGateCount> {
    rewriter: R,
    cost_metric: C,
}

impl<R, C> SuperOptimizer<R, C>
where
    C: CostMetric,
    R: Rewriter,
{
    /// Create a new superoptimizer with the given rewriter and strategy.
    pub fn new(rewriter: R) -> Self {
        Self {
            rewriter,
            cost_metric: C::new(),
        }
    }

    pub fn optimize<G: GraphLike>(&self, g: &G, opt: SuperOptOptions) -> G {
        let start_time = Instant::now();
        let mut last_best_time = Instant::now();

        let mut best_g = g.clone();
        let mut best_g_cost = self.cost_metric.cost(g);
        log::info!("Initial cost: {:?}", best_g_cost);

        // Hash of seen graphs. Dot not store graphs as this map gets huge
        let hash = weisfeiler_lehman_graph_hash(g, 3);
        let mut seen_hashes = FxHashSet::default();
        seen_hashes.insert(hash);

        // The priority queue of graphs to be processed (this should not get big)
        let mut pq = DoublePriorityQueue::with_capacity(opt.queue_size);
        pq.push(
            Entry {
                g: best_g.clone(),
                hash,
                cost: best_g_cost,
            },
            best_g_cost,
        );

        let mut g_cnt = 0;
        let mut timeout_flag = false;
        let mut last_log_time = Instant::now();
        let mut last_log_cnt = g_cnt;
        while let Some((Entry { g, cost, .. }, _)) = pq.pop_min() {
            if cost < best_g_cost {
                best_g = g.clone();
                best_g_cost = cost;
                log::info!("New best cost: {:?}", best_g_cost);
                last_best_time = Instant::now();
            }
            g_cnt += 1;

            let rewrites = self.rewriter.get_rewrites(&g);

            // Get combinations of rewrites that can be applied to the graph,
            // and filter them to keep only the ones that
            //
            // - Don't have a worse cost than the last candidate in the priority queue.
            // - Do not invalidate the graph by creating a loop.
            // - We haven't seen yet.
            for rw in rewrites {
                let r = self.rewriter.apply_rewrite(rw, &g);
                let new_g_cost = cost.saturating_add_signed(r.cost_delta);

                // Skip rewrites that have a worse cost than the last candidate in the priority queue.
                let max_cost = pq.peek_max().map(|(_, c)| *c).unwrap_or(usize::MAX);
                if new_g_cost < max_cost {
                    continue;
                }

                let new_g_hash = weisfeiler_lehman_graph_hash(&r.graph, 3);

                if !seen_hashes.insert(new_g_hash) {
                    // Ignore this graph: we've already seen it
                    continue;
                }

                // Queue the new graph.
                pq.push(
                    Entry {
                        g: r.graph,
                        hash: new_g_hash,
                        cost: new_g_cost,
                    },
                    new_g_cost,
                );
                if pq.len() > opt.queue_size {
                    pq.pop_max();
                }

                // Log progress
                if g_cnt > last_log_cnt && Instant::now() - last_log_time > Duration::from_secs(1) {
                    last_log_cnt = g_cnt;
                    last_log_time = Instant::now();
                    log::info!(
                        "Processed: {}  Queue size: {}  Seen hash count: {}",
                        g_cnt,
                        pq.len(),
                        seen_hashes.len()
                    );
                }
            }

            if let Some(timeout) = opt.timeout {
                if start_time.elapsed().as_secs() > timeout {
                    timeout_flag = true;
                    break;
                }
            }
            if let Some(p_timeout) = opt.progress_timeout {
                if last_best_time.elapsed().as_secs() > p_timeout {
                    timeout_flag = true;
                    break;
                }
            }
        }

        let timeout_str = match timeout_flag {
            true => " due to timeout",
            false => "",
        };
        log::info!("Optimization finished{timeout_str}.");
        log::info!(
            "Processed {} circuits (out of {} seen).",
            g_cnt,
            seen_hashes.len()
        );
        log::info!("Final cost: {:?}", best_g_cost);
        best_g
    }
}

/// An entry in the priority queue.
#[derive(Debug, Clone)]
struct Entry<G> {
    pub g: G,
    pub hash: u64,
    pub cost: usize,
}

impl<G> Hash for Entry<G> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hash.hash(state);
    }
}

impl<G> PartialEq for Entry<G> {
    fn eq(&self, other: &Self) -> bool {
        self.hash == other.hash
    }
}

impl<G> Eq for Entry<G> {}

#[cfg(test)]
mod test {
    use super::*;
    use crate::rewrite_sets::test::rewrite_set_2qb_lc;
    use crate::rewrite_sets::RewriteSet;
    use crate::rewriter::test::{json_simple_graph, small_graph};
    use crate::rewriter::CausalRewriter;

    use quizx::vec_graph::Graph;
    use rstest::rstest;

    #[rstest]
    #[case(small_graph())]
    #[case(json_simple_graph())]
    fn run_superopt(
        rewrite_set_2qb_lc: Vec<RewriteSet<Graph>>,
        #[case] graph: Graph,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let rewriter = CausalRewriter::from_rewrite_rules(rewrite_set_2qb_lc);
        let optimizer: SuperOptimizer<CausalRewriter<Graph>> = SuperOptimizer::new(rewriter);
        let options = SuperOptOptions {
            timeout: Some(5),
            progress_timeout: Some(1),
            ..Default::default()
        };

        let new_g = optimizer.optimize(&graph, options);

        assert!(new_g.num_edges() <= graph.num_edges());

        Ok(())
    }
}
