Causal flow preserving superoptimizer
-------------------------------------

Causal flow preserving superoptimizer for ZX diagrams.

Based on the paper "Exhaustive Generation of Quantum Circuit Rewrite Rules using ZX" [arxiv][].
Adapted from the [Badger superoptimizer] from `tket2`.

Using a rewriter and a rewrite strategy, the optimizer
will repeatedly rewrite the graph, optimizing the graph according to
the cost function provided.

Optimization is done by maintaining a priority queue of graphs and
always processing the graph with the lowest cost first. Rewrites are
computed for that graph and all new graph obtained are added to the queue.

  [arxiv]: TODO
  [Badger superoptimizer]: https://github.com/CQCL/tket2