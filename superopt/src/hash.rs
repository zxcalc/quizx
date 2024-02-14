use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use itertools::Itertools;
use quizx::graph::{GraphLike, V};
use rustc_hash::FxHashMap;

/// Computes the hash of a graph using the Weisfeiler-Lehman algorithm.
///
/// The algorithm iteratively hashes each vertex based on the hash of its neighbors.
/// This is done for `num_iterations` iterations. This is normally set to 3.
///
/// The final hash is the hash of the sorted list of vertex hashes.
pub fn weisfeiler_lehman_graph_hash(graph: &impl GraphLike, num_iterations: usize) -> u64 {
    let mut vertex_hashes: FxHashMap<V, u64> = FxHashMap::default();

    // Compute the initial hash of each vertex from its data.
    for vertex in graph.vertices() {
        let data = graph.vertex_data(vertex);
        let mut vhasher = DefaultHasher::new();
        data.hash(&mut vhasher);
        vertex_hashes.insert(vertex, vhasher.finish());
    }

    // Repeatedly compute the hash of each vertex based on the hash of its
    // neighbors.
    for _ in 0..num_iterations {
        let mut new_hashes: FxHashMap<V, u64> = FxHashMap::default();

        for vertex in graph.vertices() {
            let mut state = DefaultHasher::new();
            state.write_u64(vertex_hashes[&vertex]);

            for (neighbor, etype) in graph.incident_edges(vertex) {
                state.write_u64(vertex_hashes[&neighbor]);
                etype.hash(&mut state);
            }

            new_hashes.insert(vertex, state.finish());
        }

        vertex_hashes = new_hashes;
    }

    // The final hash is the hash of the sorted list of vertex hashes.
    let mut hasher = DefaultHasher::new();
    vertex_hashes
        .into_values()
        .sorted_unstable()
        .for_each(|h| hasher.write_u64(h));

    hasher.finish()
}
