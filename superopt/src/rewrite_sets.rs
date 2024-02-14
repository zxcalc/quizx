//! This module defines the serializable definitions for sets of causal flow
//! preserving ZX rewrite rules.
//!
//! See https://github.com/CQCL-DEV/zx-causal-flow-rewrites for a generator of
//! these sets.

use std::path::Path;

use quizx::json::{decode_graph, encode_graph};
use quizx::vec_graph::GraphLike;
use serde::{de, Deserialize, Deserializer, Serialize};

/// Reads a graph from a json-encoded list of rewrite rule sets.
pub fn read_rewrite_sets<G: GraphLike + for<'de> Deserialize<'de>>(
    filename: &Path,
) -> serde_json::Result<G> {
    let file = std::fs::File::open(filename).unwrap();
    let reader = std::io::BufReader::new(file);
    serde_json::from_reader(reader)
}

/// Writes the json-encoded representation of a list of rewrite rule sets.
pub fn write_rewrite_sets<G: GraphLike + Serialize>(
    rule_sets: &[RewriteSet<G>],
    filename: &Path,
) -> serde_json::Result<()> {
    let file = std::fs::File::create(filename).unwrap();
    let writer = std::io::BufWriter::new(file);
    serde_json::to_writer(writer, rule_sets)
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RewriteSet<G: GraphLike> {
    /// Left hand side of the rewrite rule
    pub lhs: RewriteLhs<G>,
    /// Possible input/output assignments of the boundary nodes
    pub lhs_ios: Vec<RewriteLhsIos>,
    /// List of possible right hand sides of the rewrite rule
    pub rhss: Vec<RewriteRhs<G>>,
}

/// The left hand side of a rewrite rule
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct RewriteLhs<G: GraphLike> {
    /// Replacement graph
    #[serde(deserialize_with = "de_graph_from_str")]
    #[serde(serialize_with = "ser_graph_to_str")]
    pub g: G,
}

/// Possible input/output assignments of the boundary nodes
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RewriteLhsIos(Vec<usize>, Vec<usize>);

impl RewriteLhsIos {
    pub fn new(inputs: Vec<usize>, outputs: Vec<usize>) -> Self {
        Self(inputs, outputs)
    }

    pub fn inputs(&self) -> &[usize] {
        &self.0
    }

    pub fn outputs(&self) -> &[usize] {
        &self.1
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RewriteRhs<G: GraphLike> {
    /// Two-qubit gate reduction over the LHS
    pub reduction: isize,
    /// Replacement graph
    #[serde(deserialize_with = "de_graph_from_str")]
    #[serde(serialize_with = "ser_graph_to_str")]
    pub g: G,
    /// If the rewrite is a local complementation, the list of unfused vertex indices
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(default)]
    pub unfused: Option<Vec<usize>>,
    /// If the rewrite is a pivot, the list of unfused vertex indices for the first pivot vertex
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(default)]
    pub unfused1: Option<Vec<usize>>,
    /// If the rewrite is a pivot, the list of unfused vertex indices for the second pivot vertex
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(default)]
    pub unfused2: Option<Vec<usize>>,
}

/// Deserialize a graph from a string field in the JSON.
fn de_graph_from_str<'de, D, G>(deserializer: D) -> Result<G, D::Error>
where
    D: Deserializer<'de>,
    G: GraphLike,
{
    let s = String::deserialize(deserializer)?;
    decode_graph(&s).map_err(de::Error::custom)
}

/// Serialize a graph to a string field in the JSON.
fn ser_graph_to_str<S, G>(graph: &G, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
    G: GraphLike,
{
    let s = encode_graph(graph, true).unwrap();
    s.serialize(serializer)
}

#[cfg(test)]
mod test {
    use super::*;
    use quizx::vec_graph::Graph;

    const TEST_SET: &str = include_str!("../../test_files/rewrites-2qb-lc.json");

    #[test]
    fn test_rewrite_set_serde() {
        let rewrite_sets: Vec<RewriteSet<Graph>> = serde_json::from_str(TEST_SET).unwrap();

        assert_eq!(rewrite_sets.len(), 3);
    }
}
