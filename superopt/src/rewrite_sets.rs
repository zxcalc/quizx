//! This module defines the serializable definitions for sets of causal flow
//! preserving ZX rewrite rules.
//!
//! See https://github.com/CQCL-DEV/zx-causal-flow-rewrites for a generator of
//! these sets.

use std::collections::HashMap;
use std::path::Path;

use quizx::json::{JsonGraph, VertexName};
use quizx::vec_graph::{GraphLike, V};
use serde::{Deserialize, Deserializer, Serialize};

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

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RewriteSet<G: GraphLike> {
    /// Left hand side of the rewrite rule
    pub lhs: DecodedGraph<G>,
    /// Possible input/output assignments of the boundary nodes
    pub lhs_ios: Vec<RewriteIos>,
    /// List of possible right hand sides of the rewrite rule
    pub rhss: Vec<RewriteRhs<G>>,
}

impl<G: GraphLike> RewriteSet<G> {
    /// Returns the input/output assignments of the boundary nodes of the LHS,
    /// translated to the graph indices.
    pub fn lhs_ios_translated(&self) -> impl Iterator<Item = (Vec<V>, Vec<V>)> + '_ {
        self.lhs_ios
            .iter()
            .map(move |ios| ios.translated(&self.lhs.names))
    }

    pub fn lhs_boundary(&self) -> Vec<V> {
        todo!()
    }
}

/// Possible input/output assignments of the boundary nodes
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RewriteIos(Vec<String>, Vec<String>);

impl RewriteIos {
    pub fn new(inputs: Vec<String>, outputs: Vec<String>) -> Self {
        Self(inputs, outputs)
    }

    pub fn inputs(&self) -> &[String] {
        &self.0
    }

    pub fn outputs(&self) -> &[String] {
        &self.1
    }

    pub fn translated(&self, names: &HashMap<VertexName, V>) -> (Vec<V>, Vec<V>) {
        (
            self.0.iter().map(|name| names[name]).collect(),
            self.1.iter().map(|name| names[name]).collect(),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RewriteRhs<G: GraphLike> {
    /// Two-qubit gate reduction over the LHS
    pub reduction: isize,
    /// Replacement graph
    pub g: DecodedGraph<G>,
    /// Possible input/output assignments of the boundary nodes
    pub ios: Vec<RewriteIos>,
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

impl<G: GraphLike> RewriteRhs<G> {
    pub fn boundary(&self) -> Vec<V> {
        todo!()
    }

    pub fn ios_translated(&self) -> impl Iterator<Item = (Vec<V>, Vec<V>)> + '_ {
        self.ios
            .iter()
            .map(move |ios| ios.translated(&self.g.names))
    }
}

/// A decoded graph with a map from serialized vertex names to indices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedGraph<G: GraphLike> {
    pub g: G,
    pub names: HashMap<VertexName, V>,
}

impl<'de, G: GraphLike> Deserialize<'de> for DecodedGraph<G> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s: String = Deserialize::deserialize(deserializer)?;
        let jg: JsonGraph = serde_json::from_str(&s).unwrap(); // TODO: error handling
        let (g, names) = jg.to_graph(true);
        Ok(DecodedGraph { g, names })
    }
}

impl<G: GraphLike> Serialize for DecodedGraph<G> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let jg = JsonGraph::from_graph(&self.g, true);
        let s = serde_json::to_string(&jg).map_err(serde::ser::Error::custom)?;
        s.serialize(serializer)
    }
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
