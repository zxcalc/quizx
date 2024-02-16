//! This module defines the serializable definitions for sets of causal flow
//! preserving ZX rewrite rules.
//!
//! See https://github.com/CQCL-DEV/zx-causal-flow-rewrites for a generator of
//! these sets.

use std::collections::HashMap;
use std::path::Path;

use itertools::Itertools;
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
    lhs: DecodedGraph<G>,
    /// Possible input/output assignments of the boundary nodes
    lhs_ios: Vec<RewriteIos>,
    /// List of possible right hand sides of the rewrite rule
    rhss: Vec<RewriteRhs<G>>,
}

/// Possible input/output assignments of the boundary nodes
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RewriteIos(Vec<String>, Vec<String>);

/// Auxiliary data structure for the left hand side of the rewrite rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RewriteLhs<'a, G: GraphLike> {
    /// Decoded graph representation of the left hand side of the rewrite rule
    g: &'a DecodedGraph<G>,
    /// Possible input/output assignments of the boundary nodes
    ios: &'a Vec<RewriteIos>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RewriteRhs<G: GraphLike> {
    /// Two-qubit gate reduction over the LHS
    pub reduction: isize,
    /// Replacement graph
    g: DecodedGraph<G>,
    /// Possible input/output assignments of the boundary nodes
    ios: Vec<RewriteIos>,
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

/// A decoded graph with a map from serialized vertex names to indices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedGraph<G: GraphLike> {
    pub g: G,
    names: HashMap<VertexName, V>,
}

impl<G: GraphLike> RewriteSet<G> {
    /// Returns the left hand side of the rewrite rule.
    pub fn lhs(&self) -> RewriteLhs<'_, G> {
        RewriteLhs::new(&self.lhs, &self.lhs_ios)
    }

    /// Returns the list of possible right hand sides of the rewrite rule.
    pub fn rhss(&self) -> &[RewriteRhs<G>] {
        &self.rhss
    }
}

impl<'a, G: GraphLike> RewriteLhs<'a, G> {
    pub fn new(g: &'a DecodedGraph<G>, ios: &'a Vec<RewriteIos>) -> Self {
        Self { g, ios }
    }
}

impl RewriteIos {
    pub fn new(inputs: Vec<String>, outputs: Vec<String>) -> Self {
        Self(inputs, outputs)
    }

    pub fn translated<G: GraphLike>(&self, g: &DecodedGraph<G>) -> (Vec<V>, Vec<V>) {
        let map_v = |name| {
            let v = g.from_name(name);
            unique_neighbour(&g.g, v)
        };
        (
            self.0.iter().map(map_v).collect(),
            self.1.iter().map(map_v).collect(),
        )
    }
}

impl<G: GraphLike> DecodedGraph<G> {
    pub fn name(&self, v: V) -> &VertexName {
        self.names
            .iter()
            .find(|(_, &idx)| idx == v)
            .map(|(name, _)| name)
            .unwrap_or_else(|| panic!("Vertex index {v} not found"))
    }

    pub fn from_name(&self, name: &VertexName) -> V {
        *self
            .names
            .get(name)
            .unwrap_or_else(|| panic!("Vertex name {name} not found"))
    }
}

/// Trait generalizing common operations between the LHS and RHS of a rewrite rule.
pub trait RuleSide<G: GraphLike> {
    /// The decoded graph representation of the rule side.
    fn decoded_graph(&self) -> &DecodedGraph<G>;

    /// The encoded input/output assignments of the boundary nodes.
    fn decoded_ios(&self) -> &[RewriteIos];

    /// The graph representation of the LHS/RHS
    ///
    /// Unlike the decoded graphs, these graphs have no quizx boundary vertices
    /// and are computed on the fly.
    fn graph(&self) -> G {
        let mut g = self.decoded_graph().g.clone();
        let io = g.inputs().iter().chain(g.outputs()).copied().collect_vec();
        for v in io {
            g.remove_vertex(v);
        }
        g.set_inputs(Vec::new());
        g.set_outputs(Vec::new());
        g
    }

    /// The boundary nodes of the graph in the LHS/RHS sense.
    ///
    /// These are vertices within `self.graph`, not the quizx boundary vertices.
    fn boundary<'a>(&'a self) -> impl Iterator<Item = V> + 'a
    where
        G: 'a,
    {
        let g = &self.decoded_graph().g;
        let inputs = g.inputs().as_slice();
        let outputs = g.outputs().as_slice();

        inputs
            .iter()
            .chain(outputs.iter())
            .map(|&v| unique_neighbour(g, v))
    }

    /// The input/output assignments of the boundary nodes, translated to the graph indices.
    fn ios(&self) -> impl Iterator<Item = (Vec<V>, Vec<V>)> + '_ {
        self.decoded_ios()
            .iter()
            .map(move |ios| ios.translated(self.decoded_graph()))
    }
}

impl<'a, G: GraphLike> RuleSide<G> for RewriteLhs<'a, G> {
    fn decoded_graph(&self) -> &DecodedGraph<G> {
        &self.g
    }

    fn decoded_ios(&self) -> &[RewriteIos] {
        self.ios
    }
}

impl<G: GraphLike> RuleSide<G> for RewriteRhs<G> {
    fn decoded_graph(&self) -> &DecodedGraph<G> {
        &self.g
    }

    fn decoded_ios(&self) -> &[RewriteIos] {
        &self.ios
    }
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

fn unique_neighbour(g: &impl GraphLike, v: V) -> V {
    g.neighbors(v)
        .exactly_one()
        .unwrap_or_else(|_| panic!("Boundary node {} has {} neighbors", v, g.neighbors(v).len()))
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;
    use quizx::vec_graph::Graph;
    use rstest::fixture;

    const TEST_SET: &str = include_str!("../../test_files/rewrites-2qb-lc.json");

    #[fixture]
    pub(crate) fn rewrite_set_2qb_lc() -> Vec<RewriteSet<Graph>> {
        serde_json::from_str(TEST_SET).unwrap()
    }

    #[test]
    fn test_rewrite_set_serde() {
        let rewrite_sets: Vec<RewriteSet<Graph>> = serde_json::from_str(TEST_SET).unwrap();

        assert_eq!(rewrite_sets.len(), 3);

        for set in rewrite_sets {
            let lhs = set.lhs();
            for rhs in set.rhss() {
                assert_eq!(lhs.boundary().count(), rhs.boundary().count());
            }
        }
    }
}
