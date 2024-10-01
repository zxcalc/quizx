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

//! Json encoding for interoperability with pyzx and zxlang using the .qgraph format.
//!
//! # Examples
//!
//! ```rust
//! # use quizx::graph::{EType, VType, GraphLike};
//! # use quizx::vec_graph::Graph;
//! # use quizx::tensor::ToTensor;
//! // Define a graph with 4 vertices and 3 edges.
//! let mut g = Graph::new();
//! let vs = vec![
//!     g.add_vertex(VType::B),
//!     g.add_vertex(VType::Z),
//!     g.add_vertex(VType::X),
//!     g.add_vertex(VType::B),
//! ];
//! g.set_inputs(vec![vs[0]]);
//! g.set_outputs(vec![vs[3]]);
//! g.add_edge(vs[0], vs[1]);
//! g.add_edge_with_type(vs[1], vs[2], EType::H);
//! g.add_edge(vs[2], vs[3]);
//!
//! // Encode the graph in qgraph format.
//! let json = quizx::json::encode_graph(&g).unwrap();
//!
//! // Decode the graph from the qgraph string.
//! let g2 = quizx::json::decode_graph::<Graph>(&json).unwrap();
//!
//! assert_eq!(g.to_tensor4(), g2.to_tensor4());
//! ```

mod graph;
mod phase;
mod scalar;

pub use phase::PhaseOptions;

use crate::graph::VType;
use crate::hash_graph::{EType, GraphLike};

use derive_more::{Display, Error, From};
use serde::{de, Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Returns the json-encoded representation of a graph.
pub fn encode_graph(graph: &impl crate::graph::GraphLike) -> Result<String, JsonError> {
    let jg = JsonGraph::from_graph(graph)?;
    let s = serde_json::to_string(&jg)?;
    Ok(s)
}

/// Writes the json-encoded representation of a graph to a file.
pub fn write_graph(graph: &impl crate::graph::GraphLike, filename: &Path) -> Result<(), JsonError> {
    let jg = JsonGraph::from_graph(graph)?;
    let file = std::fs::File::create(filename).unwrap();
    let writer = std::io::BufWriter::new(file);
    serde_json::to_writer(writer, &jg)?;
    Ok(())
}

/// Reads a graph from its json-encoded representation.
pub fn decode_graph<G: GraphLike>(s: &str) -> Result<G, JsonError> {
    let jg: JsonGraph = serde_json::from_str(s)?;
    jg.to_graph()
}

/// Reads a graph from a json-encoded file.
pub fn read_graph<G: GraphLike>(filename: &Path) -> Result<G, JsonError> {
    let file = std::fs::File::open(filename).unwrap();
    let reader = std::io::BufReader::new(file);
    let jg: JsonGraph = serde_json::from_reader(reader)?;
    jg.to_graph()
}

/// Identifier for an encoded vertex.
type VertexName = String;
/// Identifier for an encoded edge.
type EdgeName = String;

/// The json-encoded format for pyzx and zxlang graphs.
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
pub struct JsonGraph {
    /// Wire vertices of the graph.
    #[serde(default)]
    wire_vertices: HashMap<VertexName, VertexAttrs>,
    /// Node vertices of the graph.
    #[serde(default)]
    node_vertices: HashMap<VertexName, VertexAttrs>,
    /// Undirected edges between node vertices.
    #[serde(default)]
    undir_edges: HashMap<EdgeName, EdgeAttrs>,
    /// Types of the variables in the graph.
    ///
    /// Currently ignored by quizx.
    #[serde(default)]
    variable_types: HashMap<String, String>,
    /// The graph scalar.
    ///
    /// pyzx encodes this as a json-encoded string instead of directly embedding
    /// the dictionary.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_scalar")]
    #[serde(serialize_with = "serialize_scalar")]
    scalar: Option<JsonScalar>,
}

/// Attributes for a vertex in the json-encoded graph.
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
struct VertexAttrs {
    #[serde(default)]
    annotation: VertexAnnotations,
    #[serde(default)]
    data: VertexData,
}

/// Data for a vertex in the json-encoded graph.
#[derive(Serialize, Deserialize, Debug, Default, Clone, PartialEq, Eq)]
struct VertexData {
    /// The vertex type.
    #[serde(rename = "type")]
    typ: VType,
    /// The vertex phase.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    value: JsonPhase,
    /// A flag marking grounded nodes.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    ground: bool,
    /// Hadamard wires are encoded as nodes with this flag set,
    /// so that they can be recovered during decoding.
    ///
    /// Note that in the pyzx encoder, this is either the string "true" or "false".
    /// So we need to deserialize it into a bool manually.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_bool")]
    #[serde(serialize_with = "serialize_bool")]
    is_edge: bool,
}

/// The annotations of a vertex in the json-encoded graph.
#[derive(Serialize, Deserialize, Debug, Default, Clone, PartialEq)]
struct VertexAnnotations {
    /// This is a boundary vertex.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    boundary: bool,
    /// Coordinates of the vertex, for rendering purposes.
    #[serde(default)]
    coord: (f64, f64),
    /// The input number for the vertex.
    ///
    /// Note that in the pyzx encoder, this is either a number indicating the order in the input list, or
    /// (in older versions) a boolean flag indicating whether the vertex is an input.
    /// So we need to deserialize it manually.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_io_attribute")]
    input: Option<usize>,
    /// The output number for the vertex.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_io_attribute")]
    output: Option<usize>,
    /// A box label.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    label: Option<String>,
    /// Other arbitrary annotations associated with the vertex.
    #[serde(flatten)]
    other: HashMap<String, String>,
}

/// Attributes for an edge in the json-encoded graph.
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
struct EdgeAttrs {
    /// Edge source.
    src: VertexName,
    /// Edge target.
    tgt: VertexName,
    /// The edge type.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    #[serde(rename = "type")]
    typ: EType,
}

/// A phase, in half turns.
///
/// Encoded as string-formatted rational or a floating point number.
/// Any occurrence of "pi" or "π" in the string is ignored.
#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(transparent)]
pub struct JsonPhase(String);

/// Global scalars in a graph
#[derive(Serialize, Deserialize, Debug, Default, Clone, PartialEq)]
pub struct JsonScalar {
    /// Stores the power of sqrt(2).
    #[serde(default)]
    power2: i32,
    /// Stores complex phase of the number, in half turns.
    #[serde(default)]
    phase: JsonPhase,
    /// A floating point factor for the scalar.
    #[serde(default)]
    #[serde(skip_serializing_if = "is_default")]
    floatfactor: f64,
    /// Stores a list of legless spiders, by their phases.
    #[serde(default)]
    #[serde(skip_serializing_if = "is_default")]
    phasenodes: Vec<JsonPhase>,
    /// Whether this scalar is zero.
    #[serde(default)]
    #[serde(skip_serializing_if = "is_default")]
    is_zero: bool,
    /// Whether this represents an unknown scalar value.
    #[serde(default)]
    #[serde(skip_serializing_if = "is_default")]
    is_unknown: bool,
}

/// Helper method to skip serialization of default values in serde.
///
/// ```skip
/// use serde::Serialize;
///
/// #[derive(Serialize)]
/// struct MyStruct {
///     #[serde(skip_serializing_if = "is_default")]
///     field: i32,
/// }
/// ```
///
/// From https://github.com/serde-rs/serde/issues/818.
pub(crate) fn is_default<T: Default + PartialEq>(t: &T) -> bool {
    *t == Default::default()
}

/// Serialize a boolean to a string field.
fn serialize_bool<S>(b: &bool, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    match *b {
        true => serializer.serialize_str("true"),
        false => serializer.serialize_str("false"),
    }
}

/// Deserialize a boolean from a string field.
fn deserialize_bool<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: de::Deserializer<'de>,
{
    let s: &str = de::Deserialize::deserialize(deserializer)?;

    match s {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(de::Error::unknown_variant(s, &["true", "false"])),
    }
}

/// Deserialize the input/output attribute of a vertex.
///
/// This is either a number indicating the order in the input list, or
/// (in older versions) a boolean flag indicating whether the vertex is an input.
fn deserialize_io_attribute<'de, D>(deserializer: D) -> Result<Option<usize>, D::Error>
where
    D: de::Deserializer<'de>,
{
    let val: serde_json::Value = de::Deserialize::deserialize(deserializer)?;

    match val {
        serde_json::Value::Number(n) => Ok(Some(n.as_u64().unwrap() as usize)),
        serde_json::Value::Bool(b) => Ok(b.then_some(0)),
        _ => Err(de::Error::invalid_value(
            de::Unexpected::Str(&val.to_string()),
            &"a number or a boolean",
        )),
    }
}

/// Deserialize a scalar from a string field.
fn deserialize_scalar<'de, D>(deserializer: D) -> Result<Option<JsonScalar>, D::Error>
where
    D: de::Deserializer<'de>,
{
    let s: String = de::Deserialize::deserialize(deserializer)?;
    if s.is_empty() {
        return Ok(None);
    }
    let scalar = serde_json::from_str(&s).map_err(|_| {
        de::Error::invalid_value(
            de::Unexpected::Str(&s),
            &"a json-encoded string representing a scalar",
        )
    })?;
    Ok(Some(scalar))
}

/// Serialize a scalar to a string field.
fn serialize_scalar<S>(scalar: &Option<JsonScalar>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    if scalar.is_none() {
        return serializer.serialize_str("");
    }
    let s = serde_json::to_string(scalar.as_ref().unwrap()).map_err(serde::ser::Error::custom)?;
    serializer.serialize_str(&s)
}

/// An error that can occur when decoding a json graph.
#[derive(Debug, Display, Error, From)]
#[non_exhaustive]
pub enum JsonError {
    /// Found an invalid phase value.
    #[display("Got an invalid phase value {phase}")]
    InvalidPhase { phase: String },
    /// Found an invalid phase value in a node definition.
    #[display("Got an invalid phase value {phase} for node {name}")]
    InvalidNodePhase { name: String, phase: String },
    /// Some other serde error.
    #[from]
    SerdeError(serde_json::Error),
}

#[cfg(test)]
mod test {
    use crate::graph::GraphLike;
    use crate::vec_graph::{Graph, V};

    use super::*;

    use rstest::{fixture, rstest};

    /// Makes a simple graph.
    ///
    /// With `0` and `1` as inputs, and `6` and `7` as outputs.
    #[fixture]
    fn simple_graph() -> (Graph, Vec<V>) {
        let mut g = Graph::new();
        let vs = vec![
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::X),
            g.add_vertex(VType::X),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];

        g.set_inputs(vec![vs[0], vs[1]]);
        g.set_outputs(vec![vs[6], vs[7]]);

        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[3]);
        g.add_edge(vs[2], vs[4]);
        g.add_edge_with_type(vs[2], vs[3], EType::H);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[4], vs[6]);
        g.add_edge(vs[5], vs[7]);
        (g, vs)
    }

    //const TEST_JSON_SIMPLE: &str = include_str!("../../test_files/simple-graph.qgraph");
    const TEST_JSON_4Q_UNITARY: &str = include_str!("../../test_files/4-qubit-unitary.qgraph");

    #[rstest]
    fn json_roundtrip(simple_graph: (Graph, Vec<V>)) -> Result<(), JsonError> {
        let (g, _) = simple_graph;
        let jg = JsonGraph::from_graph(&g)?;
        let s = serde_json::to_string(&jg).unwrap();

        let g2: Graph = decode_graph(&s).unwrap();

        assert_eq!(g.num_vertices(), g2.num_vertices());
        assert_eq!(g.num_edges(), g2.num_edges());

        let inp1 = g.inputs();
        let inp2 = g2.inputs();

        assert_eq!(inp1.len(), inp2.len());
        for (v1, v2) in inp1.iter().zip(inp2.iter()) {
            assert_eq!(g.vertex_type(*v1), g2.vertex_type(*v2));

            let ns1 = g.neighbors(*v1);
            let ns2 = g2.neighbors(*v2);
            for (n1, n2) in ns1.zip(ns2) {
                assert_eq!(g.vertex_type(n1), g2.vertex_type(n2));
            }
        }

        Ok(())
    }

    #[rstest]
    //#[case::simple(TEST_JSON_SIMPLE, 9, 9)]  TODO: Parameters are not yet supported
    #[case::unitary_4q(TEST_JSON_4Q_UNITARY, 26, 30)]
    fn json_decode(#[case] json: &str, #[case] num_vertices: usize, #[case] num_edges: usize) {
        let g: Graph = decode_graph(json).unwrap();

        assert_eq!(g.num_vertices(), num_vertices);
        assert_eq!(g.num_edges(), num_edges);
    }
}
