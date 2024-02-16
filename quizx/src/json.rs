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

//! Json encoding for interoperability with pyzx.

mod graph;
mod phase;

use crate::graph::{EType, GraphLike, VType};

use serde::{de, Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Returns the json-encoded representation of a graph.
pub fn encode_graph(
    graph: &impl crate::graph::GraphLike,
    ignore_scalar: bool,
) -> serde_json::Result<String> {
    let jg = JsonGraph::from_graph(graph, ignore_scalar);
    serde_json::to_string(&jg)
}

/// Writes the json-encoded representation of a graph to a file.
pub fn write_graph(
    graph: &impl crate::graph::GraphLike,
    ignore_scalar: bool,
    filename: &Path,
) -> serde_json::Result<()> {
    let jg = JsonGraph::from_graph(graph, ignore_scalar);
    let file = std::fs::File::create(filename).unwrap();
    let writer = std::io::BufWriter::new(file);
    serde_json::to_writer(writer, &jg)
}

/// Reads a graph from its json-encoded representation.
pub fn decode_graph<G: GraphLike>(s: &str) -> serde_json::Result<G> {
    let jg: JsonGraph = serde_json::from_str(s)?;
    Ok(jg.to_graph(true).0)
}

/// Reads a graph from a json-encoded file.
pub fn read_graph<G: GraphLike>(filename: &Path) -> serde_json::Result<G> {
    let file = std::fs::File::open(filename).unwrap();
    let reader = std::io::BufReader::new(file);
    let jg: JsonGraph = serde_json::from_reader(reader)?;
    Ok(jg.to_graph(true).0)
}

/// Identifier for an encoded vertex.
pub type VertexName = String;
/// Identifier for an encoded edge.
type EdgeName = String;

/// The json-encoded format for pyzx graphs.
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
    variable_types: HashMap<String, String>,
    /// The graph scalar.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    scalar: Option<String>,
}

/// Attributes for a vertex in the json-encoded graph.
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
struct VertexAttrs {
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    annotation: VertexAnnotations,
    #[serde(skip_serializing_if = "is_default")]
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
    value: VertexPhase,
    /// A flag marking grounded nodes.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    ground: bool,
    /// Hadamard wires are encoded as nodes with this flag set,
    /// so that they can be recovered during decoding.
    ///
    /// Note that in the pyzx encoding, this is either the string "true" or "false".
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
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    coord: (f64, f64),
    /// The input number for the vertex.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
    input: Option<usize>,
    /// The output number for the vertex.
    #[serde(skip_serializing_if = "is_default")]
    #[serde(default)]
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

/// The annotations of a vertex in the json-encoded graph.
#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(transparent)]
struct VertexPhase(String);

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
    let s: String = de::Deserialize::deserialize(deserializer)?;

    match s.as_str() {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(de::Error::unknown_variant(&s, &["true", "false"])),
    }
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

    #[fixture]
    fn hadamard_edge() -> (Graph, Vec<V>) {
        let mut g = Graph::new();
        let vs = vec![
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
        ];
        g.add_edge_with_type(vs[1], vs[2], EType::H);
        g.remove_vertex(vs[0]);
        (g, vs)
    }

    const TEST_JSON: &str = include_str!("../../test_files/simple-graph.json");

    #[rstest]
    #[case::simple_graph(simple_graph())]
    #[case::hadamard(hadamard_edge())]
    fn json_roundtrip(#[case] graph: (Graph, Vec<V>)) {
        let (g, _) = graph;
        let jg = JsonGraph::from_graph(&g, true);
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
    }

    #[test]
    fn json_decode() {
        let g: Graph = decode_graph(TEST_JSON).unwrap();

        assert_eq!(g.num_vertices(), 9);
        assert_eq!(g.num_edges(), 9);
    }
}
