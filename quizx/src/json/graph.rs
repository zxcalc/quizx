// QuiZX - Rust library for quantum circuit rewriting and optimisation
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

//! Methods for converting between a `GraphLike` object and the json representation.

use num::{Rational, Zero};

use super::{
    EdgeAttrs, JsonGraph, VertexAnnotations, VertexAttrs, VertexData, VertexName, VertexPhase,
};
use crate::graph::{Coord, EType, GraphLike, VData, VType, V};

use std::collections::{BTreeMap, HashMap};

impl JsonGraph {
    /// Encode a graph using the json representation.
    pub fn from_graph(graph: &impl GraphLike, ignore_scalar: bool) -> Self {
        let mut wire_vertices = HashMap::new();
        let mut node_vertices = HashMap::new();
        let mut undir_edges = HashMap::new();

        // The encoding requires unique string names for vertices and edges.
        let mut vertex_name_gen = (0..).map(|i| format!("v{}", i));
        let mut bound_name_gen = (0..).map(|i| format!("b{}", i));
        let mut edge_name_gen = (0..).map(|i| format!("e{}", i));

        let mut v_names: HashMap<V, VertexName> = HashMap::new();

        for v in graph.vertices() {
            let typ = graph.vertex_type(v);
            let coord = graph.coord(v).to_f64();
            let v_name = match typ {
                VType::B => bound_name_gen.next(),
                _ => vertex_name_gen.next(),
            }
            .unwrap();
            v_names.insert(v, v_name.clone());

            if typ == VType::B {
                let input = graph.inputs().iter().position(|&i| i == v);
                let output = graph.outputs().iter().position(|&o| o == v);
                assert!(
                    input.is_some() || output.is_some(),
                    "Boundary vertex is not an input nor output."
                );
                let attrs = VertexAttrs {
                    annotations: VertexAnnotations {
                        boundary: true,
                        coord,
                        input,
                        output,
                        ..Default::default()
                    },
                    data: VertexData {
                        typ,
                        ..Default::default()
                    },
                };

                wire_vertices.insert(v_name, attrs);
            } else {
                let phase = graph.phase(v);
                let mut attrs = VertexAttrs {
                    annotations: VertexAnnotations {
                        coord,
                        ..Default::default()
                    },
                    data: VertexData {
                        typ,
                        value: VertexPhase::from_rational(phase, typ),
                        ..Default::default()
                    },
                };

                if typ == VType::ZBox {
                    // Data for ZBox vertices is not currently supported, so write the default.
                    attrs.annotations.label = Some("1".to_string());
                }

                node_vertices.insert(v_name, attrs);
            };
        }

        for (src, tgt, typ) in graph.edges() {
            match typ {
                EType::N | EType::Wio => {
                    let attr = EdgeAttrs {
                        src: v_names[&src].clone(),
                        tgt: v_names[&tgt].clone(),
                        typ,
                    };
                    undir_edges.insert(edge_name_gen.next().unwrap(), attr);
                }
                EType::H => {
                    // Encoded as a Hadamard node and two simple edges.
                    let h_name = vertex_name_gen.next().unwrap();
                    let coord = avg_coord(graph.coord(src), graph.coord(tgt));
                    node_vertices.insert(
                        h_name.clone(),
                        VertexAttrs {
                            annotations: VertexAnnotations {
                                coord,
                                ..Default::default()
                            },
                            data: VertexData {
                                typ: VType::H,
                                is_edge: Some(true),
                                ..Default::default()
                            },
                        },
                    );
                    undir_edges.insert(
                        edge_name_gen.next().unwrap(),
                        EdgeAttrs {
                            src: v_names[&src].clone(),
                            tgt: h_name.clone(),
                            ..Default::default()
                        },
                    );
                    undir_edges.insert(
                        edge_name_gen.next().unwrap(),
                        EdgeAttrs {
                            src: h_name,
                            tgt: v_names[&tgt].clone(),
                            ..Default::default()
                        },
                    );
                }
            };
        }

        if !ignore_scalar {
            unimplemented!("Encoding scalars is not yet supported. Use ignore_scalar=true.");
        }

        Self {
            wire_vertices,
            node_vertices,
            undir_edges,
            variable_types: vec![],
            scalar: None,
        }
    }

    /// Decode a graph from the json representation.
    pub fn to_graph<G: GraphLike>(&self) -> G {
        let mut graph = G::new();

        if !self.variable_types.is_empty() {
            unimplemented!("Variables are not currently supported.");
        }

        if self.scalar.is_some() {
            unimplemented!("Scalars are not currently supported.");
        }

        let mut names: HashMap<VertexName, V> = HashMap::new();

        // Map used to track auxiliary Hadamard nodes that should be decoded as Hadamard edges.
        // Stores the neighbor nodes of the Hadamard node.
        let mut hadamards: HashMap<&str, Vec<V>> = HashMap::new();

        for (name, attrs) in &self.node_vertices {
            if attrs.data.typ == VType::H && attrs.data.is_edge.unwrap_or_default() {
                // A virtual hadamard edge.
                hadamards.insert(name, vec![]);
                continue;
            }
            let coord = Coord::from_f64(attrs.annotations.coord);

            let v = graph.add_vertex_with_data(VData {
                ty: attrs.data.typ,
                qubit: coord.qubit(),
                row: coord.row(),
                phase: attrs.data.value.to_rational().unwrap_or(Rational::zero()),
            });
            names.insert(name.to_string(), v);
        }

        // Insert the boundary nodes, and collect the input and output vectors.
        let mut inputs: BTreeMap<usize, &str> = BTreeMap::new();
        let mut outputs: BTreeMap<usize, &str> = BTreeMap::new();
        for (name, attrs) in &self.wire_vertices {
            let coord = Coord::from_f64(attrs.annotations.coord);
            let v = graph.add_vertex_with_data(VData {
                ty: VType::B,
                qubit: coord.qubit(),
                row: coord.row(),
                phase: Rational::zero(),
            });
            names.insert(name.to_string(), v);
            if let Some(input) = attrs.annotations.input {
                inputs.insert(input, name);
            }
            if let Some(output) = attrs.annotations.output {
                outputs.insert(output, name);
            }
        }
        graph.set_inputs(inputs.into_values().map(|name| names[name]).collect());
        graph.set_outputs(outputs.into_values().map(|name| names[name]).collect());

        println!("Inputs: {:?}", graph.inputs());
        println!("Outputs: {:?}", graph.outputs());

        // Insert the edges.
        for attrs in self.undir_edges.values() {
            let src = names[&attrs.src];
            let tgt = names[&attrs.tgt];

            match (
                hadamards.contains_key(attrs.src.as_str()),
                hadamards.contains_key(attrs.tgt.as_str()),
            ) {
                (true, true) => {
                    // Both ends are virtual Hadamard nodes.
                    //
                    // Not sure how this is possible, but pyzx supports it.
                    let new_coord = Coord::from_f64(avg_coord(graph.coord(src), graph.coord(tgt)));
                    let v = graph.add_vertex_with_data(VData {
                        ty: VType::Z,
                        qubit: new_coord.qubit(),
                        row: new_coord.row(),
                        phase: Rational::zero(),
                    });
                    let name = format!("v{}", graph.num_vertices());
                    names.insert(name, v);
                    hadamards.get_mut(attrs.src.as_str()).unwrap().push(v);
                    hadamards.get_mut(attrs.tgt.as_str()).unwrap().push(v);
                    continue;
                }
                (true, false) => {
                    hadamards.get_mut(attrs.src.as_str()).unwrap().push(tgt);
                    continue;
                }
                (false, true) => {
                    hadamards.get_mut(attrs.tgt.as_str()).unwrap().push(src);
                    continue;
                }
                _ => {}
            }

            graph.add_edge_smart(src, tgt, attrs.typ);
        }

        // Add the Hadamard edges.
        for (_, neighbors) in hadamards {
            if neighbors.len() != 2 {
                panic!("Virtual Hadamard node has wrong number of neighbors.");
            }
            let (src, tgt) = (neighbors[0], neighbors[1]);
            graph.add_edge_smart(src, tgt, EType::H);
        }

        graph
    }
}

/// Returns the average of two coordinates, as a pair of f64.
///
/// Rounds the result to 3 decimal places.
fn avg_coord(a: Coord, b: Coord) -> (f64, f64) {
    let a = a.to_f64();
    let b = b.to_f64();
    (
        ((a.0 + b.0) * 1000.).round() / 2000.,
        ((a.1 + b.1) * 1000.).round() / 2000.,
    )
}