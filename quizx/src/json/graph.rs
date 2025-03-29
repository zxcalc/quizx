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

use num::{One, Rational64, Zero};

use super::phase::PhaseOptions;
use super::{
    EdgeAttrs, JsonError, JsonGraph, JsonPhase, JsonScalar, VertexAnnotations, VertexAttrs,
    VertexData, VertexName,
};
use crate::graph::{Coord, EType, GraphLike, VData, VType, V};
use crate::phase::Phase;

use std::collections::{BTreeMap, HashMap};

impl JsonGraph {
    /// Encode a graph using the json representation.
    pub fn from_graph(graph: &impl GraphLike) -> Result<Self, JsonError> {
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
            let coord = graph.coord(v);
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
                    annotation: VertexAnnotations {
                        boundary: true,
                        coord: (coord.x, coord.y),
                        input,
                        output,
                        ..Default::default()
                    },
                    ..Default::default()
                };

                wire_vertices.insert(v_name, attrs);
            } else {
                let phase = graph.phase(v);
                // Encode zero-phases as empty strings by default. If the vertex
                // is a Hadamard node, encode "1" as empty strings instead.
                let phase_options = PhaseOptions {
                    ignore_value: Some(match typ == VType::H {
                        true => Phase::one(),
                        false => Phase::zero(),
                    }),
                    ..Default::default()
                };
                let value = JsonPhase::from_phase(phase, phase_options);
                let mut attrs = VertexAttrs {
                    annotation: VertexAnnotations {
                        coord: (coord.x, coord.y),
                        ..Default::default()
                    },
                    data: VertexData {
                        typ,
                        value,
                        ..Default::default()
                    },
                };

                if typ == VType::ZBox {
                    // Data for ZBox vertices is not currently supported, so write the default.
                    attrs.annotation.label = Some("1".to_string());
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
                            annotation: VertexAnnotations {
                                coord: (coord.x, coord.y),
                                ..Default::default()
                            },
                            data: VertexData {
                                typ: VType::H,
                                is_edge: true,
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

        let scalar = graph.scalar();
        let scalar = scalar.is_one().then(|| JsonScalar::from_scalar(scalar));

        Ok(Self {
            wire_vertices,
            node_vertices,
            undir_edges,
            variable_types: Default::default(),
            scalar,
        })
    }

    /// Decode a graph from the json representation.
    pub fn to_graph<G: GraphLike>(&self) -> Result<G, JsonError> {
        let mut graph = G::new();

        if !self.variable_types.is_empty() {
            unimplemented!("Variables are not currently supported.");
        }

        let mut names: HashMap<VertexName, V> = HashMap::new();

        // Map used to track auxiliary Hadamard nodes that should be decoded as Hadamard edges.
        // Stores the neighbor nodes of the Hadamard node, and the coordinate of the Hadamard node.
        let mut hadamards: HashMap<&str, (Vec<V>, Coord)> = HashMap::new();

        for (name, attrs) in &self.node_vertices {
            let coord = Coord {
                x: attrs.annotation.coord.0,
                y: attrs.annotation.coord.1,
            };
            if attrs.data.typ == VType::H && attrs.data.is_edge {
                // A virtual hadamard edge.
                hadamards.insert(name, (vec![], coord));
                continue;
            }

            let phase = attrs
                .data
                .value
                .to_phase()
                .map_err(|_| JsonError::InvalidNodePhase {
                    name: name.to_string(),
                    phase: attrs.data.value.0.clone(),
                })?;
            let phase = match (phase, attrs.data.typ) {
                (Some(r), _) => r,
                // The phase defaults to one for Hadamard nodes,
                (None, VType::H) => Rational64::one().into(),
                // and zero for all others.
                (None, _) => Rational64::zero().into(),
            };
            let v = graph.add_vertex_with_data(VData {
                ty: attrs.data.typ,
                qubit: coord.qubit(),
                row: coord.row(),
                phase,
            });
            names.insert(name.to_string(), v);
        }

        // Insert the boundary nodes, and collect the input and output vectors.
        let mut inputs: BTreeMap<usize, &str> = BTreeMap::new();
        let mut outputs: BTreeMap<usize, &str> = BTreeMap::new();
        for (name, attrs) in &self.wire_vertices {
            let coord = Coord {
                x: attrs.annotation.coord.0,
                y: attrs.annotation.coord.1,
            };
            let v = graph.add_vertex_with_data(VData {
                ty: VType::B,
                qubit: coord.qubit(),
                row: coord.row(),
                phase: Phase::zero(),
            });
            names.insert(name.to_string(), v);
            if let Some(input) = attrs.annotation.input {
                inputs.insert(input, name);
            }
            if let Some(output) = attrs.annotation.output {
                outputs.insert(output, name);
            }
        }
        graph.set_inputs(inputs.into_values().map(|name| names[name]).collect());
        graph.set_outputs(outputs.into_values().map(|name| names[name]).collect());

        // Insert the edges.
        for attrs in self.undir_edges.values() {
            let src = || names[&attrs.src];
            let tgt = || names[&attrs.tgt];

            match (
                hadamards.get(attrs.src.as_str()),
                hadamards.get(attrs.tgt.as_str()),
            ) {
                (Some((_, src_coord)), Some((_, tgt_coord))) => {
                    // Both ends are virtual Hadamard nodes.
                    //
                    // Not sure how this is possible, but pyzx supports it.
                    let new_coord = avg_coord(*src_coord, *tgt_coord);
                    let v = graph.add_vertex_with_data(VData {
                        ty: VType::Z,
                        qubit: new_coord.qubit(),
                        row: new_coord.row(),
                        phase: Phase::zero(),
                    });
                    let name = format!("v{}", graph.num_vertices());
                    names.insert(name, v);
                    hadamards.get_mut(attrs.src.as_str()).unwrap().0.push(v);
                    hadamards.get_mut(attrs.tgt.as_str()).unwrap().0.push(v);
                    continue;
                }
                (Some(_), None) => {
                    hadamards.get_mut(attrs.src.as_str()).unwrap().0.push(tgt());
                    continue;
                }
                (None, Some(_)) => {
                    hadamards.get_mut(attrs.tgt.as_str()).unwrap().0.push(src());
                    continue;
                }
                _ => {}
            }

            graph.add_edge_smart(src(), tgt(), attrs.typ);
        }

        // Add the Hadamard edges.
        for (neighbors, _) in hadamards.values() {
            if neighbors.len() != 2 {
                panic!("Virtual Hadamard node has wrong number of neighbors.");
            }
            let (src, tgt) = (neighbors[0], neighbors[1]);
            graph.add_edge_smart(src, tgt, EType::H);
        }

        // Set the scalar.
        if let Some(scalar) = &self.scalar {
            *graph.scalar_mut() = scalar.to_scalar()?;
        }

        Ok(graph)
    }
}

/// Returns the average of two coordinates, as a pair of f64.
///
/// Rounds the result to 3 decimal places.
fn avg_coord(a: Coord, b: Coord) -> Coord {
    Coord {
        x: ((a.x + b.x) * 1000.).round() / 2000.,
        y: ((a.y + b.y) * 1000.).round() / 2000.,
    }
}
