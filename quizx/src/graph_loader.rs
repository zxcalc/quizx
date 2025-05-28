use crate::graph::{VData, VType};
use crate::hash_graph::Graph;
use crate::hash_graph::GraphLike;
use crate::phase::Phase;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
// The only reason we need fs here is because a graph loader is supposed to act on a file.
use crate::{fscalar::*, params::Parity};
use std::fs;

pub fn load_graph(path: &str) -> Result<Graph, String> {
    // Load as JSON file
    let file_content = match fs::read_to_string(path) {
        Ok(content) => content,
        Err(e) => return Err(format!("Failed to read file: {}", e)),
    };

    let data: Value = match serde_json::from_str(&file_content) {
        Ok(json) => json,
        Err(e) => return Err(format!("Failed to parse JSON: {}", e)),
    };

    // Verify required JSON structure
    let wire_vertices = data["wire_vertices"]
        .as_object()
        .ok_or("Missing or invalid wire_vertices")?;
    let node_vertices = data["node_vertices"]
        .as_object()
        .ok_or("Missing or invalid node_vertices")?;
    let _undir_edges = data["undir_edges"]
        .as_object()
        .ok_or("Missing or invalid undir_edges")?;

    let mut xcods: HashSet<i64> = HashSet::new();
    let mut ycods: HashSet<i64> = HashSet::new();

    // Collect coordinates from wire vertices
    for (_node, dets) in wire_vertices {
        let coord = match dets["annotation"].get("coord") {
            Some(coord) => coord.as_array().ok_or("Invalid coordinate format")?,
            None => {
                // Handle boundary vertices with boundary field
                let boundary = dets["annotation"]["boundary"]
                    .as_bool()
                    .ok_or("Invalid boundary field")?;
                if !boundary {
                    return Err("Invalid boundary vertex format".to_string());
                }
                continue;
            }
        };
        let x = (coord[0]
            .as_f64()
            .ok_or("Invalid x coordinate (not a number)")?
            * 1000.0) as i64;
        let y = (coord[1]
            .as_f64()
            .ok_or("Invalid y coordinate (not a number)")?
            * 1000.0) as i64;
        xcods.insert(x);
        ycods.insert(y);
    }

    // Collect coordinates from node vertices
    for (_node, dets) in node_vertices {
        let coord = dets["annotation"]["coord"]
            .as_array()
            .ok_or("Invalid coordinate format")?;
        let x = (coord[0]
            .as_f64()
            .ok_or("Invalid x coordinate (not a number)")?
            * 1000.0) as i64;
        let y = (coord[1]
            .as_f64()
            .ok_or("Invalid y coordinate (not a number)")?
            * 1000.0) as i64;
        xcods.insert(x);
        ycods.insert(y);
    }

    let mut graph = Graph::new();
    let mut id_map = HashMap::new();

    // Collect coordinates from wire vertices
    for (_node, dets) in wire_vertices {
        let coord = dets["annotation"]["coord"]
            .as_array()
            .ok_or("Invalid coordinate format")?;
        let x = (coord[0]
            .as_f64()
            .ok_or("Invalid x coordinate (not a number)")?
            * 1000.0) as i64;
        let y = (coord[1]
            .as_f64()
            .ok_or("Invalid y coordinate (not a number)")?
            * 1000.0) as i64;
        xcods.insert(x);
        ycods.insert(y);
    }

    // Collect coordinates from node vertices
    for (_node, dets) in node_vertices {
        let coord = dets["annotation"]["coord"]
            .as_array()
            .ok_or("Invalid coordinate format")?;
        let x = (coord[0]
            .as_f64()
            .ok_or("Invalid x coordinate (not a number)")?
            * 1000.0) as i64;
        let y = (coord[1]
            .as_f64()
            .ok_or("Invalid y coordinate (not a number)")?
            * 1000.0) as i64;
        xcods.insert(x);
        ycods.insert(y);
    }

    let mut x_list: Vec<_> = xcods.iter().cloned().collect();
    let mut y_list: Vec<_> = ycods.iter().cloned().collect();
    x_list.sort();
    y_list.sort();

    let x_cood_map: HashMap<i64, usize> = x_list.iter().enumerate().map(|(n, &x)| (x, n)).collect();
    let y_cood_map: HashMap<i64, usize> = y_list.iter().enumerate().map(|(n, &y)| (y, n)).collect();

    let x_cood_map_f64: HashMap<i64, f64> = x_list
        .iter()
        .enumerate()
        .map(|(_n, &x)| (x, x as f64 / 1000.0))
        .collect();
    let y_cood_map_f64: HashMap<i64, f64> = y_list
        .iter()
        .enumerate()
        .map(|(_n, &y)| (y, y as f64 / 1000.0))
        .collect();

    // Boundary vertices
    for (node, dets) in data["wire_vertices"].as_object().unwrap() {
        let coord = dets["annotation"]["coord"].as_array().unwrap();
        let row = coord[0].as_f64().unwrap();
        let qubit = coord[1].as_f64().unwrap();
        let v_val = dets["data"]["value"].as_f64().unwrap_or(0.0);
        let data: VData = VData {
            ty: VType::B,
            vars: Parity::zero(),
            phase: Phase::from_f64(v_val),
            qubit: qubit,
            row: row,
        };
        let vid = graph.add_vertex_with_data(data);
        id_map.insert(node.clone(), vid);
    }

    // Actual vertices
    for (node, dets) in data["node_vertices"].as_object().unwrap() {
        let coord = dets["annotation"]["coord"].as_array().unwrap();
        let x = (coord[0].as_f64().unwrap() * 1000.0) as i64;
        let y = (coord[1].as_f64().unwrap() * 1000.0) as i64;
        let _row = x_cood_map[&x];
        let _qubit = y_cood_map[&y];
        let v_val = dets["data"]["value"].as_f64().unwrap_or(0.0);
        let v_type = match dets["data"]["type"].as_str().unwrap() {
            "X" => VType::X,
            "Z" => VType::Z,
            _ => VType::H,
        };
        let data: VData = VData {
            ty: v_type,
            vars: Parity::zero(),
            phase: Phase::from_f64(v_val),
            qubit: y_cood_map_f64[&y],
            row: x_cood_map_f64[&x],
        };
        let vid = graph.add_vertex_with_data(data);
        id_map.insert(node.clone(), vid);
    }

    // Edges
    for (_edge, dets) in data["undir_edges"].as_object().unwrap() {
        let src = dets["src"].as_str().unwrap();
        let tgt = dets["tgt"].as_str().unwrap();
        let src_id = id_map[src];
        let tgt_id = id_map[tgt];
        graph.add_edge(src_id, tgt_id); //, ety); for now lets just do simple edges
    }

    Ok(graph)
}

// Tests
#[cfg(test)]
mod tests {

    use super::*;
    use crate::graph::GraphLike;
    use crate::graph::VType;
    use crate::phase::Phase;
    use std::collections::HashSet;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_load_graph_vertices() {
        // Create a temporary test JSON file
        let test_json = r#"
        {
            "wire_vertices": {
                "w1": {
                    "annotation": { "coord": [0, 0] }
                },
                "w2": {
                    "annotation": { "coord": [1, 0] }
                }
            },
            "node_vertices": {
                "n1": {
                    "annotation": { "coord": [0, 1] },
                    "data": {
                        "type": "X",
                        "value": 0.5
                    }
                },
                "n2": {
                    "annotation": { "coord": [1, 1] },
                    "data": {
                        "type": "Z",
                        "value": 1.0
                    }
                }
            },
            "undir_edges": {
                "e1": {
                    "src": "w1",
                    "tgt": "n1"
                },
                "e2": {
                    "src": "n1",
                    "tgt": "n2"
                },
                "e3": {
                    "src": "n2",
                    "tgt": "w2"
                }
            }
        }"#;

        // Create temporary file
        let temp_dir = tempdir().unwrap();
        let temp_file = temp_dir.path().join("test_graph.json");
        std::fs::write(&temp_file, test_json).unwrap();

        // Load the graph
        let graph = load_graph(temp_file.to_str().unwrap()).unwrap();

        // Verify vertices
        assert_eq!(graph.num_vertices(), 4);

        // Verify vertex types and phases
        for v in graph.vertices() {
            let data = graph.vertex_data(v);
            match data.ty {
                VType::X => assert_eq!(data.phase, Phase::from_f64(0.5)),
                VType::Z => assert_eq!(data.phase, Phase::from_f64(1.0)),
                _ => assert_eq!(data.phase, Phase::from_f64(0.0)),
            }
        }

        // Verify edges
        assert_eq!(graph.num_edges(), 3);
    }

    #[test]
    fn test_load_graph_coordinates() {
        let test_json = r#"
        {
            "wire_vertices": {
                "w1": {
                    "annotation": { "coord": [0, 0] }
                },
                "w2": {
                    "annotation": { "coord": [2, 0] }
                }
            },
            "node_vertices": {
                "n1": {
                    "annotation": { "coord": [1, 1] },
                    "data": {
                        "type": "X",
                        "value": 0.0
                    }
                }
            },
            "undir_edges": {
                "e1": {
                    "src": "w1",
                    "tgt": "n1"
                },
                "e2": {
                    "src": "n1",
                    "tgt": "w2"
                }
            }
        }"#;

        let temp_dir = tempfile::tempdir().unwrap();
        let temp_file = temp_dir.path().join("test_graph.json");
        fs::write(&temp_file, test_json).unwrap();

        let graph = load_graph(temp_file.to_str().unwrap());

        match graph {
            Ok(graph) => {
                // Verify coordinates are properly mapped
                let mut coords = HashSet::new();
                for v in graph.vertices() {
                    let data = graph.vertex_data(v);
                    coords.insert((data.row as i64, data.qubit as i64));
                }

                // Should have coordinates (0,0), (2,0), (1,1)
                assert_eq!(coords.len(), 3);
                assert!(coords.contains(&(0, 0)));
                assert!(coords.contains(&(2, 0)));
                assert!(coords.contains(&(1, 1)));
            }
            Err(_) => panic!("Failed to load graph"),
        }
    }

    #[test]
    fn test_load_graph_edge_types() {
        let test_json = r#"
        {
            "wire_vertices": {
                "w1": {
                    "annotation": { "coord": [0, 0] }
                },
                "w2": {
                    "annotation": { "coord": [1, 0] }
                }
            },
            "node_vertices": {
                "n1": {
                    "annotation": { "coord": [0, 1] },
                    "data": {
                        "type": "X",
                        "value": 0.0
                    }
                }
            },
            "undir_edges": {
                "e1": {
                    "src": "w1",
                    "tgt": "n1"
                },
                "e2": {
                    "src": "n1",
                    "tgt": "w2"
                }
            }
        }"#;

        let temp_dir = tempfile::tempdir().unwrap();
        let temp_file = temp_dir.path().join("test_graph.json");
        std::fs::write(&temp_file, test_json).unwrap();

        let graph = load_graph(temp_file.to_str().unwrap()).unwrap();

        // Verify edges
        assert_eq!(graph.num_edges(), 2);

        // Verify connectivity
        let mut edges = Vec::new();
        for (src, tgt, _) in graph.edge_vec() {
            edges.push((src, tgt));
        }

        assert_eq!(edges.len(), 2);
        // We can't predict exact vertex IDs, but we can verify the connectivity pattern
        assert!(edges.iter().any(|&(src, tgt)| src != tgt)); // Should have at least one edge between different vertices
        assert!(edges.iter().all(|&(src, tgt)| src != tgt)); // No self-loops
    }

    #[test]
    #[should_panic(expected = "Missing or invalid node_vertices")]
    fn test_load_graph_invalid_json() {
        // Test loading with invalid JSON
        let invalid_json = r#"{
            "wire_vertices": {
                "w1": {
                    "annotation": { "coord": [0, 0] }
                }
            }
        }"#;
        let temp_dir = tempfile::tempdir().unwrap();
        let temp_file = temp_dir.path().join("invalid.json");
        std::fs::write(&temp_file, invalid_json).unwrap();

        load_graph(temp_file.to_str().unwrap()).unwrap();
    }

    #[test]
    fn test_from_file() {
        use super::*;

        // Get the project root directory (where Cargo.toml is located)
        let manifest_dir =
            std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR should be set by Cargo");

        // Navigate to the project root and then to test_files
        let path = std::path::Path::new(&manifest_dir)
            .parent() // Go up from quizx/src to quizx
            .expect("Should have parent directory")
            .join("test_files")
            .join("xxx_final.zxg");

        // Convert to string and verify the file exists
        let path_str = path.to_str().expect("Path should be valid UTF-8");

        if !path.exists() {
            panic!("Test file not found at: {}", path.display());
        }

        // Load the graph
        let graph = load_graph(path_str).expect("Failed to load graph from file");

        // Basic validation
        assert!(
            graph.num_vertices() > 0,
            "Loaded graph should have vertices"
        );
    }
}
