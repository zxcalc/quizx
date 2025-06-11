use crate::graph::{VData, VType};
use crate::hash_graph::Graph;
use crate::hash_graph::GraphLike;
use crate::phase::Phase;
use num::Rational64;

fn parse_phase(s: &str) -> Phase {
    let rat: Rational64 = s.parse().unwrap(); // panics on malformed input
    Phase::new(rat)
}

/// Loads a quizx graph from a .zxg file (such as from zxlive)
pub fn load_graph(path: &str) -> Graph {
    let file_content = std::fs::read_to_string(path).unwrap();
    let data: serde_json::Value = serde_json::from_str(&file_content).unwrap();

    let wire_vertices = data["wire_vertices"].as_object().unwrap();
    let node_vertices = data["node_vertices"].as_object().unwrap();
    let undir_edges = data["undir_edges"].as_object().unwrap();

    let mut g = Graph::new();
    let mut id_map = std::collections::HashMap::new();

    for (id, obj) in wire_vertices {
        let coord = &obj["annotation"]["coord"];
        let qubit = coord[1].as_f64().unwrap();
        let row = coord[0].as_f64().unwrap();
        let index = g.add_vertex_with_data(VData {
            ty: VType::B,
            qubit,
            // Currently the parity parameter is not handled
            vars: crate::params::Parity::new(*Box::new([]), false),
            row,
            phase: Phase::new(0),
        });
        id_map.insert(id.clone(), index);
    }

    for (id, obj) in node_vertices {
        let coord = &obj["annotation"]["coord"];
        let qubit = coord[1].as_f64().unwrap();
        let row = coord[0].as_f64().unwrap();

        let vtype = match obj["data"]["type"].as_str().unwrap() {
            "Z" => VType::Z,
            "X" => VType::X,
            "H" => VType::H,
            "B" => VType::B,
            _ => panic!("Unknown node type"),
        };

        let phase = obj["data"]["phase"]
            .as_str()
            .map(parse_phase)
            .unwrap_or(Phase::new(0));

        let index = g.add_vertex_with_data(VData {
            ty: vtype,
            qubit,
            // currently parity is not handled
            vars: crate::params::Parity::new(*Box::new([]), false),
            row,
            phase,
        });
        id_map.insert(id.clone(), index);
    }

    for value in undir_edges.values() {
        let a_id = value["src"].as_str().expect("missing 'src'");
        let b_id = value["tgt"].as_str().expect("missing 'tgt'");

        let a = *id_map
            .get(a_id)
            .unwrap_or_else(|| panic!("Unknown vertex ID: {}", a_id));

        let b = *id_map
            .get(b_id)
            .unwrap_or_else(|| panic!("Unknown vertex ID: {}", b_id));

        g.add_edge(a, b);
    }

    g
}

#[cfg(test)]
pub mod tests {
    use super::*;
    fn test_file(name: &str) -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../")
            .join("test_files")
            .join(name)
            .to_str()
            .unwrap()
            .to_string()
    }
    #[test]
    fn test_load_graph() {
        let g = load_graph(&test_file("xx_stab.zxg"));
        println!(
            "Graph loaded with {} vertices and {} edges",
            g.num_vertices(),
            g.num_edges()
        );
    }
}
