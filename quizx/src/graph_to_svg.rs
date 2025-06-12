use crate::detection_webs::{Pauli, PauliWeb};
use crate::graph::GraphLike;
use crate::hash_graph::{Graph, VType};
use num::{FromPrimitive, Rational64};
use std::fmt::Write;

/// Helper function
/// Formats phase to string
fn format_phase(phase: f64) -> Option<String> {
    let tol = 1e-6;

    if phase.abs() < tol {
        return None;
    }

    let rat = Rational64::from_f64(phase).unwrap_or_else(|| {
        Rational64::from_f64((phase * 10.0).round() / 10.0).unwrap_or(Rational64::from_integer(0))
    });

    let numer = *rat.numer();
    let denom = *rat.denom();

    let label = match (numer, denom) {
        (1, 1) => "π".to_string(),
        (-1, 1) => "-π".to_string(),
        (n, 1) => format!("{n}π"),
        (1, d) => format!("π/{d}"),
        (-1, d) => format!("-π/{d}"),
        (n, d) => format!("{n}π/{d}"),
    };

    Some(label)
}

/// Optional overlay edge with custom color and width
pub struct OverlayEdge {
    pub from: usize,
    pub to: usize,
    pub color: &'static str,
    pub opacity: f32,
    pub width: f32,
}

/// Maps Pauli operator to a highlight color
fn pauli_color(pauli: Pauli) -> &'static str {
    match pauli {
        Pauli::X => "red",
        Pauli::Y => "blue",
        Pauli::Z => "green",
    }
}

/// Renders a quizx graph to SVG with optional overlayed translucent colored edges
/// (The translucent colored edges are for pauliwebs)
pub fn graph_to_svg(graph: &Graph, overlays: &[OverlayEdge]) -> String {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let scale = 80.0;

    for v in graph.vertices() {
        let d = graph.vertex_data(v);
        let x = d.row * scale;
        let y = d.qubit * scale;
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
    }

    let padding = 100.0;
    let width = (max_x - min_x) + 2.0 * padding;
    let height = (max_y - min_y) + 2.0 * padding;
    let offset_x = -min_x + padding;
    let offset_y = -min_y + padding;
    let mut svg = String::new();
    let radius = 15.0;

    writeln!(
        &mut svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"#
    ).unwrap();
    writeln!(
        &mut svg,
        r#"<rect width="100%" height="100%" fill="white"/>"#
    )
    .unwrap();

    // Highlighted edges
    for edge in overlays {
        let (x1, y1) = {
            let d = graph.vertex_data(edge.from);
            (d.row * scale + offset_x, d.qubit * scale + offset_y)
        };
        let (x2, y2) = {
            let d = graph.vertex_data(edge.to);
            (d.row * scale + offset_x, d.qubit * scale + offset_y)
        };
        writeln!(
            &mut svg,
            r#"<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-opacity="{opacity}" stroke-width="{width}" />"#,
            color = edge.color,
            opacity = edge.opacity,
            width = edge.width
        ).unwrap();
    }

    // Standard edges
    for (a, b, _) in graph.edges() {
        let pa = graph.vertex_data(a);
        let pb = graph.vertex_data(b);
        let (x1, y1) = (pa.row * scale + offset_x, pa.qubit * scale + offset_y);
        let (x2, y2) = (pb.row * scale + offset_x, pb.qubit * scale + offset_y);
        writeln!(
            &mut svg,
            r#"<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="black" stroke-width="2"/>"#
        )
        .unwrap();
    }

    // Nodes
    for v in graph.vertices() {
        let data = graph.vertex_data(v);
        let (x, y) = (data.row * scale + offset_x, data.qubit * scale + offset_y);
        let color = match graph.vertex_type(v) {
            VType::Z => "green",
            VType::X => "red",
            VType::H => "yellow",
            VType::B => "black",
            _ => "gray",
        };

        writeln!(
            &mut svg,
            r#"<circle cx="{x}" cy="{y}" r="{radius}" fill="{color}" stroke="black"/>"#
        )
        .unwrap();

        // Node ID above the circle
        writeln!(
            &mut svg,
            r#"<text x="{x}" y="{y_offset}" text-anchor="middle" dominant-baseline="middle" fill="black" font-size="12">{}</text>"#,
            v,
            y_offset = y - radius - 8.0
        ).unwrap();

        // Phase inside the circle (only if nonzero)
        if let Some(label) = format_phase(data.phase.to_f64()) {
            writeln!(
                &mut svg,
                r#"<text x="{x}" y="{y}" text-anchor="middle" dominant-baseline="middle" fill="white" font-size="12">{label}</text>"#
            ).unwrap();
        }
    }

    svg.push_str("</svg>");
    svg
}

/// Converts a PauliWeb to overlay edges and renders everything
pub fn graph_to_svg_with_pauliweb(graph: &Graph, pauliweb: Option<&PauliWeb>) -> String {
    let mut overlays = vec![];

    if let Some(web) = pauliweb {
        for ((from, to), pauli) in &web.edge_operators {
            overlays.push(OverlayEdge {
                from: *from,
                to: *to,
                color: pauli_color(*pauli),
                opacity: 0.4,
                width: 9.0,
            });
        }
    }

    graph_to_svg(graph, &overlays)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec_graph::VData;
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
    fn test_graph_to_svg() {
        let mut g = Graph::new();

        let vdata_a = VData {
            ty: VType::Z,
            qubit: 0.0,
            row: 0.0,
            ..Default::default()
        };
        let a = g.add_vertex_with_data(vdata_a);

        let vdata_b = VData {
            ty: VType::X,
            qubit: 1.0,
            row: 1.0,
            phase: crate::phase::Phase::new(Rational64::from_f64(-1.5).unwrap()),
            ..Default::default()
        };
        let b = g.add_vertex_with_data(vdata_b);

        g.add_edge(a, b);

        let mut web = PauliWeb::new();
        web.set_edge(a, b, Pauli::X);

        let svg = graph_to_svg_with_pauliweb(&g, Some(&web));
        let output_path = test_file("test_graph_to_svg.svg");

        std::fs::write(output_path, svg).unwrap();
    }
}
