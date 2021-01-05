use crate::graph::*;
use std::iter::FromIterator;

pub fn check_spider_fusion(g: &impl IsGraph, v0: V, v1: V) -> bool {
    ((g.vertex_type(v0) == VType::Z && g.vertex_type(v1) == VType::Z) ||
     (g.vertex_type(v0) == VType::X && g.vertex_type(v1) == VType::X)) &&
        g.connected(v0,v1) && g.edge_type(v0,v1) == EType::N
}

pub fn spider_fusion_unsafe(g: &mut impl IsGraph, v0: V, v1: V) {
    for (v,et) in Vec::from_iter(g.incident_edges(v1)) {
        if v != v0 {
            g.add_edge_smart(v0, v, et);
        }
    }

    g.add_to_phase(v0, g.phase(v1));
    g.remove_vertex(v1);
}

pub fn spider_fusion(g: &mut impl IsGraph, v0: V, v1: V) -> bool {
    if check_spider_fusion(g, v0, v1) {
        spider_fusion_unsafe(g, v0, v1); true
    } else { false }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Scalar;
    use crate::vec_graph::Graph;
    use num::Rational;

    #[test]
    fn spider_fusion_simple() {
        let mut g = Graph::new();
        let vs = [
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];

        g.set_phase(vs[2], Rational::new(1,2));
        g.set_phase(vs[3], Rational::new(1,4));

        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[2]);
        g.add_edge(vs[2], vs[3]);
        g.add_edge(vs[3], vs[4]);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[3], vs[6]);

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 6);

        spider_fusion(&mut g, vs[2], vs[3]);

        assert_eq!(g.num_vertices(), 6);
        assert_eq!(g.num_edges(), 5);
        assert_eq!(g.degree(vs[2]), 5);

        assert_eq!(g.phase(vs[2]), Rational::new(3,4));
    }

    #[test]
    fn spider_fusion_smart() {
        // a spider fusion that creates parallel edges, which get
        // removed by complementarity
        let mut g = Graph::new();
        let vs = [
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::Z),
            g.add_vertex(VType::X),
            g.add_vertex(VType::B),
            g.add_vertex(VType::B),
        ];

        g.set_phase(vs[2], Rational::new(1,2));
        g.set_phase(vs[3], Rational::new(1,4));

        g.add_edge(vs[0], vs[2]);
        g.add_edge(vs[1], vs[2]);
        g.add_edge(vs[2], vs[3]);
        g.add_edge(vs[2], vs[4]);
        g.add_edge(vs[3], vs[4]);
        g.add_edge(vs[3], vs[5]);
        g.add_edge(vs[4], vs[6]);

        assert_eq!(g.num_vertices(), 7);
        assert_eq!(g.num_edges(), 7);
        assert_eq!(g.degree(vs[2]), 4);
        assert_eq!(g.degree(vs[4]), 3);
        assert_eq!(g.scalar, Scalar::one());

        spider_fusion(&mut g, vs[2], vs[3]);

        assert_eq!(g.num_vertices(), 6);
        assert_eq!(g.num_edges(), 4);
        assert_eq!(g.degree(vs[2]), 3);
        assert_eq!(g.degree(vs[4]), 1);
        assert_eq!(g.scalar, Scalar::rt2_pow(-2));

        assert_eq!(g.phase(vs[2]), Rational::new(3,4));
    }
}
