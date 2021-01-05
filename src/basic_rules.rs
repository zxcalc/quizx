use crate::graph::*;
use std::iter::FromIterator;
use num::Rational;

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

pub fn check_local_comp(g: &impl IsGraph, v: V) -> bool {
    g.vertex_type(v) == VType::Z &&
    *g.phase(v).denom() == 2 &&
    g.incident_edges(v).all(|(v0,et)|
        g.vertex_type(v0) == VType::Z && et == EType::H)
}

pub fn local_comp_unsafe(g: &mut impl IsGraph, v: V) {
    let p = g.phase(v);

    // add a totally connected graph of the nhd of v
    let ns: Vec<V> = g.neighbors(v).collect();
    for &n0 in &ns {
        g.add_to_phase(n0, -p);
        for &n1 in &ns {
            if n0 < n1 { // avoids double-counting
                g.add_edge_smart(n0, n1, EType::H);
            }
        }
    }
    g.remove_vertex(v);

    let x = ns.len() as i32;
    g.scalar().mul_rt2_pow(((x-1)*(x-2))/2);
    g.scalar().mul_phase(Rational::new(*p.numer(), 4));
}

pub fn local_comp(g: &mut impl IsGraph, v: V) -> bool {
    if check_local_comp(g, v) {
        local_comp_unsafe(g, v); true
    } else { false }
}

pub fn check_pivot(g: &impl IsGraph, v0: V, v1: V) -> bool {
    g.vertex_type(v0) == VType::Z &&
    g.vertex_type(v1) == VType::Z &&
    g.phase(v0).is_integer() &&
    g.phase(v1).is_integer() &&
    g.incident_edges(v0).all(|(w,et)|
        g.vertex_type(w) == VType::Z && et == EType::H)
}

pub fn pivot_unsafe(g: &mut impl IsGraph, v0: V, v1: V) {
    let p0 = g.phase(v0);
    let p1 = g.phase(v1);

    // add a complete bipartite graph between the neighbors of v0
    // and the neighbors of v1
    let ns0: Vec<V> = g.neighbors(v0).collect();
    let ns1: Vec<V> = g.neighbors(v1).collect();
    let mut z: i32 = 0; // the number of neighbors of v0 and v1
    for &n0 in &ns0 {
        g.add_to_phase(v0, p1);
        for &n1 in &ns1 {
            if n0 != v1 && n1 != v0 {
                if n0 != n1 {
                    g.add_edge_smart(n0, n1, EType::H);
                } else {
                    g.add_to_phase(n0, Rational::new(1, 1));
                    z += 1;
                }
            }
        }
    }

    for &n1 in &ns1 {
        g.add_to_phase(n1, p0);
    }

    g.remove_vertex(v0);
    g.remove_vertex(v1);

    let x = (ns0.len() as i32) - z; // the number of neighbors of just v0
    let y = (ns1.len() as i32) - z; // the number of neighbors of just v1
    g.scalar().mul_rt2_pow(x*y + y*z + x*z);

    if *p0.numer() != 0 && *p1.numer() != 0 {
        g.scalar().mul_phase(Rational::new(1,1));
    }
}

pub fn pivot(g: &mut impl IsGraph, v0: V, v1: V) -> bool {
    if check_pivot(g, v0, v1) {
        pivot_unsafe(g, v0, v1); true
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

        let success = spider_fusion(&mut g, vs[2], vs[3]);
        assert!(success, "Spider fusion should match");

        assert_eq!(g.num_vertices(), 6);
        assert_eq!(g.num_edges(), 4);
        assert_eq!(g.degree(vs[2]), 3);
        assert_eq!(g.degree(vs[4]), 1);
        assert_eq!(g.scalar, Scalar::rt2_pow(-2));

        assert_eq!(g.phase(vs[2]), Rational::new(3,4));
    }

    #[test]
    fn local_comp_1() {
        let mut g = Graph::new();
        g.add_vertex(VType::Z);
        g.set_phase(0, Rational::new(1,2));
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_edge_with_type(0,1,EType::H);
        g.add_edge_with_type(0,2,EType::H);
        g.add_edge_with_type(0,3,EType::H);
        g.add_edge_with_type(0,4,EType::H);

        g.add_vertex(VType::B);
        g.add_edge(1,5);

        assert_eq!(g.num_vertices(), 6);
        assert_eq!(g.num_edges(), 5);

        let success = local_comp(&mut g, 0);
        assert!(success, "Local comp should match");

        assert_eq!(g.num_vertices(), 5);
        assert_eq!(g.num_edges(), 7);

        for i in 1..5 {
            assert_eq!(g.phase(i), Rational::new(-1,2));
        }

        assert_eq!(*g.scalar(),
            Scalar::rt2_pow((4-1)*(4-2)/2) *
            Scalar::phase(Rational::new(1,4)));

        let h = g.clone();
        let fail = local_comp(&mut g, 1);
        assert!(!fail, "Local comp should not match");
        assert_eq!(g,h);
    }
}
