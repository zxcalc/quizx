use crate::graph::*;
use crate::basic_rules::*;

pub fn vertex_simp<G: GraphLike>(
    g: &mut G,
    check: fn(&G, V) -> bool,
    rule: fn(&mut G, V) -> ()
    ) -> bool
{
    let mut got_match = false;
    loop {
        match g.find_vertex(|v| check(g, v)) {
            Some(v) => rule(g, v),
            None => break,
        };
        got_match = true;
    }

    got_match
}

pub fn edge_simp<G: GraphLike>(
    g: &mut G,
    check: fn(&G, V, V) -> bool,
    rule: fn(&mut G, V, V) -> ()
    ) -> bool
{
    let mut got_match = false;
    loop {
        match g.find_edge(|v0, v1, _| check(g, v0, v1)) {
            Some((v0, v1, _)) => rule(g, v0, v1),
            None => break,
        };
        got_match = true;
    }

    got_match
}

pub fn id_simp(g: &mut impl GraphLike) -> bool {
    vertex_simp(g, check_remove_id, remove_id_unchecked)
}

pub fn local_comp_simp(g: &mut impl GraphLike) -> bool {
    vertex_simp(g, check_local_comp, local_comp_unchecked)
}

pub fn spider_simp(g: &mut impl GraphLike) -> bool {
    edge_simp(g, check_spider_fusion, spider_fusion_unchecked)
}

pub fn pivot_simp(g: &mut impl GraphLike) -> bool {
    edge_simp(g, check_pivot, pivot_unchecked)
}

pub fn gen_pivot_simp(g: &mut impl GraphLike) -> bool {
    edge_simp(g, check_gen_pivot_reduce, gen_pivot_unchecked)
}

pub fn interior_clifford_simp(g: &mut impl GraphLike) -> bool {
    spider_simp(g);
    g.x_to_z();
    let mut got_match = false;
    let mut m = true;
    while m {
        m = false;
        m = m || id_simp(g);
        m = m || spider_simp(g);
        m = m || pivot_simp(g);
        m = m || local_comp_simp(g);
        if m { got_match = true; }
    }

    got_match
}

pub fn clifford_simp(g: &mut impl GraphLike) -> bool {
    let mut got_match = false;
    let mut m = true;
    while m {
        m = false;
        m = m || interior_clifford_simp(g);
        m = m || gen_pivot_simp(g);
        if m { got_match = true; }
    }

    got_match
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::*;
    use crate::vec_graph::Graph;
    use crate::tensor::ToTensor;

    #[test]
    fn cliff_simp() {
        let c = Circuit::from_qasm(r#"
            qreg q[3];
            h q[1];
            h q[1];
            s q[0];
            cx q[0], q[1];
            cx q[1], q[0];
            cx q[0], q[1];
            s q[1];
            cx q[1], q[2];
            cx q[0], q[2];
            h q[2];
            z q[0];
        "#).unwrap();

        let g: Graph = c.to_graph();

        let mut h = g.clone();
        assert!(interior_clifford_simp(&mut h));
        assert_eq!(g.to_tensor4(), h.to_tensor4());

        let mut h = g.clone();
        assert!(clifford_simp(&mut h));
        println!("{}", g.to_dot());
        println!("{}", h.to_dot());
        assert_eq!(g.to_tensor4(), h.to_tensor4());
    }
}
