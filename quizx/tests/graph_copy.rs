use quizx::graph::{EType, GraphLike, VType};
use quizx::params::Expr;
use quizx::scalar::{FromPhase, Scalar4};
use quizx::tensor::ToTensor;

fn copy_preserves_map<G: GraphLike>() {
    let mut g = G::new();
    let hole = g.add_vertex(VType::Z);
    let input = g.add_vertex(VType::B);
    let spider = g.add_vertex_with_phase(VType::Z, (1, 4));
    let output = g.add_vertex(VType::B);
    g.remove_vertex(hole);
    g.add_edge(input, spider);
    g.add_edge_with_type(spider, output, EType::H);
    g.set_inputs(vec![input]);
    g.set_outputs(vec![output]);
    *g.scalar_mut() = Scalar4::from_phase((1, 4));
    let copy = g.copy(false);
    assert_eq!(copy.vertices().max(), Some(copy.num_vertices() - 1));
    assert_eq!(copy.inputs().len(), 1);
    assert_eq!(copy.outputs().len(), 1);
    assert_eq!(copy.to_tensor4(), g.to_tensor4());
    assert_eq!(g.copy(true).to_tensor4(), g.to_adjoint().to_tensor4());

    let other_input = g.add_vertex(VType::B);
    let other_output = g.add_vertex(VType::B);
    g.add_edge(other_input, other_output);
    g.set_inputs(vec![other_input, input]);
    g.set_outputs(vec![output, other_output]);
    assert_eq!(g.copy(false).to_tensor4(), g.to_tensor4());

    let expr = Expr::linear(vec![0]);
    g.mul_scalar_factor(expr.clone(), Scalar4::from_phase((1, 2)));
    let copy = g.copy(false);
    assert_eq!(copy.get_scalar_factor(&expr), g.get_scalar_factor(&expr));
    assert_eq!(copy.scalar_factors().count(), g.scalar_factors().count());
}

#[test]
fn vec_graph_copy() {
    copy_preserves_map::<quizx::vec_graph::Graph>();
}

#[test]
fn hash_graph_copy() {
    copy_preserves_map::<quizx::hash_graph::Graph>();
}
