use quizx::circuit::Circuit;
use quizx::equality::equal_circuit_tensor;
use quizx::gate::{GType, Gate};

/// Regression test for https://github.com/zxcalc/quizx/issues/137
#[test]
fn cccx_decompositions_with_and_without_ancilla_agree() {
    let noancilla = Circuit::from_qasm(include_str!("fixtures/cccx_noancilla.qasm")).unwrap();
    let mut with_ancilla = Circuit::from_qasm(include_str!("fixtures/cccx_v_chain.qasm")).unwrap();
    with_ancilla.push_front(Gate::new(GType::InitAncilla, vec![4]));
    with_ancilla.push_back(Gate::new(GType::PostSelect, vec![4]));
    assert!(equal_circuit_tensor(&noancilla, &with_ancilla));
}
