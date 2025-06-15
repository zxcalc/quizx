use approx::abs_diff_eq;
use num::Complex;

use crate::circuit::Circuit;
use crate::graph::GraphLike;
use crate::simplify::full_simp;
use crate::tensor::ToTensor;
use crate::vec_graph::Graph;

/// Checks if two graphs have the same number of input qubits and output qubits respectively.
/// This check is computationally inexpensive and suitable for any number of qubits.
pub fn equal_graph_dim(g1: &Graph, g2: &Graph) -> bool {
    if g1.inputs().len() != g2.inputs().len() {
        // both graphs have an unequal number of input qubits
        return false;
    }
    if g1.outputs().len() != g2.outputs().len() {
        // both graphs have an unequal number of output qubits
        return false;
    }
    true
}

/// Checks if two circuits have the same number of input qubits and output qubits respectively.
/// This check is computationally inexpensive and suitable for any number of qubits.
pub fn equal_circuit_dim(c1: &Circuit, c2: &Circuit) -> bool {
    let g1: Graph = c1.to_graph();
    let g2: Graph = c2.to_graph();
    equal_graph_dim(&g1, &g2)
}

/// Checks the equality of two circuit graphs by comparing the linear maps they represent.
/// This approach is only feasible for a small number of qubits. (up to 7)
pub fn equal_graph_tensor(g1: &Graph, g2: &Graph) -> bool {
    // First, quickly check if both tensors have the same dimension
    if !equal_graph_dim(g1, g2) {
        return false;
    }
    g1.to_tensor4() == g2.to_tensor4()
}

/// Checks the equality of two circuits by comparing the linear maps they represent.
/// This approach is only feasible for a small number of qubits. (up to 7)
pub fn equal_circuit_tensor(c1: &Circuit, c2: &Circuit) -> bool {
    let g1: Graph = c1.to_graph();
    let g2: Graph = c2.to_graph();
    equal_graph_tensor(&g1, &g2)
}

/// Implements `Circuit.verify_equality` from pyzx.
///
/// Verifies the equality of two circuit graphs by investigating whether they "cancel each other out".
/// This is done by composing one circuit graph with the adjoint of the other. If simplifying the
/// result yields the identity, then the two circuit graphs are verifiably equal.
///
/// Note that the simplification may not yield the identity even if both circuit graphs are equal,
/// which is why this approch gives an inconclusive answer if the resulting circuit is not
/// the identity. In general, this approach can't verify that two circuits are unequal.
///
/// This approach is feasible even for a high number of qubits.
pub fn equal_graph_with_options(g1: &Graph, g2: &Graph, up_to_global_phase: bool) -> Option<bool> {
    if !equal_graph_dim(g1, g2) {
        // both graphs are verifiably unequal due to an unequal number of input qubits or output qubits
        return Some(false);
    }
    let mut g = g1.to_adjoint();
    g.plug(g2);
    full_simp(&mut g);
    if g.is_identity() {
        if !up_to_global_phase {
            // both graphs are verifiably equal if the resulting global phase is zero
            // otherwise, they are verifiably unequal / only equal up to global phase
            let c: Complex<f64> = g.scalar().complex_value();
            return Some(abs_diff_eq!(c.arg(), 0.0));
        }
        // both graphs are verifiably equal up to global phase
        return Some(true);
    }
    // both graphs are neither verifiably equal nor verifiably unequal
    None
}

/// Verifies the equality of two graphs up to global phase.
pub fn equal_graph(g1: &Graph, g2: &Graph) -> Option<bool> {
    equal_graph_with_options(g1, g2, true)
}

/// Verifies the equality of two circuits.
pub fn equal_circuit_with_options(
    c1: &Circuit,
    c2: &Circuit,
    up_to_global_phase: bool,
) -> Option<bool> {
    let g1: Graph = c1.to_graph();
    let g2: Graph = c2.to_graph();
    equal_graph_with_options(&g1, &g2, up_to_global_phase)
}

/// Verifies the equality of two circuits up to global phase.
pub fn equal_circuit(c1: &Circuit, c2: &Circuit) -> Option<bool> {
    equal_circuit_with_options(c1, c2, true)
}

#[cfg(test)]
mod tests {
    use num::Rational64;

    use super::equal_circuit_tensor;
    use super::equal_circuit_with_options;
    use crate::circuit::Circuit;

    /// Inspired by `BothCircuitsEmptyZXChecker` found in `test_equality.cpp` from mqt-qcec
    #[test]
    fn both_circuits_empty() {
        let c1 = Circuit::new(1);
        let c2 = Circuit::new(1);

        // c1 and c2 are equal
        assert!(equal_circuit_tensor(&c1, &c2));
        assert!(equal_circuit_with_options(&c1, &c2, false).unwrap());
    }

    /// Inspired by `CloseButNotEqualConstruction` found in `test_equality.cpp` from mqt-qcec
    #[test]
    fn close_but_not_equal() {
        let mut c1 = Circuit::new(1);
        let mut c2 = Circuit::new(1);

        c1.add_gate("x", vec![0]);
        c2.add_gate("x", vec![0]);

        // slightly change c2
        c2.add_gate_with_phase("rz", vec![0], Rational64::new(1, 1024));

        // c1 and c2 are unequal
        assert!(!equal_circuit_tensor(&c1, &c2));
        // c1 and c2 are not verifiably unequal with a simplification-based approach
        assert!(equal_circuit_with_options(&c1, &c2, true).is_none());
    }

    #[test]
    fn cx_with_ancilla_as_x() {
        let mut c1 = Circuit::new(1);
        c1.add_gate("x", vec![0]);

        let mut c2 = Circuit::new(2);
        c2.add_gate("init_anc", vec![1]); // initiualize ancilla: |0⟩
        c2.add_gate("x", vec![1]); // flip ancilla: |1⟩
        c2.add_gate("cx", vec![1, 0]); // CX controlling for ancilla now behaves like X
        c2.add_gate("x", vec![1]); // flip ancilla back: |0⟩ (otherwise the tensor zeros out!)
        c2.add_gate("post_sel", vec![1]); // remove ancilla

        assert!(equal_circuit_tensor(&c1, &c2));
        assert!(equal_circuit_with_options(&c1, &c2, true).unwrap());
    }

    /// Inspired by `GlobalPhase` found in `test_equality.cpp` from mqt-qcec
    #[test]
    fn global_phase() {
        let mut c1 = Circuit::new(1);
        let mut c2 = Circuit::new(1);

        c1.add_gate("x", vec![0]);
        c2.add_gate("x", vec![0]);

        // flip the global phase of c2
        c2.add_gate("z", vec![0]);
        c2.add_gate("x", vec![0]);
        c2.add_gate("z", vec![0]);
        c2.add_gate("x", vec![0]);

        // c1 and c2 are equal up to global phase
        assert!(equal_circuit_with_options(&c1, &c2, true).unwrap());
        // c1 and c2 are unequal
        assert!(!equal_circuit_tensor(&c1, &c2));
        assert!(!equal_circuit_with_options(&c1, &c2, false).unwrap());
    }

    /// Inspired by `SimulationMoreThan64Qubits` found in `test_equality.cpp` from mqt-qcec
    #[test]
    fn more_than_64_qubits() {
        let mut c1 = Circuit::new(65);
        c1.add_gate("h", vec![0]);
        for i in 1..65 {
            c1.add_gate("cx", vec![0, i]);
        }
        let c2 = c1.clone();

        // comparing the tensors of c1 and c2 is infeasible due to high number of qubits

        // c1 and c2 are verifiably equal
        assert!(equal_circuit_with_options(&c1, &c2, false).unwrap());
    }
}
