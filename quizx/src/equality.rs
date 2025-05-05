use approx::abs_diff_eq;
use num::Complex;

use crate::circuit::Circuit;
use crate::graph::GraphLike;
use crate::simplify::full_simp;
use crate::tensor::ToTensor;
use crate::vec_graph::Graph;

/// Checks the equality of two circuits by comparing the linear maps they represent.
/// This approach is only feasible for a small number of qubits.
pub fn compare_tensors(c1: &Circuit, c2: &Circuit) -> bool {
    c1.to_tensorf() == c2.to_tensorf()
}

/// Implements `Circuit.verify_equality` from pyzx.
///
/// Verifies the equality of two circuits by investigating whether they "cancel each other out".
/// This is done by composing one circuit with the adjoint of the other. If simplifying the
/// result yields the identity, then the two circuits are verifiably equal.
///
/// Note that the simplification may not yield the identity even if both circuits are equal,
/// which is why this approch gives an inconclusive answer if the resulting circuit is not
/// the identity. In general, this approach can't verify that two circuits are unequal.
///
/// This approach is feasible for a high number of qubits.
pub fn verify_equality_with_options(
    c1: &Circuit,
    c2: &Circuit,
    up_to_global_phase: bool,
) -> Option<bool> {
    if c1.num_qubits() != c2.num_qubits() {
        // both circuits are verifiably unequal due to an unequal number of qubits
        return Some(false);
    }
    let c = c1.to_adjoint() + c2;
    let mut g: Graph = c.to_graph();
    full_simp(&mut g);
    if g.is_identity() {
        if !up_to_global_phase {
            // both circuits are verifiably equal if the resulting global phase is zero
            // otherwise, they are verifiably unequal / only equal up to global phase
            let c: Complex<f64> = g.scalar().into();
            return Some(abs_diff_eq!(c.arg(), 0.0));
        }
        // both circuits are verifiably equal up to global phase
        return Some(true);
    }
    // both circuits are neither verifiably equal nor verifiably unequal
    None
}

/// Verifies the equality of two circuits up to global phase.
pub fn verify_equality(c1: &Circuit, c2: &Circuit) -> Option<bool> {
    verify_equality_with_options(c1, c2, true)
}

#[cfg(test)]
mod tests {
    use num::Rational64;

    use super::compare_tensors;
    use super::verify_equality_with_options;
    use crate::circuit::Circuit;

    /// Inspired by `BothCircuitsEmptyZXChecker` found in `test_equality.cpp` from mqt-qcec
    #[test]
    fn both_circuits_empty() {
        let c1 = Circuit::new(1);
        let c2 = Circuit::new(1);

        // c1 and c2 are equal
        assert!(compare_tensors(&c1, &c2));
        assert!(verify_equality_with_options(&c1, &c2, false).unwrap());
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
        assert!(verify_equality_with_options(&c1, &c2, true).unwrap());
        // c1 and c2 are unequal
        assert!(!compare_tensors(&c1, &c2));
        assert!(!verify_equality_with_options(&c1, &c2, false).unwrap());
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
        assert!(!compare_tensors(&c1, &c2));
        // c1 and c2 are not verifiably unequal with a simplification-based approach
        assert!(verify_equality_with_options(&c1, &c2, true).is_none());
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
        assert!(verify_equality_with_options(&c1, &c2, false).unwrap());
    }
}
