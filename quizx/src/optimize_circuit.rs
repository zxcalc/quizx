//! Basic gate-level circuit optimisation (issue #157).
//!
//! Fills the previously-stubbed [`Circuit::optimize`] with a lightweight
//! forward-scanning pass that cancels inverse gate pairs and fuses adjacent
//! phase gates, sliding gates past one another wherever they commute.
//! Measurement, reset, ancilla-init and post-selection gates act as hard
//! barriers. This is the "middle ground" the issue asks for: it does the
//! obvious cancellations and commutations and leaves heavier rewriting to the
//! ZX level.

use crate::circuit::Circuit;
use crate::gate::GType::*;
use crate::gate::{GType, Gate};
use crate::params::Parity;
use num::Zero;

impl Circuit {
    /// Implements (roughly) `basic_optimization` from pyzx.
    ///
    /// Repeatedly scans the gate list, combining each gate with the first
    /// later gate it can cancel or fuse with, provided every gate in between
    /// commutes with it. Runs to a fixpoint.
    pub fn optimize(&mut self) {
        let mut gates: Vec<Gate> = self.gates.iter().cloned().collect();
        loop {
            let mut changed = false;
            let mut i = 0;
            while i < gates.len() {
                let mut j = i + 1;
                let mut combined = false;
                while j < gates.len() {
                    if let Some(res) = try_combine(&gates[i], &gates[j]) {
                        gates.remove(j);
                        match res {
                            None => {
                                gates.remove(i);
                            }
                            Some(g) => {
                                gates[i] = g;
                            }
                        }
                        changed = true;
                        combined = true;
                        break;
                    } else if gates_commute(&gates[i], &gates[j]) {
                        j += 1;
                    } else {
                        break;
                    }
                }
                if !combined {
                    i += 1;
                }
            }
            if !changed {
                break;
            }
        }
        self.gates = gates.into();
    }
}

/// Gates past which nothing may be reordered.
fn is_barrier(t: GType) -> bool {
    matches!(
        t,
        Measure | MeasureReset | InitAncilla | PostSelect | UnknownGate
    )
}

/// Gates diagonal in the computational basis. Any two of these commute, even
/// when they share qubits.
fn is_diagonal(t: GType) -> bool {
    matches!(t, ZPhase | Z | S | T | Sdg | Tdg | CZ)
}

fn shares_qubit(a: &Gate, b: &Gate) -> bool {
    a.qs.iter().any(|q| b.qs.contains(q))
}

/// Conservative commutation test: gates on disjoint qubits always commute;
/// gates sharing qubits commute only when both are diagonal. Barriers never
/// commute. Being conservative here keeps the pass provably correct: it never
/// reorders gates that do not genuinely commute.
fn gates_commute(a: &Gate, b: &Gate) -> bool {
    if is_barrier(a.t) || is_barrier(b.t) {
        return false;
    }
    if !shares_qubit(a, b) {
        return true;
    }
    is_diagonal(a.t) && is_diagonal(b.t)
}

/// Try to combine two gates acting on the same qubits.
/// `None` => not combinable; `Some(None)` => they annihilate;
/// `Some(Some(g))` => they fuse into `g`.
fn try_combine(a: &Gate, b: &Gate) -> Option<Option<Gate>> {
    if a.qs != b.qs {
        return None;
    }
    // Fuse adjacent non-symbolic phase gates of the same axis.
    if a.t == b.t
        && (a.t == ZPhase || a.t == XPhase)
        && a.vars == Parity::zero()
        && b.vars == Parity::zero()
    {
        let sum = a.phase + b.phase;
        return if sum.is_zero() {
            Some(None)
        } else {
            Some(Some(Gate {
                t: a.t,
                qs: a.qs.clone(),
                phase: sum,
                vars: Parity::zero(),
            }))
        };
    }
    // Cancel a gate against its adjoint. Covers self-inverse gates
    // (HAD, NOT, Z, CNOT, CZ, SWAP, TOFF, CCZ, ...) and the S/Sdg, T/Tdg pairs.
    let mut adj = a.clone();
    adj.adjoint();
    if adj == *b {
        Some(None)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::circuit::Circuit;
    use crate::gate::GType::*;
    use crate::gate::Gate;
    use crate::phase::Phase;

    fn run(gs: Vec<Gate>, nq: usize) -> Vec<Gate> {
        let mut c = Circuit::new(nq);
        for g in gs {
            c.push(g);
        }
        c.optimize();
        c.gates.into_iter().collect()
    }

    #[test]
    fn had_had_cancels() {
        assert!(run(vec![Gate::new(HAD, vec![0]), Gate::new(HAD, vec![0])], 1).is_empty());
    }

    #[test]
    fn cnot_cnot_cancels() {
        assert!(run(
            vec![Gate::new(CNOT, vec![0, 1]), Gate::new(CNOT, vec![0, 1])],
            2
        )
        .is_empty());
    }

    #[test]
    fn s_sdg_cancels() {
        assert!(run(vec![Gate::new(S, vec![0]), Gate::new(Sdg, vec![0])], 1).is_empty());
    }

    #[test]
    fn zphase_fuses() {
        let out = run(
            vec![
                Gate::new_with_phase(ZPhase, vec![0], (1, 2)),
                Gate::new_with_phase(ZPhase, vec![0], (1, 2)),
            ],
            1,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].phase, Phase::from((1, 1)));
    }

    #[test]
    fn zphase_fuses_to_zero() {
        assert!(run(
            vec![
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
                Gate::new_with_phase(ZPhase, vec![0], (-1, 4)),
            ],
            1,
        )
        .is_empty());
    }

    #[test]
    fn fuse_across_commuting_gate() {
        // ZPhase q0 ; HAD q1 ; ZPhase q0  ->  fused ZPhase(1/2) q0 ; HAD q1
        let out = run(
            vec![
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
                Gate::new(HAD, vec![1]),
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
            ],
            2,
        );
        assert_eq!(out.len(), 2);
        assert!(out
            .iter()
            .any(|g| g.t == ZPhase && g.phase == Phase::from((1, 2))));
        assert!(out.iter().any(|g| g.t == HAD));
    }

    #[test]
    fn barrier_blocks_fusion() {
        let out = run(
            vec![
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
                Gate::new(Measure, vec![0]),
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
            ],
            1,
        );
        assert_eq!(out.len(), 3);
    }

    #[test]
    fn no_false_cancel_overlap() {
        // HAD q0 ; CNOT q0,q1 : share q0, not both diagonal -> blocked, unchanged.
        let out = run(
            vec![Gate::new(HAD, vec![0]), Gate::new(CNOT, vec![0, 1])],
            2,
        );
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn chain_of_quarters() {
        // four 1/4 phases -> a single phase of 1
        let out = run(
            vec![
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
                Gate::new_with_phase(ZPhase, vec![0], (1, 4)),
            ],
            1,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].phase, Phase::from((1, 1)));
    }
}

#[cfg(test)]
mod semantic_tests {
    use crate::circuit::Circuit;
    use crate::equality::equal_circuit_tensor;

    /// `optimize` must never change the linear map a circuit represents, and
    /// must never increase the gate count. Checked on random Clifford+T
    /// circuits via quizx's own tensor comparison. Kept small (3 qubits) so it
    /// runs quickly in CI.
    #[test]
    fn optimize_preserves_semantics() {
        for seed in 0..30u64 {
            let mut b = Circuit::random();
            b.seed(seed).qubits(3).depth(30).clifford_t(0.3);
            let c = b.build();
            let mut opt = c.clone();
            opt.optimize();
            assert!(
                equal_circuit_tensor(&c, &opt),
                "seed {seed}: linear map changed"
            );
            assert!(
                opt.num_gates() <= c.num_gates(),
                "seed {seed}: gate count grew"
            );
        }
    }
}
