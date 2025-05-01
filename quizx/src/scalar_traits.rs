use crate::phase::Phase;

/// Produce a number from rational root of -1.
pub trait FromPhase {
    /// Returns a number from a rational phase.
    fn from_phase(p: impl Into<Phase>) -> Self;
    /// Returns the number -1.
    fn minus_one() -> Self;
}

/// Contains the numbers sqrt(2) and 1/sqrt(2), often used for
/// renormalisation of qubit tensors and matrices.
pub trait Sqrt2: Sized {
    /// Return the number sqrt(2).
    fn sqrt2() -> Self {
        Self::sqrt2_pow(1)
    }
    /// Return the number 1/sqrt(2).
    fn one_over_sqrt2() -> Self {
        Self::sqrt2_pow(-1)
    }
    /// Return the p-th power of sqrt(2).
    fn sqrt2_pow(p: i32) -> Self;
}
