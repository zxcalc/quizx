//! Utility functions for phase values.

use num::rational::Ratio;
use num::Integer;

/// Approximate a fraction to a Rational64 number with a small denominator.
///
/// Mimics python's
/// <https://docs.python.org/3/library/fractions.html#fractions.Fraction.limit_denominator>
///
/// # Panics
///
/// Panics if `max_denom` is 0.
///
/// # Example
/// ```
/// # use quizx::phase::utils::limit_denominator;
/// # use num::FromPrimitive;
/// # use num::rational::Rational64;
/// assert_eq!(
///     limit_denominator(Rational64::from_f64(3.141592653589793).unwrap(), 10),
///     Rational64::new(22, 7)
/// );
/// assert_eq!(
///     limit_denominator(Rational64::from_f64(3.141592653589793).unwrap(), 100),
///     Rational64::new(311, 99)
/// );
/// assert_eq!(
///     limit_denominator(Rational64::new(4321, 8765), 10000),
///     Rational64::new(4321, 8765)
/// );
/// ```
pub fn limit_denominator<T>(fraction: Ratio<T>, max_denom: T) -> Ratio<T>
where
    T: Clone + Integer,
{
    // Using cpython's limit_denominator as a reference.
    //
    // https://github.com/python/cpython/blob/bfc57d43d8766120ba0c8f3f6d7b2ac681a81d8a/Lib/fractions.py#L357
    if max_denom <= T::one() {
        panic!("max_denom must be greater than 1");
    }

    let mut numer = fraction.numer().clone();
    let mut denom = fraction.denom().clone();
    if denom <= max_denom {
        return fraction;
    }

    // Upper and lower bounds for the approximation.
    let mut numer_0 = T::zero();
    let mut denom_0 = T::one();
    let mut numer_1 = T::one();
    let mut denom_1 = T::zero();
    loop {
        let a = numer.div_floor(&denom);
        let new_numer = denom_0.clone() + a.clone() * denom_1.clone();
        if new_numer > max_denom {
            break;
        }
        let new_denom = numer_0 + a.clone() * numer_1.clone();
        numer_0 = numer_1;
        denom_0 = denom_1;
        numer_1 = new_denom;
        denom_1 = new_numer;

        let tmp_denom = numer - a * denom.clone();
        numer = denom;
        denom = tmp_denom;
    }
    let k = (max_denom - denom_0.clone()).div_floor(&denom_1);

    if (T::one() + T::one()) * denom * (denom_0.clone() + k.clone() * denom_1.clone())
        <= *fraction.denom()
    {
        Ratio::new_raw(numer_1, denom_1)
    } else {
        Ratio::new_raw(numer_0 + k.clone() * numer_1, denom_0 + k * denom_1)
    }
}
