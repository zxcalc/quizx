//! Json encoding for scalar values.
//!
//! This definition is compatible with the `pyzx` JSON format for scalars.

use std::f64::consts::PI;

use num::{One, Zero};

use crate::phase::Phase;
use crate::scalar::{Coeffs, FromPhase, Scalar};

use super::{JsonError, JsonPhase, JsonScalar};

impl JsonScalar {
    /// Encode a scalar.
    pub fn from_scalar<C: Coeffs>(scalar: &Scalar<C>) -> Self {
        match scalar {
            Scalar::Float(complex) => {
                let (r, theta) = complex.to_polar();
                // Encoding `theta` as a `Phase` here converts it to a fractional value,
                // which may cause a loss of precision.
                let phase = JsonPhase::from_phase(theta / PI, false);
                JsonScalar {
                    phase,
                    floatfactor: r,
                    is_zero: scalar.is_zero(),
                    ..Default::default()
                }
            }
            Scalar::Exact(pow, _) => JsonScalar {
                power2: *pow,
                phase: JsonPhase::from_phase(scalar.phase(), false),
                is_zero: scalar.is_zero(),
                ..Default::default()
            },
        }
    }

    /// Returns an scalar marked as "unknown".
    pub fn unknown() -> Self {
        JsonScalar {
            is_unknown: true,
            ..Default::default()
        }
    }

    /// Decode the json into a [`Scalar`].
    pub fn to_scalar<C: Coeffs>(&self) -> Result<Scalar<C>, JsonError> {
        if self.is_unknown {
            // TODO: Unknown scalar flag?
            return Ok(Scalar::one());
        }

        if self.is_zero {
            return Ok(Scalar::zero());
        }

        let phase = self.phase.to_phase()?.unwrap_or(Phase::zero());
        let mut s = Scalar::from_phase(phase);

        if self.power2 != 0 {
            s.mul_sqrt2_pow(self.power2);
        }

        if !self.floatfactor.is_zero() {
            s *= Scalar::real(self.floatfactor);
        }

        for p in &self.phasenodes {
            let p = p.to_phase()?.unwrap_or(Phase::zero());
            s *= Scalar::one_plus_phase(p);
        }

        Ok(s)
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::scalar::ScalarN;

    use super::*;

    #[rstest]
    #[case(ScalarN::zero())]
    #[case(ScalarN::one())]
    #[case(ScalarN::from_phase(1))]
    #[case(ScalarN::from_phase((1,2)))]
    #[case(ScalarN::from_phase((-1,2)))]
    #[case(ScalarN::real(2.0))]
    #[case(ScalarN::complex(1.0, 1.0))]
    fn scalar_roundtrip(#[case] scalar: ScalarN) -> Result<(), JsonError> {
        println!(
            "initial: {scalar:?}. As complex: {}",
            scalar.complex_value()
        );
        let json_scalar = JsonScalar::from_scalar(&scalar);
        println!("encoded: {json_scalar:?}");
        let decoded: ScalarN = json_scalar.to_scalar()?;
        println!(
            "decoded: {decoded:?}. As complex: {}",
            decoded.complex_value()
        );
        assert!(decoded.approx_eq(&scalar, 1e-6));

        Ok(())
    }
}
