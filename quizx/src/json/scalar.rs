//! Json encoding for scalar values.
//!
//! This definition is compatible with the `pyzx` JSON format for scalars.

use std::f64::consts::PI;

use num::{Complex, One, Zero};

use super::phase::PhaseOptions;
use super::{JsonError, JsonPhase, JsonScalar};
use crate::fscalar::*;
use crate::phase::Phase;

impl From<&FScalar> for JsonScalar {
    fn from(value: &FScalar) -> Self {
        let phase_options = PhaseOptions {
            ignore_approx: true,
            ignore_pi: true,
            limit_denom: Some(256),
            ..Default::default()
        };

        match value.exact_phase_and_sqrt2_pow() {
            Some((phase, pow)) => JsonScalar {
                power2: pow as i32,
                phase: JsonPhase::from_phase(phase, phase_options),
                floatfactor: 1.0,
                is_zero: false,
                ..Default::default()
            },
            None => {
                let complex: Complex<f64> = value.into();
                let (r, theta) = complex.to_polar();
                JsonScalar {
                    power2: 0,
                    phase: JsonPhase::from_phase(theta / PI, phase_options),
                    floatfactor: r,
                    is_zero: value.is_zero(),
                    ..Default::default()
                }
            }
        }
    }
}

impl From<FScalar> for JsonScalar {
    fn from(value: FScalar) -> Self {
        JsonScalar::from(&value)
    }
}

impl TryFrom<&JsonScalar> for FScalar {
    type Error = JsonError;

    fn try_from(value: &JsonScalar) -> Result<Self, Self::Error> {
        if value.is_unknown {
            // TODO: Unknown scalar flag?
            return Ok(FScalar::one());
        }

        if value.is_zero {
            return Ok(FScalar::zero());
        }

        let phase = value.phase.to_phase()?.unwrap_or(Phase::zero());
        let mut s = FScalar::from(phase);

        if value.power2 != 0 {
            s.mul_sqrt2_pow(value.power2);
        }

        if !value.floatfactor.is_zero() {
            s *= FScalar::from(value.floatfactor);
        }

        for p in &value.phasenodes {
            let p = p.to_phase()?.unwrap_or(Phase::zero());
            s *= FScalar::one_plus_phase(p);
        }

        Ok(s)
    }
}

impl TryFrom<JsonScalar> for FScalar {
    type Error = JsonError;

    fn try_from(value: JsonScalar) -> Result<Self, Self::Error> {
        FScalar::try_from(&value)
    }
}

impl JsonScalar {
    /// Returns an scalar marked as "unknown".
    pub fn unknown() -> Self {
        JsonScalar {
            is_unknown: true,
            ..Default::default()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use rstest::rstest;

    #[rstest]
    #[case(FScalar::zero())]
    #[case(FScalar::one())]
    #[case(FScalar::from_phase(1))]
    #[case(FScalar::from_phase((1,2)))]
    #[case(FScalar::from_phase((-1,2)))]
    #[case(FScalar::real(2.0))]
    #[case(FScalar::complex(1.0, 1.0))]
    #[case(FScalar::dyadic(3, [0, 1, 0, -1]))]
    #[case(FScalar::dyadic(-2, [0, 7, 0, 7]))]
    #[case(FScalar::dyadic(0, [-2, 0, -2, 0]))]
    #[case(FScalar::dyadic(30, [2, 0, -2, 0]))]
    #[case(FScalar::dyadic(-10, [2, 0, 0, 0]))]
    fn scalar_roundtrip(#[case] scalar: FScalar) -> Result<(), JsonError> {
        let json_scalar = JsonScalar::from(&scalar);
        let decoded: FScalar = FScalar::try_from(&json_scalar)?;
        assert_abs_diff_eq!(scalar, decoded);

        Ok(())
    }
}
