//! Json encoding for scalar values.
//!
//! This definition is compatible with the `pyzx` JSON format for scalars.

use std::f64::consts::PI;

use num::{Complex, One, Zero};

use super::phase::PhaseOptions;
use super::{JsonError, JsonPhase, JsonScalar};
use crate::phase::Phase;
use crate::scalar::*;

impl From<&Scalar4> for JsonScalar {
    fn from(value: &Scalar4) -> Self {
        let phase_options = PhaseOptions {
            ignore_approx: true,
            ignore_pi: true,
            limit_denom: Some(256),
            ..Default::default()
        };

        match value.exact_phase_and_sqrt2_pow() {
            Some((phase, pow)) => JsonScalar {
                power2: pow,
                phase: JsonPhase::from_phase(phase, phase_options),
                floatfactor: 1.0,
                is_zero: false,
                ..Default::default()
            },
            None => {
                let complex: Complex<f64> = value.complex_value();
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

impl From<Scalar4> for JsonScalar {
    fn from(value: Scalar4) -> Self {
        JsonScalar::from(&value)
    }
}

impl TryFrom<&JsonScalar> for Scalar4 {
    type Error = JsonError;

    fn try_from(value: &JsonScalar) -> Result<Self, Self::Error> {
        if value.is_unknown {
            // TODO: Unknown scalar flag?
            return Ok(Scalar4::one());
        }

        if value.is_zero {
            return Ok(Scalar4::zero());
        }

        let phase = value.phase.to_phase()?.unwrap_or(Phase::zero());
        let mut s = Scalar4::from(phase);

        if value.power2 != 0 {
            s.mul_sqrt2_pow(value.power2);
        }

        if !value.floatfactor.is_zero() {
            s *= Scalar4::from(value.floatfactor);
        }

        for p in &value.phasenodes {
            let p = p.to_phase()?.unwrap_or(Phase::zero());
            s *= Scalar4::one_plus_phase(p);
        }

        Ok(s)
    }
}

impl TryFrom<JsonScalar> for Scalar4 {
    type Error = JsonError;

    fn try_from(value: JsonScalar) -> Result<Self, Self::Error> {
        Scalar4::try_from(&value)
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
    #[case(Scalar4::zero())]
    #[case(Scalar4::one())]
    #[case(Scalar4::from_phase(1))]
    #[case(Scalar4::from_phase((1,2)))]
    #[case(Scalar4::from_phase((-1,2)))]
    #[case(Scalar4::real(2.0))]
    #[case(Scalar4::complex(1.0, 1.0))]
    #[case(Scalar4::new([0, 1, 0, -1], 3))]
    #[case(Scalar4::new([0, 7, 0, 7], -2))]
    #[case(Scalar4::new([-2, 0, -2, 0], 0))]
    #[case(Scalar4::new([2, 0, -2, 0], 30))]
    #[case(Scalar4::new([2, 0, 0, 0], -10))]
    fn scalar_roundtrip(#[case] scalar: Scalar4) -> Result<(), JsonError> {
        let json_scalar = JsonScalar::from(&scalar);
        let decoded: Scalar4 = Scalar4::try_from(&json_scalar)?;
        assert_abs_diff_eq!(scalar, decoded);

        Ok(())
    }
}
