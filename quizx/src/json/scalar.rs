//! Json encoding for scalar values.
//!
//! This definition is compatible with the `pyzx` JSON format for scalars.

use std::f64::consts::PI;

use num::{One, Zero};

use crate::phase::Phase;
use crate::scalar::{Coeffs, FromPhase, Scalar};

use super::{JsonPhase, JsonScalar};

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
            Scalar::Exact(pow, _) => {
                JsonScalar {
                    // The json format encodes a power of sqrt(2), whereas `pow` here is a power of 2.
                    power2: pow + 2,
                    phase: JsonPhase::from_phase(scalar.phase(), false),
                    is_zero: scalar.is_zero(),
                    ..Default::default()
                }
            }
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
    pub fn to_scalar<C: Coeffs>(&self) -> Scalar<C> {
        if self.is_unknown {
            // TODO: Unknown scalar flag?
            return Scalar::one();
        }

        if self.is_zero {
            return Scalar::zero();
        }

        let phase = self.phase.to_phase().unwrap_or(Phase::zero());
        let mut s = Scalar::from_phase(phase);

        if self.power2 != 0 {
            s.mul_sqrt2_pow(self.power2);
        }

        if self.floatfactor != 0.0 {
            s *= Scalar::real(self.floatfactor);
        }

        for p in &self.phasenodes {
            let p = p.to_phase().unwrap_or(Phase::zero());
            s *= Scalar::one_plus_phase(p);
        }

        s
    }
}
