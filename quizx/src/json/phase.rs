// QuiZX - Rust library for quantum circuit rewriting and optimisation
//         using the ZX-calculus
// Copyright (C) 2021 - Aleks Kissinger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! Methods for converting phases.

use super::{JsonError, JsonPhase};
use crate::phase::utils::limit_denominator;
use crate::phase::Phase;

use num::{FromPrimitive, One, Rational64, Zero};

/// A set of options for encoding and decoding phases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhaseOptions {
    /// Ignore this value when encoding phases.
    ///
    /// If set to `Some`, then the value will be encoded as an empty string.
    pub ignore_value: Option<Phase>,
    /// If set to `true`, then encoded phases will not include the '~' symbol for approximate values.
    pub ignore_approx: bool,
    /// If set to `true`, then encoded phases will not include the 'pi' symbol.
    pub ignore_pi: bool,
    /// Limit the denominator of the phase to this value.
    pub limit_denom: Option<i64>,
}

impl Default for PhaseOptions {
    fn default() -> Self {
        Self {
            ignore_value: None,
            ignore_approx: false,
            ignore_pi: false,
            limit_denom: Some(256),
        }
    }
}

impl JsonPhase {
    /// Encode a vertex phase.
    ///
    /// See [`PhaseOptions`] for encoding options.
    pub fn from_phase(phase: impl Into<Phase>, options: PhaseOptions) -> Self {
        // This is directly ported from pyzx,
        // trying to match its behaviour as closely as possible.
        let phase = phase.into();

        if let Some(ignore_value) = options.ignore_value {
            if phase == ignore_value {
                return Self("".to_string());
            }
        }

        let phase: Rational64 = phase.to_rational();

        if phase.is_zero() {
            return Self("0".to_string());
        }

        let mut simstr: &str = "";
        let mut phase = phase;
        if let Some(limit_denom) = options.limit_denom {
            if *phase.denom() > limit_denom {
                if !options.ignore_approx {
                    simstr = "~";
                }
                phase = limit_denominator(phase, limit_denom);
            }
        };

        // NOTE: We could insert π instead of "pi" here, but
        // `pyzx.json_to_graph` panics when it sees π.
        //
        // That decoder method is deprecated in `pyzx 0.8.0` (replaced by
        // `pyzx.Graph.from_json`), so we should be able to remove this
        // workaround in the future.
        let numer = match (options.ignore_pi, *phase.numer()) {
            (false, 1) => "pi".to_string(),
            (false, -1) => "-pi".to_string(),
            (false, n) => format!("{}*pi", n),
            (true, n) => format!("{}", n),
        };

        let denom = match *phase.denom() {
            1 => "".to_string(),
            d => format!("/{}", d),
        };

        Self(format!("{simstr}{numer}{denom}"))
    }

    /// Decode a vertex phase.
    ///
    /// Variables are not currently supported.
    ///
    /// Returns `None` if the string is empty or if it contains an invalid value.
    pub fn to_phase(&self) -> Result<Option<Phase>, JsonError> {
        if self.0.is_empty() {
            return Ok(None);
        }

        // Helper function to return when the phase is invalid.
        let phase_error = || JsonError::InvalidPhase {
            phase: self.0.clone(),
        };

        //  // trying to match its behaviour as closely as possible.
        let s: String = self
            .0
            .chars()
            .map(|c| c.to_ascii_lowercase())
            .filter(|c| !c.is_whitespace())
            .filter(|&c| c != 'π')
            .filter(|&c| c != '~')
            .collect::<String>();
        let s = s.replace(r#"\pi"#, "");
        let s = s.replace("pi", "");
        let s = s.as_str();

        // Drop dangling '*' from removing the "pi".
        let s = s.trim_start_matches('*').trim_end_matches('*');

        if s.is_empty() {
            // The phase was just "pi"
            return Ok(Some(Phase::one()));
        }
        if s == "-" {
            return Ok(Some(-Phase::one()));
        }
        if s.contains('.') || s.contains('e') {
            let f: f64 = s.parse().map_err(|_| phase_error())?;
            let phase: Phase = f.into();
            let phase = phase.limit_denominator(256);
            return Ok(Some(phase));
        }
        if s.contains('/') {
            let mut parts = s.split('/');
            let num: &str = parts.next().unwrap().trim_end_matches('*');
            let den: i64 = parts.next().unwrap().parse().map_err(|_| phase_error())?;
            return Ok(Some(
                match num {
                    "" => Rational64::new(1, den),
                    "-" => Rational64::new(-1, den),
                    _ => Rational64::new(num.parse().map_err(|_| phase_error())?, den),
                }
                .into(),
            ));
        }

        let n: i64 = s.parse().map_err(|_| phase_error())?;
        let r: Rational64 = Rational64::from_i64(n).ok_or_else(phase_error)?;
        Ok(Some(r.into()))
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(0, "0")]
    #[case(1, "pi")]
    #[case((1, 2), "pi/2")]
    #[case((1, 3), "pi/3")]
    #[case((-1, 2), "-pi/2")]
    #[case((-1, 3), "-pi/3")]
    #[case((2, 3), "2*pi/3")]
    fn test_from_phase(#[case] phase: impl Into<Phase>, #[case] expected: &str) {
        let phase = phase.into();
        let json_phase = JsonPhase::from_phase(phase, Default::default());
        assert_eq!(json_phase.0, expected);
    }

    #[rstest]
    #[case("0", 0)]
    #[case("1", 1)]
    #[case("1/2", (1, 2))]
    #[case("1/3", (1, 3))]
    #[case("-1/2", (-1, 2))]
    #[case("-1", 1)]
    #[case("pi", 1)]
    #[case("-pi", 1)]
    #[case("pi/3", (1, 3))]
    #[case("-pi/3", (-1, 3))]
    #[case("1/3 * pi", (1, 3))]
    #[case("2*pi/3", (2, 3))]
    #[case("-0.3333333333333333*pi", (-1, 3))]
    #[case("1*π", 1)]
    #[case("π", 1)]
    #[case("~-pi/2", (-1, 2))]
    #[case(r#"7\pi/4"#, (7, 4))]
    fn test_to_phase(#[case] s: &str, #[case] expected: impl Into<Phase>) {
        let expected = expected.into();
        let json_phase = JsonPhase(s.to_string());
        let phase = json_phase.to_phase().unwrap().unwrap_or(Phase::zero());
        assert_eq!(phase, expected);
    }
}
