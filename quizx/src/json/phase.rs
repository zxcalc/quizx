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

use super::JsonPhase;
use crate::phase::utils::limit_denominator;
use crate::phase::Phase;

use num::{FromPrimitive, One, Rational64, Zero};

impl JsonPhase {
    /// Encode a vertex phase.
    ///
    /// By default, zero-values are encoded as an empty string.
    /// If `ignore_one` is set to `true`, then phases with value one are encoded as an empty instead.
    pub fn from_phase(phase: impl Into<Phase>, ignore_one: bool) -> Self {
        // This is directly ported from pyzx,
        // trying to match its behaviour as closely as possible.
        let phase = phase.into();

        if !ignore_one && phase.is_zero() {
            return Self("".to_string());
        }
        if ignore_one && phase.is_one() {
            return Self("".to_string());
        }

        Self::from_rational(phase.to_rational())
    }

    /// Encode a phase expressed as a rational number.
    fn from_rational(phase: Rational64) -> Self {
        if phase.is_zero() {
            return Self("0".to_string());
        }

        let (phase, simstr) = if *phase.denom() > 256 {
            let phase = limit_denominator(phase, 256);
            (phase, "~")
        } else {
            (phase, "")
        };

        let numer = match *phase.numer() {
            1 => "".to_string(),
            -1 => "-".to_string(),
            n => format!("{}*", n),
        };

        let denom = match *phase.denom() {
            1 => "".to_string(),
            d => format!("/{}", d),
        };

        // NOTE: We could insert π instead of "pi" here, but
        // `pyzx.json_to_graph` panics when it sees π.
        //
        // That decoder method is deprecated in `pyzx 0.8.0` (replaced by
        // `pyzx.Graph.from_json`), so we should be able to remove this
        // workaround in the future.
        Self(format!("{simstr}{numer}pi{denom}"))
    }

    /// Decode a vertex phase.
    ///
    /// Variables are not currently supported.
    ///
    /// Returns `None` if the string is empty or if it contains an invalid value.
    pub fn to_phase(&self) -> Option<Phase> {
        if self.0.is_empty() {
            return None;
        }

        // This is directly ported from pyzx,
        // trying to match its behaviour as closely as possible.
        let s: String = self
            .0
            .chars()
            .map(|c| c.to_ascii_lowercase())
            .filter(|c| !c.is_whitespace())
            .filter(|&c| c != 'π')
            .collect::<String>();
        let s = s.replace("pi", "");
        let s = s.as_str();

        // Drop dangling '*' from removing the "pi".
        let s = s.trim_start_matches('*').trim_end_matches('*');

        if s.is_empty() {
            // The phase was just "pi"
            return Some(Phase::one());
        }
        if s == "-" {
            return Some(-Phase::one());
        }
        if s.contains('.') || s.contains('e') {
            let f: f64 = s.parse().ok()?;
            let phase: Phase = f.into();
            println!("Limiting denominator of {phase:?}");
            return Some(phase.limit_denominator(256));
        }
        if s.contains('/') {
            let mut parts = s.split('/');
            let num: &str = parts.next()?.trim_end_matches('*');
            let den: i64 = parts.next()?.parse().ok()?;
            return Some(
                match num {
                    "" => Rational64::new(1, den),
                    "-" => Rational64::new(-1, den),
                    _ => Rational64::new(num.parse().ok()?, den),
                }
                .into(),
            );
        }

        let n: i64 = s.parse().ok()?;
        Rational64::from_i64(n).map(Into::into)
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(0, "")]
    #[case(1, "pi")]
    #[case((1, 2), "pi/2")]
    #[case((1, 3), "pi/3")]
    #[case((-1, 2), "-pi/2")]
    #[case((-1, 3), "-pi/3")]
    #[case((2, 3), "2*pi/3")]
    fn test_from_phase(#[case] phase: impl Into<Phase>, #[case] expected: &str) {
        let phase = phase.into();
        println!("phase: {phase:?} expected: {expected}");
        let json_phase = JsonPhase::from_phase(phase, false);
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
    fn test_to_phase(#[case] s: &str, #[case] expected: impl Into<Phase>) {
        let expected = expected.into();
        println!("encoded: {s:?} expected: {expected}");
        let json_phase = JsonPhase(s.to_string());
        let phase = json_phase.to_phase().unwrap();
        assert_eq!(phase, expected);
    }
}
