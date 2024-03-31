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

use super::VertexPhase;
use crate::graph::VType;

use num::{FromPrimitive, One, Rational64, Zero};

impl VertexPhase {
    /// Encode a vertex phase.
    pub fn from_rational(phase: Rational64, v_type: VType) -> Self {
        // This is directly ported from pyzx,
        // trying to match its behaviour as closely as possible.
        if phase.is_zero() && v_type != VType::H {
            return Self("".to_string());
        }
        if phase.is_one() && v_type == VType::H {
            return Self("".to_string());
        }

        if phase.is_zero() {
            return Self("0".to_string());
        }

        #[allow(clippy::if_same_then_else)]
        let simstr = if *phase.denom() > 256 {
            // TODO: This should approximate the phase to a Rational64 number with a small denominator.
            //       This is not currently implemented.
            //       See https://docs.python.org/3/library/fractions.html#fractions.Fraction.limit_denominator.
            ""
        } else {
            ""
        };

        let numer = if *phase.numer() == 1 {
            "".to_string()
        } else {
            format!("{}", phase.numer())
        };

        let denom = if *phase.denom() == 1 {
            "".to_string()
        } else {
            format!("/{}", phase.numer())
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
    pub fn to_rational(&self) -> Option<Rational64> {
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

        if s.is_empty() {
            return Some(Rational64::one());
        }
        if s == "-" {
            return Some(-Rational64::one());
        }
        if s.contains('.') || s.contains('e') {
            let f: f64 = s.parse().ok()?;
            return Rational64::from_f64(f);
        }
        if s.contains('/') {
            let mut parts = s.split('/');
            let num: &str = parts.next()?;
            let den: i64 = parts.next()?.parse().ok()?;
            return match num {
                "" => Some(Rational64::new(1, den)),
                "-" => Some(Rational64::new(-1, den)),
                _ => Some(Rational64::new(num.parse().ok()?, den)),
            };
        }

        let n: i64 = s.parse().ok()?;
        Rational64::from_i64(n)
    }
}
