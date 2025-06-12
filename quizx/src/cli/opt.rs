//! The `opt` CLI subcommand.

use clap::{Args, Parser};
use std::fs;
use std::path::PathBuf;

use crate::circuit::Circuit;
use crate::extract::ToCircuit;
use crate::simplify;
use crate::vec_graph::Graph;

use super::CliError;

/// Run the circuit optimizer.
#[derive(Parser, Debug)]
pub struct OptArgs {
    /// QASM file to optimize.
    input: PathBuf,

    /// Output to a file instead of printing the result.
    #[arg(long, short)]
    out: Option<PathBuf>,

    /// Switch to select the optimization method. Defaults to `--full`.
    #[command(flatten)]
    method: Option<OptMethod>,
}

impl OptArgs {
    /// Run the `opt` command using the provided arguments.
    pub fn run(self) -> Result<(), CliError> {
        let circ = Circuit::from_file(self.input.to_str().unwrap())?;
        let mut g = circ.to_graph();
        self.method.unwrap_or_default().simp(&mut g);
        let qasm = g
            .to_circuit()
            .expect("Extraction should succeed since we start from a circuit")
            .to_qasm();
        if let Some(out_path) = self.out {
            fs::write(out_path, qasm)?;
        } else {
            println!("{qasm}");
        }
        Ok(())
    }
}

/// Optimization method.
#[derive(Args, Debug)]
#[group(multiple = false)]
// Ideally, this should be an enum. Unfortunately, clap doesn't support enum arg
// groups yet, so we have to use a struct where only one field is allowed to be
// set to `true`. See https://github.com/clap-rs/clap/issues/2621.
pub struct OptMethod {
    /// Optimize using the `full_simp` method (default).
    #[arg(long)]
    full: bool,

    /// Optimize using the `flow_simp` method.
    #[arg(long)]
    flow: bool,

    /// Optimize using the `clifford_simp` method.
    #[arg(long)]
    clifford: bool,
}

impl Default for OptMethod {
    fn default() -> Self {
        OptMethod {
            full: true,
            flow: false,
            clifford: false,
        }
    }
}

impl OptMethod {
    fn simp(&self, g: &mut Graph) {
        if self.full {
            simplify::full_simp(g);
        } else if self.flow {
            simplify::flow_simp(g);
        } else if self.clifford {
            simplify::clifford_simp(g);
        }
    }
}
