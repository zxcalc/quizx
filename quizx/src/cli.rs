//! The QuiZX command line interface.

use clap::{crate_version, Parser};

pub mod opt;
pub mod sim;

/// CLI arguments.
#[derive(Parser, Debug)]
#[clap(version = crate_version!(), long_about = None)]
#[clap(about = "QuiZX command line interface")]
pub enum Cli {
    /// Run the circuit optimizer.
    Opt(opt::OptArgs),
    /// Run the circuit simulator.
    Sim(sim::SimArgs),
}

/// Error type for the CLI.
#[derive(Debug, derive_more::Display, derive_more::From)]
pub enum CliError {
    /// Error reading or writing files.
    #[display("IO error: {_0}")]
    IO(std::io::Error),
    /// Error parsing a QASM file.
    #[display("Error parsing input circuit: {_0}")]
    CircuitParse(String),
    /// Provided bit/Pauli string has the wrong length
    #[display("Circuit has {_0} qubits, but the provided {_2} string has length {_1}")]
    StringWrongLen(usize, usize, String),
}

impl Cli {
    pub fn run(self) -> Result<(), CliError> {
        match self {
            Cli::Opt(args) => args.run(),
            Cli::Sim(args) => args.run(),
        }
    }
}
