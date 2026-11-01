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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_error_display_io() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let err = CliError::IO(io_err);
        let msg = format!("{}", err);
        assert!(msg.contains("IO error"));
        assert!(msg.contains("file not found"));
    }

    #[test]
    fn cli_error_display_parse() {
        let err = CliError::CircuitParse("invalid gate".to_string());
        let msg = format!("{}", err);
        assert!(msg.contains("Error parsing"));
        assert!(msg.contains("invalid gate"));
    }

    #[test]
    fn cli_error_display_string_len() {
        let err = CliError::StringWrongLen(5, 3, "Pauli".to_string());
        let msg = format!("{}", err);
        assert!(msg.contains("5 qubits"));
        assert!(msg.contains("length 3"));
        assert!(msg.contains("Pauli"));
    }

    #[test]
    fn cli_error_from_io() {
        let io_err = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "denied");
        let cli_err: CliError = io_err.into();
        assert!(matches!(cli_err, CliError::IO(_)));
    }
}
