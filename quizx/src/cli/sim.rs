//! The `sim` CLI subcommand.

use clap::{Args, Parser};
use itertools::Itertools;
use num::rational::Ratio;
// use num::Zero;
use rand::{thread_rng, Rng};
use std::error::Error;
use std::fs;
use std::path::PathBuf;

use crate::circuit::Circuit;
use crate::decompose::{Decomposer, Driver};
use crate::fscalar::FScalar;
use crate::graph::{BasisElem, GraphLike, VType};
use crate::simplify;
use crate::vec_graph::Graph;

use super::CliError;

/// Run the circuit simulator.
#[derive(Parser, Debug)]
pub struct SimArgs {
    /// QASM file to simulate.
    input: PathBuf,

    /// Output to a file instead of printing the results.
    #[arg(long, short)]
    out: Option<PathBuf>,

    /// Switch to select the decomposition method. Defaults to `--cats`.
    #[command(flatten)]
    method: Option<SimMethod>,

    /// Switch to select the simulation task. Defaults to `--shots 1`.
    #[command(flatten)]
    task: Option<SimTask>,

    /// Distribute computation across available CPU cores up to the given depth.
    #[arg(long, short)]
    parallel: Option<usize>,
}

impl SimArgs {
    /// Run the `sim` command using the provided arguments.
    pub fn run(self) -> Result<(), CliError> {
        let circ = Circuit::from_file(self.input.to_str().unwrap())?;
        let mut d = self.method.unwrap_or_default().build_decomposer();
        let result = self
            .task
            .unwrap_or_default()
            .run(&circ, &mut d, self.parallel)?;

        if let Some(out_path) = self.out {
            fs::write(out_path, result)?;
        } else {
            println!("{result}");
        }
        Ok(())
    }
}

/// Simulation methods.
#[derive(Args, Debug)]
#[group(multiple = false)]
// Ideally, this should be an enum. Unfortunately, clap doesn't support enum arg
// groups yet, so we have to use a struct where only one field is allowed to be
// set to `true`. See https://github.com/clap-rs/clap/issues/2621.
pub struct SimMethod {
    /// Use the cat state decomposition strategy (default).
    #[arg(long)]
    cats: bool,

    /// Use the BSS decomposition strategy.
    #[arg(long)]
    bss: bool,
}

impl Default for SimMethod {
    fn default() -> Self {
        SimMethod {
            bss: false,
            cats: true,
        }
    }
}

impl SimMethod {
    fn build_decomposer(&self) -> Decomposer<Graph> {
        let driver: Driver = if self.cats {
            Driver::BssWithCats(true)
        } else {
            Driver::BssTOnly(true)
        };
        let mut decomposer = Decomposer::empty();
        decomposer.with_full_simp().with_driver(driver);
        decomposer
    }
}

/// Simulation tasks.
#[derive(Args, Debug)]
#[group(multiple = false)]
// Ideally, this should be an enum. Unfortunately, clap doesn't support enum arg
// groups yet, so we have to use a struct where only one field is allowed to be
// set to `Some`. See https://github.com/clap-rs/clap/issues/2621.
pub struct SimTask {
    /// Sample from the circuit for a given number of shots.
    #[arg(long, short)]
    shots: Option<usize>,

    /// Compute an amplitude.
    #[arg(long = "amplitude", short = 'a', value_parser = parse_bit_string)]
    bit_string: Option<BitString>,

    /// Compute an expectation value.
    #[arg(long = "expval", short = 'e', value_parser = parse_pauli_string)]
    pauli_string: Option<PauliString>,
}

impl Default for SimTask {
    fn default() -> Self {
        SimTask {
            shots: Some(1),
            bit_string: None,
            pauli_string: None,
        }
    }
}

impl SimTask {
    pub fn run(
        &self,
        circ: &Circuit,
        decomposer: &mut Decomposer<Graph>,
        parallel: Option<usize>,
    ) -> Result<String, CliError> {
        if let Some(shots) = self.shots {
            Ok((0..shots)
                .map(|_| sample(circ, decomposer, parallel))
                .join("\n")
                .to_string())
        } else if let Some(ref bit_str) = self.bit_string {
            Ok(format!("{}", amplitude(circ, decomposer, bit_str, parallel)?).to_string())
        } else if let Some(ref pauli_str) = self.pauli_string {
            Ok(format!(
                "{}",
                expectation_value(circ, decomposer, pauli_str, parallel)?
            )
            .to_string())
        } else {
            unreachable!()
        }
    }
}

// Need to wrap the vector into a type alias, otherwise clap tries to do some
// varag parsing magic that breaks our custom parser below.
type BitString = Vec<bool>;

#[derive(Debug, derive_more::Display)]
#[display("'{_0}' is not a valid bit. Expected sequence of 0s and 1s.")]
struct BitStringParseError(char);

impl Error for BitStringParseError {}

fn parse_bit_string(s: &str) -> Result<BitString, BitStringParseError> {
    s.chars()
        .map(|c| match c.to_ascii_uppercase() {
            '0' => Ok(false),
            '1' => Ok(true),
            _ => Err(BitStringParseError(c)),
        })
        .collect()
}

#[derive(Clone, Copy, Debug)]
enum Pauli {
    I,
    X,
    Y,
    Z,
}

// Need to wrap the vector into a type alias, otherwise clap tries to do some
// varag parsing magic that breaks our custom parser below.
type PauliString = Vec<Pauli>;

#[derive(Debug, derive_more::Display)]
#[display("'{_0}' is not a Pauli. Expected one of 'I', 'X', 'Y', 'Z'.")]
struct PauliStringParseError(char);

impl Error for PauliStringParseError {}

fn parse_pauli_string(s: &str) -> Result<PauliString, PauliStringParseError> {
    s.chars()
        .map(|c| match c.to_ascii_uppercase() {
            'I' => Ok(Pauli::I),
            'X' => Ok(Pauli::X),
            'Y' => Ok(Pauli::Y),
            'Z' => Ok(Pauli::Z),
            _ => Err(PauliStringParseError(c)),
        })
        .collect()
}

/// Sample from a circuit by computing marginals via doubling of the diagram.
fn sample(circ: &Circuit, decomposer: &mut Decomposer<Graph>, parallel: Option<usize>) -> String {
    let qs = circ.num_qubits();
    let mut xs: Vec<bool> = vec![];
    let mut rng = thread_rng();
    for _ in 0..qs {
        let mut g: Graph = circ.to_graph();
        g.plug_inputs(&vec![BasisElem::Z0; qs]);
        for x in &xs {
            // Plug removes the output, so we have to keep using index 0
            g.plug_output(0, if *x { BasisElem::Z1 } else { BasisElem::Z0 });
        }
        g.plug_output(0, BasisElem::Z1);
        g.plug(&g.to_adjoint());

        let scalar = decomp_graph(g, decomposer, parallel);
        xs.push(rng.gen_bool(scalar.complex_value().re));
    }
    xs.iter().map(|x| if *x { '1' } else { '0' }).join("")
}

/// Compute an amplitude.
fn amplitude(
    circ: &Circuit,
    decomposer: &mut Decomposer<Graph>,
    bit_str: &BitString,
    parallel: Option<usize>,
) -> Result<f64, CliError> {
    let qs = circ.num_qubits();
    let bit_str = match bit_str.as_slice() {
        [b] => &vec![*b; qs],
        bs if bs.len() == qs => bit_str,
        _ => {
            return Err(CliError::StringWrongLen(
                qs,
                bit_str.len(),
                "bit".to_string(),
            ))
        }
    };

    let mut g: Graph = circ.to_graph();
    g.plug_inputs(&vec![BasisElem::Z0; qs]);
    g.plug_outputs(
        &bit_str
            .iter()
            .map(|x| if *x { BasisElem::Z1 } else { BasisElem::Z0 })
            .collect_vec(),
    );

    let scalar = decomp_graph(g, decomposer, parallel);
    let amp = scalar * scalar.conj();
    Ok(amp.complex_value().re)
}

/// Computes an expectation value by doubling the diagram.
fn expectation_value(
    circ: &Circuit,
    decomposer: &mut Decomposer<Graph>,
    pauli_str: &PauliString,
    parallel: Option<usize>,
) -> Result<f64, CliError> {
    let qs = circ.num_qubits();
    let pauli_str = match pauli_str.as_slice() {
        [p] => &vec![*p; qs],
        ps if ps.len() == qs => pauli_str,
        _ => {
            return Err(CliError::StringWrongLen(
                qs,
                pauli_str.len(),
                "Pauli".to_string(),
            ))
        }
    };

    let mut g: Graph = circ.to_graph();
    g.plug_inputs(&vec![BasisElem::Z0; qs]);
    let g_adj = g.to_adjoint();
    for (i, p) in pauli_str.iter().enumerate() {
        let b = g.outputs()[i];
        let [(v, _)] = g.incident_edge_vec(b).try_into().unwrap();
        match p {
            Pauli::I => {}
            Pauli::X => {
                let x = g.add_vertex_with_phase(VType::X, 1);
                g.remove_edge(v, b);
                g.add_edge(v, x);
                g.add_edge(x, b);
            }
            Pauli::Y => {
                let x = g.add_vertex_with_phase(VType::X, 1);
                let z = g.add_vertex_with_phase(VType::Z, 1);
                g.remove_edge(v, b);
                g.add_edge(v, z);
                g.add_edge(z, x);
                g.add_edge(x, b);
                g.scalar_mut().mul_phase(Ratio::new(1, 2));
            }
            Pauli::Z => {
                let z = g.add_vertex_with_phase(VType::Z, 1);
                g.remove_edge(v, b);
                g.add_edge(v, z);
                g.add_edge(z, b);
            }
        }
    }
    g.plug(&g_adj);

    let scalar = decomp_graph(g, decomposer, parallel);
    Ok(scalar.complex_value().re)
}

/// Run the provided decomposer on a graph.
fn decomp_graph(
    mut g: Graph,
    decomposer: &mut Decomposer<Graph>,
    parallel: Option<usize>,
) -> FScalar {
    simplify::full_simp(&mut g);
    decomposer.set_target(g);
    if let Some(_depth) = parallel {
        decomposer.decompose_parallel().scalar()
    } else {
        decomposer.decompose().scalar()
    }
}

#[cfg(test)]
mod test {
    use assert_cmd::Command;
    use predicates::{ord::eq, str::contains};
    use rstest::{fixture, rstest};

    const CIRC: &str = "../circuits/small/mod5_4.qasm";
    const SAMPLE: &str = "00001\n";

    #[fixture]
    fn cmd() -> Command {
        let mut cmd = Command::cargo_bin("quizx").unwrap();
        cmd.arg("sim");
        cmd
    }

    #[rstest]
    fn default(mut cmd: Command) {
        cmd.arg(CIRC).assert().success().stdout(eq(SAMPLE));
    }

    #[rstest]
    fn samples(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--shots")
            .arg("2")
            .assert()
            .success()
            .stdout(eq([SAMPLE, SAMPLE].join("")));
    }

    #[rstest]
    fn amplitude_single(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--amplitude")
            .arg("0")
            .assert()
            .success()
            .stdout(eq("0\n"));
    }

    #[rstest]
    fn amplitude_long(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--amplitude")
            .arg("00001")
            .assert()
            .success()
            .stdout(eq("1\n"));
    }

    #[rstest]
    fn expectation_single(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--expval")
            .arg("Z")
            .assert()
            .success()
            .stdout(eq("-1\n"));
    }

    #[rstest]
    fn expectation_explicit(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--expval")
            .arg("IXYZZ")
            .assert()
            .success()
            .stdout(eq("0\n"));
    }

    #[rstest]
    fn expectation_lower(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--expval")
            .arg("ixzyy")
            .assert()
            .success()
            .stdout(eq("0\n"));
    }

    #[rstest]
    fn bss(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--bss")
            .assert()
            .success()
            .stdout(eq(SAMPLE));
    }

    #[rstest]
    fn cats(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--cats")
            .assert()
            .success()
            .stdout(eq(SAMPLE));
    }

    #[rstest]
    fn parallel_sample(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--parallel")
            .arg("2")
            .assert()
            .success()
            .stdout(eq(SAMPLE));
    }

    #[rstest]
    fn parallel_amplitude(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--expval")
            .arg("z")
            .arg("--parallel")
            .arg("2")
            .assert()
            .success()
            .stdout(eq("-1\n"));
    }

    #[rstest]
    fn parallel_expectation(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--expval")
            .arg("z")
            .arg("--parallel")
            .arg("2")
            .assert()
            .success()
            .stdout(eq("-1\n"));
    }

    #[rstest]
    fn doesnt_exist(mut cmd: Command) {
        cmd.arg("blah")
            .assert()
            .failure()
            .stderr(contains("Error parsing input circuit: can't read file"));
    }

    #[rstest]
    fn multiple_methods(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--bss")
            .arg("--cats")
            .assert()
            .failure()
            .stderr(contains(
                "the argument '--bss' cannot be used with '--cats'",
            ));
    }

    #[rstest]
    fn multiple_tasks(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--shots")
            .arg("2")
            .arg("--expval")
            .arg("Z")
            .assert()
            .failure()
            .stderr(contains(
                "the argument '--shots <SHOTS>' cannot be used with '--expval <PAULI_STRING>'",
            ));
    }

    #[rstest]
    fn bad_bit(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--amplitude")
            .arg("00210")
            .assert()
            .failure()
            .stderr(contains("invalid value '00210' for '--amplitude <BIT_STRING>': '2' is not a valid bit. Expected sequence of 0s and 1s."));
    }

    #[rstest]
    fn bit_string_len(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--amplitude")
            .arg("010")
            .assert()
            .failure()
            .stderr(contains(
                "Circuit has 5 qubits, but the provided bit string has length 3",
            ));
    }

    #[rstest]
    fn bad_pauli(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--expval")
            .arg("ZXYAI")
            .assert()
            .failure()
            .stderr(contains("invalid value 'ZXYAI' for '--expval <PAULI_STRING>': 'A' is not a Pauli. Expected one of 'I', 'X', 'Y', 'Z'."));
    }

    #[rstest]
    fn pauli_string_len(mut cmd: Command) {
        cmd.arg(CIRC)
            .arg("--expval")
            .arg("ZZZ")
            .assert()
            .failure()
            .stderr(contains(
                "Circuit has 5 qubits, but the provided Pauli string has length 3",
            ));
    }
}
