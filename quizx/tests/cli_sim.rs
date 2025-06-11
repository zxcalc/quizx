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
