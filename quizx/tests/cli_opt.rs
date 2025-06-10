#[cfg(test)]
mod test {
    use assert_cmd::Command;
    use predicates::str::contains;
    use rstest::{fixture, rstest};

    const CIRC: &str = "../circuits/small/mod5_4.qasm";

    #[fixture]
    fn cmd() -> Command {
        let mut cmd = Command::cargo_bin("quizx").unwrap();
        cmd.arg("opt");
        cmd
    }

    #[rstest]
    fn default(mut cmd: Command) {
        cmd.arg(CIRC).assert().success();
    }

    #[rstest]
    fn full(mut cmd: Command) {
        cmd.arg(CIRC).arg("--full").assert().success();
    }

    #[rstest]
    fn flow(mut cmd: Command) {
        cmd.arg(CIRC).arg("--flow").assert().success();
    }

    #[rstest]
    fn clifford(mut cmd: Command) {
        cmd.arg(CIRC).arg("--clifford").assert().success();
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
        cmd.arg("--full")
            .arg("--clifford")
            .assert()
            .failure()
            .stderr(contains(
                "the argument '--full' cannot be used with '--clifford'",
            ));
    }
}
