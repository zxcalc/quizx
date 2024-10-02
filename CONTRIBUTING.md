# Welcome to the quizx development guide <!-- omit in toc -->

This guide is intended to help you get started with developing quizx.

If you find any errors or omissions in this document, please [open an issue](https://github.com/zxlang/quizx/issues/new)!

## #️⃣ Setting up the development environment

To develop the rust library, you will need to have the rust toolchain installed. You can install it by following the instructions at [rust-lang.org/tools/install](https://www.rust-lang.org/tools/install).

If you are using VSCode, you can install the `rust-analyzer` extension to get code completion and other features.

To develop the python library, you will additionally need the `uv` package manager. You can install it by following the instructions at [docs.astral.sh/uv/](https://docs.astral.sh/uv/).

Finally, we provide a `just` file to help manage the common development workflow. You can install `just` by following the instructions at [just.systems](https://just.systems/).

Once you have these installed, run `just setup` to download the necessary dependencies and set up some pre-commit hooks.

Run `just` to see all available commands.

## 🚀 Building the project

There is a miscellaneous collection of rust programs written using the library,
found in `examples/`. To execute these programs, run:

```bash
cargo run --release --example <program_name>
```

To build the python library, run:

```bash
# Setup the dependencies and build the python library
uv run maturin develop
# The library will now be available for import in python
uv run python -c "import quizx"
```

## 🏃 Running the tests

To compile and test the code, run:

```bash
just test
# or, to test only the rust code or the python code
just test rust
just test python
```

## 💅 Coding Style

The `rustfmt` and `ruff` tools are used to enforce a consistent coding for rust
and python code, respectively. The CI will fail if the code is not formatted
correctly.

To format your code, run:

```bash
just format
```

We also use various linters to catch common mistakes and enforce best practices. To run these, use:

```bash
just check
```

To quickly fix common issues, run:

```bash
just fix
# or, to fix only the rust code or the python code
just fix rust
just fix python
```

## 🌐 Contributing to quizx

We welcome contributions to quizx! Please open [an issue](https://github.com/CQCL/hugr/issues/new) or [pull request](https://github.com/CQCL/hugr/compare) if you have any questions or suggestions.

PRs should be made against the `master` branch, and should pass all CI checks before being merged.

Please title your PRs using the [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) format;
```
<type>(<scope>)!: <description>
```
Where the scope is optional, and the `!` is only included if this is a semver breaking change that requires a major version bump.

We accept the following contribution types:

- feat: New features.
- fix: Bug fixes.
- docs: Improvements to the documentation.
- style: Formatting, missing semi colons, etc; no code change.
- refactor: Refactoring code without changing behaviour.
- perf: Code refactoring focused on improving performance.
- test: Adding missing tests, refactoring tests; no production code change.
- ci: CI related changes. These changes are not published in the changelog.
- chore: Updating build tasks, package manager configs, etc. These changes are not published in the changelog.
- revert: Reverting previous commits.
 