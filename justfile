# List the available commands
help:
    @just --list --justfile {{justfile()}}

# Prepare the environment for development, installing all the dependencies and
# setting up the pre-commit hooks.
setup:
    uv run pre-commit install -t pre-commit

# Run the pre-commit checks.
check:
    uv run pre-commit run --all-files

# Build the project.
build:
    cargo build
    uv run maturin develop

# Run the tests.
test-rust:
    cargo test --all-features
# Run the tests.
test-python:
    uv run maturin develop && uv run pytest

# Auto-fix all lint warnings.
fix-rust:
    cargo clippy --all-targets --all-features --workspace --fix --allow-staged --allow-dirty
# Auto-fix all lint warnings.
fix-python:
    uv run ruff check --fix

# Format the rust code.
format-rust:
    cargo fmt
# Format the rust code.
format-python:
    uv run ruff format
