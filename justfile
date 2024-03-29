# List the available commands
help:
    @just --list --justfile {{justfile()}}

# Prepare the environment for development, installing all the dependencies and
# setting up the pre-commit hooks.
setup:
    poetry install
    poetry run pre-commit install -t pre-commit

# Run the pre-commit checks.
check:
    poetry run pre-commit run --all-files

# Build the project.
build language="[rust|python]": (_run_lang language \
        "cargo build" \
        "poetry run maturin develop"
    )

# Run all the tests.
test language="[rust|python]": (_run_lang language \
        "cargo test --all-features" \
        "poetry run maturin develop && poetry run pytest"
    )

# Auto-fix all clippy warnings.
fix language="[rust|python]": (_run_lang language \
        "cargo clippy --all-targets --all-features --workspace --fix --allow-staged --allow-dirty" \
        "poetry run ruff check --fix"
    )

# Format the code.
format language="[rust|python]": (_run_lang language \
        "cargo fmt" \
        "poetry run ruff format"
    )

# Load a shell with all the dependencies installed
shell:
    poetry shell


# Runs a rust and a python command, depending on the `language` variable.
#
# If `language` is set to `rust` or `python`, only run the command for that language.
# Otherwise, run both commands.
_run_lang language rust_cmd python_cmd:
    #!/usr/bin/env bash
    set -euo pipefail
    if [ "{{ language }}" = "rust" ]; then
        set -x
        {{ rust_cmd }}
    elif [ "{{ language }}" = "python" ]; then
        set -x
        {{ python_cmd }}
    else
        set -x
        {{ rust_cmd }}
        {{ python_cmd }}
    fi
