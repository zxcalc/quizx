# Reproducing the stabiliser rank data

The data from the paper "Simulating quantum circuits with ZX-calculus reduced stabiliser decompositions" was produced using the `quizx` library and the shell scripts in the [scripts](scripts) directory. You'll to have Rust installed. A convenient way to get it on most systems is via [rustup](https://rustup.rs/).

To reproduce, first grab a version of `quizx`. Either get the latest version with:

    git clone https://github.com/zxlang/quizx.git

or the exact version used in the paper with:

    git clone --branch stabrank-v1 https://github.com/zxlang/quizx.git


In either case, build `quizx` via:

    cd quizx/quizx
    cargo build --release

The data itself is generated using two binaries, `pauli_gadget_stabrank` and `hidden_shift_stabrank`. These take several parameters as input, including a random seed, and produce a single row of CSV data (in its own file) as output. The easiest way to run these is by using the shell scripts found in the [scripts](scripts) folder. To produce the exact data used in the paper, you can run these commands from the root of git repo:

    cd scripts
    ./pauli_gadget_data.sh
    ./hidden_shift_data.sh

These take a long time to run, and benefit from having many cores available.
