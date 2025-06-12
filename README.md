# QuiZX: a quick Rust port of PyZX

[![pypi][]](https://pypi.org/project/quizx/)
[![py-version][]](https://pypi.org/project/quizx/)
[![crates][]](https://crates.io/crates/quizx)
[![msrv][]](https://github.com/zxlang/quizx)
[![rs-docs][]](https://docs.rs/quizx)

  [pypi]: https://img.shields.io/pypi/v/quizx
  [py-version]: https://img.shields.io/pypi/pyversions/quizx
  [crates]: https://img.shields.io/crates/v/quizx
  [msrv]: https://img.shields.io/crates/msrv/quizx
  [rs-docs]: https://img.shields.io/docsrs/quizx?label=rust%20docs


QuiZX is a command line tool and [Rust](https://www.rust-lang.org/) library for circuit optimisation and classical simulation using the [ZX-calculus](https://zxcalculus.com). It began life as high-performance Rust port of the Python library [PyZX](https://github.com/zxlang/pyzx) and is largely inter-operable with that library using QuiZX's Python bindings.

## Getting started

There are three ways to use QuiZX: as a command-line tool, a Rust library, or via its Python bindings. Install the command line tool by running:

```bash
cargo install quizx
```

The `quizx` command currently has two sub-commands `quizx opt` for circuit optimisation and `quizx sim` for classical simulation.  For example, to run the circuit optimiser on a QASM file with default options, run:

```bash
quizx opt circuit.qasm -o circuit_opt.qasm 
```

To run the classical simulator for 16 shots, run:

```bash
quizx opt circuit.qasm -s 16
```

Run `quizx help` to get an overview of the available commands then `quizx COMMAND --help` for command-specific help.


## Rust and Python API

To use QuiZX from your Rust code, add it in the usual way to `Cargo.toml`. See the [docs](https://docs.rs/quizx/latest/quizx/) for an overview of the Rust API, as well as mixed bag of examples in the [examples](https://github.com/zxcalc/quizx/tree/master/quizx/examples) folder.

To use QuiZX via Python bindings, simply `pip install quizx`. You may wish to also run `pip install pyzx` to use various functions from PyZX on QuiZX graphs, most notably drawing functions. See [getting_started.ipynb](https://github.com/zxcalc/quizx/blob/master/demos/getting_started.ipynb) for basic usage. QuiZX graphs expose the same public methods as PyZX graphs, so any function in PyZX that expects a `pyzx.graph.BaseGraph` should also work fine (in theory) on QuiZX graphs. This is not thoroughly tested, so if you run into problems, please file an issue.


## Performance

We are in the process of more rigorously benchmarking the performance of QuiZX. You can see the results so far in the [benches](https://github.com/zxcalc/quizx/tree/master/quizx/benches) folder. These are written using the [Criterion](https://docs.rs/criterion/latest/criterion/) benchmarking library and are runnable from the command line using `cargo bench`.

As a very anecdotal example of the performance difference with PyZX, the program `spider_chain` builds a chain of 1 million Z spiders and fuses them all. In PyZX, you can fuse all the spiders in a ZX-diagram as follows:

```python
from pyzx.basicrules import *

success = True
while success:
    success = any(fuse(g, g.edge_s(e), g.edge_t(e)) for e in g.edges())
```

In QuiZX, the Rust code is slightly more verbose, but similar in spirit:
```rust
use quizx::basic_rules::*;

loop {
    match g.find_edge(|v0,v1,_| check_spider_fusion(&g, v0, v1)) {
        Some((v0,v1,_)) => spider_fusion_unchecked(&mut g, v0, v1),
        None => break,
    };
}
```

On my laptop, the PyZX code takes about 98 seconds to fuse 1 million spiders, whereas the QuiZX code takes 17 milliseconds.


## Development

The library is available both as a Rust crate and a Python package. See the [Python Changelog](https://github.com/zxcalc/quizx/blob/master/pybindings/CHANGELOG.md) and the [Rust Changelog](https://github.com/zxcalc/quizx/blob/master/quizx/CHANGELOG.md) for the latest updates.

See the [CONTRIBUTING.md](https://github.com/zxcalc/quizx/blob/master/CONTRIBUTING.md) file for detailed instructions for building the libraries from source and contributing to the project.
Pull requests are welcome!


## License

This project is licensed under Apache License, Version 2.0 ([LICENSE][] or http://www.apache.org/licenses/LICENSE-2.0).
