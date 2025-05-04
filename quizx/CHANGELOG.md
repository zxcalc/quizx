# Changelog

This is the changelog for the `quizx` rust library.
For the changelog of the `quizx` python library, see the separate [`CHANGELOG.md`](https://github.com/zxcalc/quizx/blob/master/pybindings/CHANGELOG.md) file.


## [0.0.0](https://github.com/zxcalc/quizx/releases/tag/quizx@v0.0.0) - 2025-05-04

### Bug Fixes

- fixed some bugs and tested decomp up to size 10
- fixed seed for candidates
- fixed bug with gadget fusion
- fixed toffoli bug
- fixed bug in test
- fixed overflow problem in tensor contraction
- fixed bad scalar in fuse_gadgets
- Lints and warnings ([#16](https://github.com/zxcalc/quizx/pull/16))
- Move rust bindings to a `quizx._quizx` submodule ([#20](https://github.com/zxcalc/quizx/pull/20))
- gadget fusion ([#40](https://github.com/zxcalc/quizx/pull/40))
- sdist build config ([#73](https://github.com/zxcalc/quizx/pull/73))
- Fix bug in JsonScalar::from_scalar ([#74](https://github.com/zxcalc/quizx/pull/74))
- make qubit and row into floats (closes #93)
- pick qubit number for unfuse_gadget to match pyzx convention
- incorrect matching of cat states on ZX diagrams with boundary (fixes #64)
- satisfy autoformatter
- bug in `Decomposer::push_cat_decomp` ([#107](https://github.com/zxcalc/quizx/pull/107))
- stop skipping json serialisation of graph scalars ([#108](https://github.com/zxcalc/quizx/pull/108))

### Documentation

- Update links to the renamed organization ([#45](https://github.com/zxcalc/quizx/pull/45))

### New Features

- Rust decoder for `.qgraph` files (pyzx/quantomatic graph format) ([#22](https://github.com/zxcalc/quizx/pull/22))
- `Phase` struct ([#31](https://github.com/zxcalc/quizx/pull/31))
- Scalar serialization and pyzx interop ([#33](https://github.com/zxcalc/quizx/pull/33))
- Small update to the Python API ([#76](https://github.com/zxcalc/quizx/pull/76))

### Performance

- minor improvements to cat decomp

### Refactor

- simplified and updated Python/PyZX bindings ([#106](https://github.com/zxcalc/quizx/pull/106))
- simplified scalar implementation ([#109](https://github.com/zxcalc/quizx/pull/109))
- use DecompFunc enum rather than boolean flag

### Testing

- tested BSS decomp and it works
# Changelog

This is the changelog for the `quizx` rust library.
For the changelog of the `quizx` python library, see the separate [`CHANGELOG.md`](https://github.com/zxcalc/quizx/blob/master/pybindings/CHANGELOG.md) file.
