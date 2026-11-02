# Changelog

This is the changelog for the `quizx` rust library.
For the changelog of the `quizx` python library, see the separate [`CHANGELOG.md`](https://github.com/zxcalc/quizx/blob/master/pybindings/CHANGELOG.md) file.


## [0.4.0](https://github.com/zxcalc/quizx/compare/quizx@v0.3.0...quizx@v0.4.0) - 2026-11-01

### Breaking Changes

- Graphs now use `Scalar4` with four `Dyadic` coefficients instead of `FScalar`. Replace `FScalar::dyadic(pow, coeffs)` with `Scalar4::new(coeffs, pow)`, and use `TryFrom`/`TryInto` to convert to `Complex<f64>`. ([#158](https://github.com/zxcalc/quizx/pull/158))
- Use `Tensor4`/`to_tensor4()` for `Scalar4` tensors. `TensorF`/`to_tensorf()` now represent `Complex<f64>` tensors, replacing `Tensor64`/`to_tensor64()`. Custom `TensorElem` implementations must provide `TryFrom<Scalar4>` with a `Debug` error type instead of `From<FScalar>`. ([#158](https://github.com/zxcalc/quizx/pull/158))
- Upgrade `rand` to 0.9. Downstream code passing random number generators into QuiZX APIs must use the corresponding `rand` 0.9 types and traits. ([#166](https://github.com/zxcalc/quizx/pull/166))

### Bug Fixes

- Fix a typo in the README
- Fix cut-rank caching ([#173](https://github.com/zxcalc/quizx/pull/173))
- Remove stale edges from the rank cache
- Fix another rankwidth caching bug
- Preserve boundaries and scalars when copying graphs ([#199](https://github.com/zxcalc/quizx/pull/199))
- Parse Qiskit phase gates in custom QASM definitions ([#200](https://github.com/zxcalc/quizx/pull/200))
- Avoid adding two edges for a self-loop ([#187](https://github.com/zxcalc/quizx/pull/187))
- Support ancillae and post-selection in circuit tensors ([#203](https://github.com/zxcalc/quizx/pull/203))
- Fix mantissa overflows ([#204](https://github.com/zxcalc/quizx/pull/204))

### New Features

- [**breaking**] Add a dyadic rational type and scalars based on it ([#158](https://github.com/zxcalc/quizx/pull/158))
- Use `bitgauss` for fast GF(2) linear algebra ([#161](https://github.com/zxcalc/quizx/pull/161))
- [**breaking**] Add low-rankwidth decompositions of ZX diagrams ([#166](https://github.com/zxcalc/quizx/pull/166))
- [**breaking**] Add Python bindings for rankwidth functions ([#168](https://github.com/zxcalc/quizx/pull/168))
- Add keyword arguments and properties to Python bindings ([#169](https://github.com/zxcalc/quizx/pull/169))
- Extend PyZX compatibility to `>=0.9,<0.11` ([#198](https://github.com/zxcalc/quizx/pull/198))
- Implement approximate equality of tensors with floating-point scalars ([#137](https://github.com/zxcalc/quizx/issues/137), [#182](https://github.com/zxcalc/quizx/pull/182))
- Implement `Circuit::optimize` (basic gate-level optimisation) ([#193](https://github.com/zxcalc/quizx/pull/193))
- Implement the Z–X bialgebra rewrite rule ([#189](https://github.com/zxcalc/quizx/pull/189))

### Testing

- Add unit tests for `gate`, `phase`, `util`, `random_graph`, `annealer`, `cli`, and `json` ([#188](https://github.com/zxcalc/quizx/pull/188))

## [0.3.0](https://github.com/zxcalc/quizx/compare/quizx@v0.2.0...quizx@v0.3.0) - 2025-06-12

### Bug Fixes

- bug in pack() not setting boundary names

### New Features

- added method to pack VecGraph vertices, removing holes ([#126](https://github.com/zxcalc/quizx/pull/126))
- [**breaking**] Generator for surface code circuits ([#128](https://github.com/zxcalc/quizx/pull/128))
- [**breaking**] converting a circuit to a graph while performing Clifford simplifications ([#133](https://github.com/zxcalc/quizx/pull/133))
- [**breaking**] better equality checker ([#136](https://github.com/zxcalc/quizx/pull/136))
- Add command line interface ([#139](https://github.com/zxcalc/quizx/pull/139))
- converting Clifford circuits to graphs has the option to maintain AP form ([#141](https://github.com/zxcalc/quizx/pull/141))
- [**breaking**] Rework of Decomposer to allow for speperation of disconnected components and exchange the driving decomposition selector ([#147](https://github.com/zxcalc/quizx/pull/147))
- [**breaking**] Add spider-cutting decomposition ([#143](https://github.com/zxcalc/quizx/pull/143))
- detection web module ([#142](https://github.com/zxcalc/quizx/pull/142))
- [**breaking**] Advanced drivers and some refactoring ([#149](https://github.com/zxcalc/quizx/pull/149))

### Testing

- add benchmarks ([#148](https://github.com/zxcalc/quizx/pull/148))
- Make CLI tests integration tests instead of unit tests ([#153](https://github.com/zxcalc/quizx/pull/153))

## [0.2.0](https://github.com/zxcalc/quizx/compare/quizx@v0.1.0...quizx@v0.2.0) - 2025-05-12

### Documentation

- *(rs)* Add crate documentation ([#122](https://github.com/zxcalc/quizx/pull/122))

### New Features

- Simplification-Based Equality Checker ([#110](https://github.com/zxcalc/quizx/pull/110))
- [**breaking**] Support for boolean variables on spiders and scalars ([#125](https://github.com/zxcalc/quizx/pull/125))

## [0.1.0](https://github.com/zxcalc/quizx/releases/tag/quizx@v0.1.0) - 2025-05-04

First crate release
