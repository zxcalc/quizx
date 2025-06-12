# Changelog

This is the changelog for the `quizx` rust library.
For the changelog of the `quizx` python library, see the separate [`CHANGELOG.md`](https://github.com/zxcalc/quizx/blob/master/pybindings/CHANGELOG.md) file.


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
