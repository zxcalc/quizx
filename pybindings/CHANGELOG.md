# Changelog

## Unreleased

## Added
- `add_edge` method to `VecGraph` class.

## [0.3.0](https://github.com/zxcalc/quizx/compare/quizx-py-v0.2.0...quizx-py-v0.3.0) (2025-06-12)


### ⚠ BREAKING CHANGES

* New Drivers in the driver enum
* Added rand as dependency
* The decomposer works different now. I replaced the stack with a recursive approach
* 

### Features

* Advanced drivers and some refactoring ([#149](https://github.com/zxcalc/quizx/issues/149)) ([96e5eb9](https://github.com/zxcalc/quizx/commit/96e5eb9a1b7fceb9aff2c6c2e3f90e0dbcc95c4d))
* Generator for surface code circuits ([#128](https://github.com/zxcalc/quizx/issues/128)) ([4da7be4](https://github.com/zxcalc/quizx/commit/4da7be4c4e0becbb7f11f345eebb425890369f7c))
* Improve Python bindings for the Decomposer ([#151](https://github.com/zxcalc/quizx/issues/151)) ([a53273f](https://github.com/zxcalc/quizx/commit/a53273f5508e0148dbdc031bc66e51c7154a94bb))
* Rework of Decomposer to allow for speperation of disconnected components and exchange the driving decomposition selector ([#147](https://github.com/zxcalc/quizx/issues/147)) ([efe19c5](https://github.com/zxcalc/quizx/commit/efe19c51efb6d7406dd1c696c1464ae65fbaa7c8))

## [0.2.0](https://github.com/zxcalc/quizx/compare/quizx-py-v0.1.1...quizx-py-v0.2.0) (2025-05-13)


### ⚠ BREAKING CHANGES

* 

### Features

* simplified and more comprehensive Python bindings for graphs, simplification, and extraction ([eb8ed58](https://github.com/zxcalc/quizx/commit/eb8ed589cc75c8c1efd67e6dfa39dced379f1611))
* Support for boolean variables on spiders and scalars ([#125](https://github.com/zxcalc/quizx/issues/125)) ([47920e1](https://github.com/zxcalc/quizx/commit/47920e1be2356a3b47f7ef7b861dc5ff8c1413a3))


### Bug Fixes

* make qubit and row into floats (closes [#93](https://github.com/zxcalc/quizx/issues/93)) ([0da3b64](https://github.com/zxcalc/quizx/commit/0da3b64b6a7ad63a690c9c389a04a4e140ec3b55))
* **python:** Implement add_edge in VecGraph ([#102](https://github.com/zxcalc/quizx/issues/102)) ([9718a97](https://github.com/zxcalc/quizx/commit/9718a973ed647cd6ec6961029d42f51fb3a24112))


### Documentation

* **rs:** Add crate documentation ([#122](https://github.com/zxcalc/quizx/issues/122)) ([b7415e7](https://github.com/zxcalc/quizx/commit/b7415e701c37013965944de378ded66a49f8e8e4))


### Miscellaneous Chores

* **py:** release 0.2.0 ([d25e3a9](https://github.com/zxcalc/quizx/commit/d25e3a9f8416a8a7574a422d6f66bb3ee755cf9c))
* release 0.1.2 ([13d5dc7](https://github.com/zxcalc/quizx/commit/13d5dc7b458f6e21f9c9c0ce357abda6b075864a))

## [0.1.1](https://github.com/zxcalc/quizx/compare/quizx-py-v0.1.0...quizx-py-v0.1.1) (2025-02-21)


### Features

* Bumped `pyzx` version to `0.9.0` ([#95](https://github.com/zxcalc/quizx/issues/95))

### Documentation

* Add pypi badges to the readmes ([#84](https://github.com/zxcalc/quizx/issues/84)) ([2d03128](https://github.com/zxcalc/quizx/commit/2d031280d630ebb68b0bc97bd8e71c6629d1319d))

## 0.1.0 (2024-10-30)


This is the initial experimental release of `quizx`.


----

This is the changelog for the `quizx` python library.
For the changelog of the `quizx` rust library, see the separate [`CHANGELOG.md`](https://github.com/zxcalc/quizx/blob/master/CHANGELOG.md) file.
