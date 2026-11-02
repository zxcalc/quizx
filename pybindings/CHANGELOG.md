# Changelog

## [0.4.0](https://github.com/zxcalc/quizx/compare/quizx-py-v0.3.0...quizx-py-v0.4.0) (2026-11-02)

### Breaking Changes

- Rename `quizx.Scalar` to `quizx.Scalar4`. `VecGraph.scalar` now uses PyZX's `Scalar`; use `VecGraph.scalar4` for native scalars without conversion. The scalar getters return copies, so assign modified scalars back to the graph. ([#170](https://github.com/zxcalc/quizx/pull/170))
- `VecGraph.adjoint()` now returns a new graph instead of modifying the original. Use `g = g.adjoint()` to replace a graph with its adjoint. ([#198](https://github.com/zxcalc/quizx/pull/198))

### Features

- Add a dyadic rational type and scalars based on it ([#158](https://github.com/zxcalc/quizx/pull/158))
- Add keyword arguments and properties to bindings ([#169](https://github.com/zxcalc/quizx/pull/169))
- Extend PyZX compatibility to `>=0.9,<0.11` ([#198](https://github.com/zxcalc/quizx/pull/198))
- Add Python bindings for rankwidth functions ([#168](https://github.com/zxcalc/quizx/pull/168))
- Use `bitgauss` for faster GF(2) linear algebra in circuit extraction ([#161](https://github.com/zxcalc/quizx/pull/161))
- Add wheel builds for free-threaded Python 3.14, Windows ARM64, and WebAssembly. WebAssembly wheels are distributed as GitHub release assets. ([#191](https://github.com/zxcalc/quizx/pull/191), [#196](https://github.com/zxcalc/quizx/pull/196))

### Bug Fixes

- Add missing methods to `VecGraph` bindings ([#170](https://github.com/zxcalc/quizx/pull/170))
- Fix conversion from PyZX scalars ([81e198e](https://github.com/zxcalc/quizx/commit/81e198e242aa161ec03ef9a6ba5d87fd1d5e9005))
- Fix cut-rank caching ([#173](https://github.com/zxcalc/quizx/pull/173))
- Remove stale edges from the rank cache ([7ab591f](https://github.com/zxcalc/quizx/commit/7ab591f4802f581a57c024126465fb086d2ce433))
- Fix rank-cache invalidation when swapping leaves or moving subtrees ([1acc52a](https://github.com/zxcalc/quizx/commit/1acc52a17d1b0dd76ab49fbbad6e81d6602f49ea))
- Preserve inputs, outputs, and scalars in `VecGraph.copy()` ([#199](https://github.com/zxcalc/quizx/pull/199))
- Avoid adding two edges for a self-loop in `VecGraph` ([#187](https://github.com/zxcalc/quizx/pull/187))
- Parse Qiskit phase gates in custom QASM definitions passed to `quizx.qasm()` ([#200](https://github.com/zxcalc/quizx/pull/200))
- Fix mantissa overflows when converting and formatting scalars ([#204](https://github.com/zxcalc/quizx/pull/204))

### Documentation

- Correct repository metadata and contribution links ([#197](https://github.com/zxcalc/quizx/pull/197))

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
For the changelog of the `quizx` rust library, see the separate [`CHANGELOG.md`](https://github.com/zxcalc/quizx/blob/master/quizx/CHANGELOG.md) file.
