// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(clippy::useless_conversion)]

pub mod circuit;
pub mod decompose;
pub mod rankwidth;
pub mod scalar;
pub mod util;
pub mod vec_graph;

use crate::circuit::to_pyzx_circuit;
use crate::decompose::{PyDecomposer, SimpFunc};
use crate::rankwidth::PyDecompTree;
use crate::scalar::PyScalar;
use crate::vec_graph::PyVecGraph;

use ::quizx::extract::ExtractError;
use ::quizx::extract::ToCircuit;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pymodule]
fn quizx(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(interior_clifford_simp, m)?)?;
    m.add_function(wrap_pyfunction!(clifford_simp, m)?)?;
    m.add_function(wrap_pyfunction!(fuse_gadgets, m)?)?;
    m.add_function(wrap_pyfunction!(full_simp, m)?)?;
    m.add_function(wrap_pyfunction!(extract_circuit, m)?)?;
    m.add_function(wrap_pyfunction!(qasm, m)?)?;
    m.add_function(wrap_pyfunction!(new_graph, m)?)?;
    m.add_class::<PyVecGraph>()?;
    m.add_class::<PyDecomposer>()?;
    m.add_class::<SimpFunc>()?;
    m.add_class::<PyScalar>()?;
    m.add_class::<PyDecompTree>()?;
    Ok(())
}

#[pyfunction(name = "Graph")]
fn new_graph() -> PyVecGraph {
    PyVecGraph::new()
}

#[pyfunction]
fn qasm(source: &str) -> PyResult<PyVecGraph> {
    let c = ::quizx::circuit::Circuit::from_qasm(source).map_err(PyValueError::new_err)?;
    Ok(PyVecGraph { g: c.to_graph() })
}

#[pyfunction]
fn interior_clifford_simp(g: &mut PyVecGraph) {
    ::quizx::simplify::interior_clifford_simp(&mut g.g);
}

#[pyfunction]
fn clifford_simp(g: &mut PyVecGraph) {
    ::quizx::simplify::clifford_simp(&mut g.g);
}

#[pyfunction]
fn fuse_gadgets(g: &mut PyVecGraph) {
    ::quizx::simplify::fuse_gadgets(&mut g.g);
}

#[pyfunction]
fn full_simp(g: &mut PyVecGraph) {
    ::quizx::simplify::full_simp(&mut g.g);
}

#[pyfunction]
fn extract_circuit(py: Python<'_>, g: &mut PyVecGraph) -> PyResult<PyObject> {
    match g.g.to_circuit() {
        Ok(c) => to_pyzx_circuit(py, c),
        Err(ExtractError(s, _, _)) => Err(PyValueError::new_err(s)),
    }
}
