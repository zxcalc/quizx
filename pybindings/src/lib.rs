// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(clippy::useless_conversion)]

pub mod circuit;
pub mod decompose;
pub mod scalar;
pub mod util;
pub mod vec_graph;

use crate::circuit::to_pyzx_circuit;
use crate::decompose::Decomposer;
use crate::scalar::Scalar;
use crate::vec_graph::VecGraph;

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
    m.add_class::<VecGraph>()?;
    m.add_class::<Decomposer>()?;
    m.add_class::<Scalar>()?;
    Ok(())
}

#[pyfunction]
fn qasm(source: &str) -> PyResult<VecGraph> {
    let c = ::quizx::circuit::Circuit::from_qasm(source).map_err(PyValueError::new_err)?;
    Ok(VecGraph { g: c.to_graph() })
}

#[pyfunction]
fn interior_clifford_simp(g: &mut VecGraph) {
    ::quizx::simplify::interior_clifford_simp(&mut g.g);
}

#[pyfunction]
fn clifford_simp(g: &mut VecGraph) {
    ::quizx::simplify::clifford_simp(&mut g.g);
}

#[pyfunction]
fn fuse_gadgets(g: &mut VecGraph) {
    ::quizx::simplify::fuse_gadgets(&mut g.g);
}

#[pyfunction]
fn full_simp(g: &mut VecGraph) {
    ::quizx::simplify::full_simp(&mut g.g);
}

#[pyfunction]
fn extract_circuit(py: Python<'_>, g: &mut VecGraph) -> PyResult<PyObject> {
    match g.g.to_circuit() {
        Ok(c) => to_pyzx_circuit(py, c),
        Err(ExtractError(s, _, _)) => Err(PyValueError::new_err(s)),
    }
}
