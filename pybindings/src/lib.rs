// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(clippy::useless_conversion)]

pub mod circuit;
pub mod scalar;
pub mod vec_graph;

use crate::circuit::to_pyzx_circuit;
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
    m.add_class::<VecGraph>()?;
    // m.add_class::<Circuit>()?;
    // m.add_class::<CircuitStats>()?;
    m.add_class::<Decomposer>()?;
    m.add_class::<Scalar>()?;
    Ok(())
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

#[pyclass]
struct CircuitStats {
    s: ::quizx::circuit::CircuitStats,
}

/// A (mostly) opaque wrapper for quizx circuits
#[pyclass]
struct Circuit {
    c: ::quizx::circuit::Circuit,
    s: Option<::quizx::circuit::CircuitStats>,
}

#[pymethods]
impl Circuit {
    #[staticmethod]
    fn from_qasm(qasm: String) -> Circuit {
        Circuit {
            c: ::quizx::circuit::Circuit::from_qasm(&qasm).unwrap(),
            s: None,
        }
    }

    #[staticmethod]
    fn load(file: String) -> Circuit {
        Circuit {
            c: ::quizx::circuit::Circuit::from_file(&file).unwrap(),
            s: None,
        }
    }

    fn to_qasm(&self) -> String {
        self.c.to_qasm()
    }
    fn to_graph(&self) -> VecGraph {
        VecGraph {
            g: self.c.to_graph(),
        }
    }

    fn num_gates(&self) -> usize {
        self.c.num_gates()
    }
    fn stats(&mut self) -> CircuitStats {
        // generate stats the first time this method is called
        if self.s.is_none() {
            self.s = Some(self.c.stats());
        }
        CircuitStats { s: self.s.unwrap() }
    }
}

#[pymethods]
impl CircuitStats {
    fn qubits(&self) -> usize {
        self.s.qubits
    }
    fn total(&self) -> usize {
        self.s.total
    }
    fn oneq(&self) -> usize {
        self.s.oneq
    }
    fn twoq(&self) -> usize {
        self.s.twoq
    }
    fn moreq(&self) -> usize {
        self.s.moreq
    }
    fn cliff(&self) -> usize {
        self.s.cliff
    }
    fn non_cliff(&self) -> usize {
        self.s.non_cliff
    }
    #[allow(clippy::inherent_to_string)]
    fn to_string(&self) -> String {
        self.s.to_string()
    }
}

#[pyclass]
struct Decomposer {
    d: ::quizx::decompose::Decomposer<::quizx::vec_graph::Graph>,
}

#[pymethods]
impl Decomposer {
    #[staticmethod]
    fn empty() -> Decomposer {
        Decomposer {
            d: ::quizx::decompose::Decomposer::empty(),
        }
    }

    #[new]
    fn new(g: &VecGraph) -> Decomposer {
        Decomposer {
            d: ::quizx::decompose::Decomposer::new(&g.g),
        }
    }

    #[getter]
    fn get_scalar(&self) -> Scalar {
        self.d.scalar.into()
    }

    #[setter]
    fn set_scalar(&mut self, scalar: Scalar) {
        self.d.scalar = scalar.into();
    }

    fn graphs(&self) -> PyResult<Vec<VecGraph>> {
        let mut gs = vec![];
        for (_a, g) in &self.d.stack {
            gs.push(VecGraph { g: g.clone() });
        }
        Ok(gs)
    }

    fn done(&self) -> PyResult<Vec<VecGraph>> {
        let mut gs = vec![];
        for g in &self.d.done {
            gs.push(VecGraph { g: g.clone() });
        }
        Ok(gs)
    }

    fn save(&mut self, b: bool) {
        self.d.save(b);
    }

    fn apply_optimizations(&mut self, b: bool) {
        if b {
            self.d.with_simp(::quizx::decompose::SimpFunc::FullSimp);
        } else {
            self.d.with_simp(::quizx::decompose::SimpFunc::NoSimp);
        }
    }

    fn max_terms(&self) -> f64 {
        self.d.max_terms()
    }
    fn decomp_top(&mut self) {
        self.d.decomp_top();
    }
    fn decomp_all(&mut self) {
        self.d.decomp_all();
    }
    fn decomp_until_depth(&mut self, depth: usize) {
        self.d.decomp_until_depth(depth);
    }
    fn decomp_parallel(&mut self, depth: usize) {
        self.d = self.d.clone().decomp_parallel(depth);
    }
    fn use_cats(&mut self, b: bool) {
        self.d.use_cats(b);
    }
    fn get_nterms(&self) -> usize {
        self.d.nterms
    }
}
