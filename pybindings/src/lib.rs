// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(clippy::useless_conversion)]

pub mod scalar;

use crate::scalar::Scalar;

use num::Rational64;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use quizx::extract::ToCircuit;
use quizx::graph::*;
use quizx::phase::Phase;

#[pymodule]
fn _quizx(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dummy, m)?)?;
    m.add_function(wrap_pyfunction!(interior_clifford_simp, m)?)?;
    m.add_function(wrap_pyfunction!(clifford_simp, m)?)?;
    m.add_function(wrap_pyfunction!(fuse_gadgets, m)?)?;
    m.add_function(wrap_pyfunction!(full_simp, m)?)?;
    m.add_function(wrap_pyfunction!(extract_circuit, m)?)?;
    m.add_class::<VecGraph>()?;
    m.add_class::<Circuit>()?;
    m.add_class::<CircuitStats>()?;
    m.add_class::<Decomposer>()?;
    m.add_class::<Scalar>()?;
    Ok(())
}

#[pyfunction]
fn dummy(a: i64) -> String {
    format!("FOO! {}", a)
}

#[pyfunction]
fn interior_clifford_simp(g: &mut VecGraph) {
    quizx::simplify::interior_clifford_simp(&mut g.g);
}

#[pyfunction]
fn clifford_simp(g: &mut VecGraph) {
    quizx::simplify::clifford_simp(&mut g.g);
}

#[pyfunction]
fn fuse_gadgets(g: &mut VecGraph) {
    quizx::simplify::fuse_gadgets(&mut g.g);
}

#[pyfunction]
fn full_simp(g: &mut VecGraph) {
    quizx::simplify::full_simp(&mut g.g);
}

#[pyfunction]
fn extract_circuit(g: &mut VecGraph) -> Circuit {
    Circuit {
        c: g.g.to_circuit().unwrap(),
        s: None,
    }
}

#[pyclass]
struct CircuitStats {
    s: quizx::circuit::CircuitStats,
}

/// A (mostly) opaque wrapper for quizx circuits
#[pyclass]
struct Circuit {
    c: quizx::circuit::Circuit,
    s: Option<quizx::circuit::CircuitStats>,
}

#[pymethods]
impl Circuit {
    #[staticmethod]
    fn from_qasm(qasm: String) -> Circuit {
        Circuit {
            c: quizx::circuit::Circuit::from_qasm(&qasm).unwrap(),
            s: None,
        }
    }

    #[staticmethod]
    fn load(file: String) -> Circuit {
        Circuit {
            c: quizx::circuit::Circuit::from_file(&file).unwrap(),
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

/// Wrapper for quizx::vec_graph::Graph
#[pyclass]
struct VecGraph {
    pub g: quizx::vec_graph::Graph,
}

#[pymethods]
impl VecGraph {
    #[new]
    fn new() -> VecGraph {
        VecGraph {
            g: quizx::vec_graph::Graph::new(),
        }
    }

    fn vindex(&self) -> usize {
        self.g.vindex()
    }
    fn neighbor_at(&self, v: usize, n: usize) -> usize {
        self.g.neighbor_at(v, n)
    }
    fn num_vertices(&self) -> usize {
        self.g.num_vertices()
    }
    fn num_edges(&self) -> usize {
        self.g.num_edges()
    }
    fn add_vertex(&mut self, ty_num: u8, qubit: f64, row: f64, phase: (i64, i64)) -> usize {
        let ty = match ty_num {
            1 => VType::Z,
            2 => VType::X,
            3 => VType::H,
            _ => VType::B,
        };
        let phase = Phase::new((phase.0, phase.1));
        self.g.add_vertex_with_data(VData {
            ty,
            phase,
            qubit,
            row,
        })
    }

    fn contains_vertex(&self, v: usize) -> bool {
        self.g.contains_vertex(v)
    }

    fn add_edge(&mut self, e: (usize, usize), et_num: u8) {
        let et = match et_num {
            2 => EType::H,
            _ => EType::N,
        };
        self.g.add_edge_with_type(e.0, e.1, et)
    }

    fn add_edge_smart(&mut self, e: (usize, usize), et_num: u8) {
        let et = match et_num {
            1 => EType::H,
            _ => EType::N,
        };
        self.g.add_edge_smart(e.0, e.1, et)
    }

    fn remove_vertex(&mut self, v: usize) {
        self.g.remove_vertex(v)
    }
    fn remove_edge(&mut self, e: (usize, usize)) {
        self.g.remove_edge(e.0, e.1)
    }
    fn degree(&self, v: usize) -> usize {
        self.g.degree(v)
    }
    fn connected(&self, s: usize, t: usize) -> bool {
        self.g.connected(s, t)
    }

    fn vertex_type(&self, v: usize) -> u8 {
        match self.g.vertex_type(v) {
            VType::B => 0,
            VType::Z => 1,
            VType::X => 2,
            VType::H => 3,
            VType::WInput => 4,
            VType::WOutput => 5,
            VType::ZBox => 6,
        }
    }

    fn set_vertex_type(&mut self, v: usize, ty_num: u8) {
        let ty = match ty_num {
            1 => VType::Z,
            2 => VType::X,
            3 => VType::H,
            4 => VType::WInput,
            5 => VType::WOutput,
            6 => VType::ZBox,
            _ => VType::B,
        };
        self.g.set_vertex_type(v, ty);
    }

    fn edge_type(&self, e: (usize, usize)) -> u8 {
        match self.g.edge_type_opt(e.0, e.1) {
            Some(EType::N) => 1,
            Some(EType::H) => 2,
            Some(EType::Wio) => 3,
            None => 0,
        }
    }

    fn set_edge_type(&mut self, e: (usize, usize), et_num: u8) {
        let et = match et_num {
            2 => EType::H,
            3 => EType::Wio,
            _ => EType::N,
        };
        self.g.set_edge_type(e.0, e.1, et);
    }

    fn phase(&self, v: usize) -> (i64, i64) {
        let p = self.g.phase(v).to_rational();
        (*p.numer(), *p.denom())
    }

    fn set_phase(&mut self, v: usize, phase: (i64, i64)) {
        self.g.set_phase(v, Rational64::new(phase.0, phase.1));
    }

    fn add_to_phase(&mut self, v: usize, phase: (i64, i64)) {
        self.g.add_to_phase(v, Rational64::new(phase.0, phase.1));
    }

    fn qubit(&mut self, v: usize) -> f64 {
        self.g.qubit(v)
    }
    fn set_qubit(&mut self, v: usize, q: f64) {
        self.g.set_qubit(v, q);
    }
    fn row(&mut self, v: usize) -> f64 {
        self.g.row(v)
    }
    fn set_row(&mut self, v: usize, r: f64) {
        self.g.set_row(v, r);
    }

    fn inputs(&self) -> Vec<V> {
        self.g.inputs().clone()
    }
    fn num_inputs(&self) -> usize {
        self.g.inputs().len()
    }
    fn set_inputs(&mut self, inputs: Vec<V>) {
        self.g.set_inputs(inputs)
    }
    fn outputs(&self) -> Vec<V> {
        self.g.outputs().clone()
    }
    fn num_outputs(&self) -> usize {
        self.g.outputs().len()
    }
    fn set_outputs(&mut self, outputs: Vec<V>) {
        self.g.set_outputs(outputs)
    }

    /// Returns the graph scalar.
    #[getter]
    fn get_scalar(&self) -> Scalar {
        self.g.scalar().clone().into()
    }

    /// Sets the graph scalar.
    #[setter]
    fn set_scalar(&mut self, scalar: Scalar) {
        *self.g.scalar_mut() = scalar.into();
    }

    fn adjoint(&mut self) {
        self.g.adjoint()
    }

    fn plug(&mut self, other: &VecGraph) {
        self.g.plug(&other.g);
    }

    fn clone(&self) -> VecGraph {
        VecGraph { g: self.g.clone() }
    }
}

#[pyclass]
struct Decomposer {
    d: quizx::decompose::Decomposer<quizx::vec_graph::Graph>,
}

#[pymethods]
impl Decomposer {
    #[staticmethod]
    fn empty() -> Decomposer {
        Decomposer {
            d: quizx::decompose::Decomposer::empty(),
        }
    }

    #[new]
    fn new(g: &VecGraph) -> Decomposer {
        Decomposer {
            d: quizx::decompose::Decomposer::new(&g.g),
        }
    }

    #[getter]
    fn get_scalar(&self) -> Scalar {
        self.d.scalar.clone().into()
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
            self.d.with_simp(quizx::decompose::SimpFunc::FullSimp);
        } else {
            self.d.with_simp(quizx::decompose::SimpFunc::NoSimp);
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
