// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(clippy::useless_conversion)]

pub mod scalar;

use crate::scalar::Scalar;

use num::Rational64;
use pyo3::exceptions::PyNotImplementedError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use ::quizx::extract::ToCircuit;
use ::quizx::graph::*;
use ::quizx::phase::Phase;

type E = (V,V);

#[pymodule]
fn quizx(m: &Bound<'_, PyModule>) -> PyResult<()> {
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
fn extract_circuit(g: &mut VecGraph) -> Circuit {
    Circuit {
        c: g.g.to_circuit().unwrap(),
        s: None,
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

/// Wrapper for quizx::vec_graph::Graph
#[pyclass]
struct VecGraph {
    pub g: ::quizx::vec_graph::Graph,
}

#[pymethods]
impl VecGraph {
    #[new]
    fn new() -> VecGraph {
        VecGraph {
            g: ::quizx::vec_graph::Graph::new(),
        }
    }

    // TODO: should return a PyZX-compatible scalar

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

    fn clone(&self) -> VecGraph {
        VecGraph { g: self.g.clone() }
    }

    fn inputs(&self) -> Vec<V> {
        self.g.inputs().clone()
    }

    fn set_inputs(&mut self, inputs: Vec<V>) {
        self.g.set_inputs(inputs)
    }

    fn outputs(&self) -> Vec<V> {
        self.g.outputs().clone()
    }

    fn set_outputs(&mut self, outputs: Vec<V>) {
        self.g.set_outputs(outputs)
    }

    fn add_vertices(&mut self, amount: usize) -> Vec<V> {
        Vec::from_iter(std::iter::repeat_n(self.g.add_vertex(VType::B), amount))
    }

    fn add_vertex_indexed(&mut self, v: V) -> PyResult<()> {
        match self.g.add_named_vertex_with_data(v, VData::default()) {
            Ok(_) => Ok(()),
            Err(_) => Err(PyValueError::new_err("Vertex already exists")),
        }
    }

    #[pyo3(signature = (edge_pair, edgetype=1))]
    fn add_edge(&mut self, edge_pair: (usize, usize), edgetype: u8) {
        let et = match edgetype {
            2 => EType::H,
            _ => EType::N,
        };
        self.g.add_edge_smart(edge_pair.0, edge_pair.1, et)
    }

    fn remove_vertices(&mut self, vertices: PyObject) -> PyResult<()> {
        Python::with_gil(|py| -> PyResult<()> {
            for vo in vertices.bind(py).try_iter()? {
                let v = vo?.extract::<V>()?;
                self.g.remove_vertex(v);
            }
            Ok(())
        })
    }

    fn remove_edges(&mut self, edges: PyObject) -> PyResult<()> {
        Python::with_gil(|py| -> PyResult<()> {
            for eo in edges.bind(py).try_iter()? {
                let e = eo?.extract::<E>()?;
                self.g.remove_edge(e.0, e.1);
            }
            Ok(())
        })
    }

    fn num_vertices(&self) -> usize {
        self.g.num_vertices()
    }

    fn num_edges(&self) -> usize {
        self.g.num_edges()
    }

    fn vertices(&self) -> Vec<V> {
        Vec::from_iter(self.g.vertices())
    }

    #[pyo3(signature = (s=None, t=None))]
    fn edges(&self, s: Option<V>, t: Option<V>) -> Vec<E> {
        // TODO: handle s and t args
        Vec::from_iter(self.g.edges().map(|(s,t,_)| (s,t)))
    }

    fn edge_st(&self, edge: E) -> (V, V) {
        edge
    }

    fn incident_edges(&self, vertex: V) -> Vec<E> {
        Vec::from_iter(self.g.neighbors(vertex).map(|w|
            if vertex < w { (vertex, w) } else { (w, vertex) }))
    }

    fn edge_type(&self, e: E) -> u8 {
        match self.g.edge_type_opt(e.0, e.1) {
            Some(EType::N) => 1,
            Some(EType::H) => 2,
            Some(EType::Wio) => 3,
            None => 0,
        }
    }

    fn set_edge_type(&mut self, e: E, t: u8) {
        let et = match t {
            2 => EType::H,
            3 => EType::Wio,
            _ => EType::N,
        };
        self.g.set_edge_type(e.0, e.1, et);
    }

    #[pyo3(name="type")]
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

    fn set_type(&mut self, vertex: V, t: u8) {
        let ty = match t {
            1 => VType::Z,
            2 => VType::X,
            3 => VType::H,
            4 => VType::WInput,
            5 => VType::WOutput,
            6 => VType::ZBox,
            _ => VType::B,
        };
        self.g.set_vertex_type(vertex, ty);
    }

    fn phase(&self, v: usize) -> (i64, i64) {
        // TODO: should return Fraction
        let p = self.g.phase(v).to_rational();
        (*p.numer(), *p.denom())
    }

    fn set_phase(&mut self, v: usize, phase: (i64, i64)) {
        // TODO: should accept Fraction
        self.g.set_phase(v, Rational64::new(phase.0, phase.1));
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

    fn clear_vdata(&self, vertex: V) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn vdata_keys(&self, vertex: V) -> PyResult<Vec<String>> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn vdata(&self, vertex: V, key: String, default: PyObject) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn set_vdata(&self, vertex: V, key: String, val: PyObject) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    // TODO: fix signatures below
    fn is_ground(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn grounds(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn set_ground(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn is_hybrid(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn multigraph(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }


    fn phases(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn types(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn qubits(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn rows(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn depth(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn edge(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn connected(&self, s: V, t: V) -> bool {
        self.g.connected(s, t)
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

    fn add_edges(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn remove_vertex(&mut self, v: V) {
        self.g.remove_vertex(v)
    }

    fn remove_edge(&mut self, e: (usize, usize)) {
        self.g.remove_edge(e.0, e.1)
    }

    fn add_to_phase(&mut self, v: usize, phase: (i64, i64)) {
        // TODO: should accept Fraction
        self.g.add_to_phase(v, Rational64::new(phase.0, phase.1));
    }


    fn num_inputs(&self) -> usize {
        self.g.inputs().len()
    }

    fn num_outputs(&self) -> usize {
        self.g.outputs().len()
    }

    fn set_position(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn neighbors(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn vertex_degree(&self, v: V) -> usize {
        self.g.degree(v)
    }

    fn edge_s(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn edge_t(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn vertex_set(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn edge_set(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn stats(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn copy(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn adjoint(&mut self) {
        self.g.adjoint()
    }

    fn map_qubits(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn compose(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn tensor(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn __iadd__(&self, other: PyObject) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn __add__(&self, other: PyObject) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn __mul__(&self, other: PyObject) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn __matmul__(&self, other: PyObject) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn merge(&self, other: PyObject) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn subgraph_from_vertices(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn apply_state(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn apply_effect(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn to_tensor(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn to_matrix(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn to_dict(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn to_json(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn to_graphml(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn to_tikz(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    // fn from_json(cls, js:Union[str,Dict[str,Any]]) -> VecGraph: ...
    // fn from_tikz(cls, tikz: str, warn_overlap:bool= True, fuse_overlap:bool = True, ignore_nonzx:bool = False) -> VecGraph: ...

    fn is_id(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn pack_circuit_rows(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn qubit_count(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn auto_detect_io(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn normalize(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn translate(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn add_edge_table(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn set_phase_master(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn update_phase_index(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn fuse_phases(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn phase_negate(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn vertex_from_phase_index(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn remove_isolated_vertices(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn vdata_dict(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn set_vdata_dict(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn is_well_formed(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn get_auto_simplify(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn set_auto_simplify(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
    }

    fn is_phase_gadget(&self) -> PyResult<()> {
        return Err(PyNotImplementedError::new_err("Not implemented on backend: quizx_vec"));
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
