// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(clippy::useless_conversion)]

pub mod scalar;

use std::collections::HashMap;
use std::collections::HashSet;

use crate::scalar::Scalar;

use ::quizx::extract::ToCircuit;
use ::quizx::graph::*;
use ::quizx::phase::*;
use num::Rational64;
use num::Zero;
use pyo3::exceptions::*;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

type E = (V, V);

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
        match self.g.add_named_vertex_with_data(v, VData::empty()) {
            Ok(_) => Ok(()),
            Err(_) => Err(PyValueError::new_err("Vertex already exists")),
        }
    }

    #[pyo3(signature = (edge_pair, edgetype=1))]
    fn add_edge(&mut self, edge_pair: (usize, usize), edgetype: u8) {
        let et = match edgetype {
            2 => EType::H,
            3 => EType::Wio,
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
        match (s, t) {
            (Some(s), Some(t)) => {
                if self.g.connected(s, t) {
                    vec![self.edge(s, t)]
                } else {
                    vec![]
                }
            }
            (Some(s), None) => self.incident_edges(s),
            _ => Vec::from_iter(self.g.edges().map(|(s, t, _)| (s, t))),
        }
    }

    fn edge_st(&self, edge: E) -> (V, V) {
        edge
    }

    fn incident_edges(&self, vertex: V) -> Vec<E> {
        Vec::from_iter(self.g.neighbors(vertex).map(|w| {
            if vertex < w {
                (vertex, w)
            } else {
                (w, vertex)
            }
        }))
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

    #[pyo3(name = "type")]
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

    fn phase(&self, v: usize) -> PyResult<PyObject> {
        let p = self.g.phase(v).to_rational();
        // (*p.numer(), *p.denom())
        Python::with_gil(|py| -> PyResult<PyObject> {
            let fractions = PyModule::import(py, "fractions")?;
            let f = fractions.getattr("Fraction")?;
            Ok(f.call((*p.numer(), *p.denom()), None)?.unbind())
        })
    }

    fn set_phase(&mut self, v: usize, phase: PyObject) -> PyResult<()> {
        Python::with_gil(|py| -> PyResult<()> {
            let phase1 = pyfraction_to_phase(phase.bind(py))?;
            self.g.set_phase(v, phase1);
            Ok(())
        })
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
        let _vertex = vertex;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn vdata_keys(&self, vertex: V) -> PyResult<Vec<String>> {
        let _vertex = vertex;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn vdata(&self, vertex: V, key: String, default: PyObject) -> PyResult<()> {
        let _vertex = vertex;
        let _key = key;
        let _default = default;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn set_vdata(&self, vertex: V, key: String, val: PyObject) -> PyResult<()> {
        let _vertex = vertex;
        let _key = key;
        let _val = val;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    // TODO: fix signatures below
    #[pyo3(signature = (vertex))]
    fn is_ground(&self, vertex: V) -> bool {
        let _vertex = vertex;
        false
    }

    fn grounds(&self) -> Vec<V> {
        vec![]
    }

    #[pyo3(signature = (vertex, flag=true))]
    fn set_ground(&self, vertex: V, flag: bool) -> PyResult<()> {
        let _vertex = vertex;
        let _flag = flag;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn is_hybrid(&self) -> bool {
        false
    }

    fn multigraph(&self) -> bool {
        false
    }

    fn phases(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn types(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn qubits(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn rows(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn depth(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn edge(&self, s: V, t: V) -> E {
        if s < t {
            (s, t)
        } else {
            (t, s)
        }
    }

    fn connected(&self, s: V, t: V) -> bool {
        self.g.connected(s, t)
    }

    #[pyo3(signature = (ty=0, qubit=-1.0, row=-1.0, phase=None, ground=false))]
    fn add_vertex(
        &mut self,
        ty: u8,
        qubit: f64,
        row: f64,
        phase: Option<PyObject>,
        ground: bool,
    ) -> PyResult<usize> {
        let _ground = ground;
        // TODO: should accept Fraction for phase
        let ty1 = match ty {
            1 => VType::Z,
            2 => VType::X,
            3 => VType::H,
            _ => VType::B,
        };

        let phase1 = match phase {
            Some(p) => Python::with_gil(|py| pyfraction_to_phase(p.bind(py)))?,
            None => Phase::zero(),
        };

        Ok(self.g.add_vertex_with_data(VData {
            ty: ty1,
            phase: phase1,
            qubit,
            row,
        }))
    }

    #[pyo3(signature = (edge_pairs, edgetype=1))]
    fn add_edges(&mut self, edge_pairs: PyObject, edgetype: u8) -> PyResult<()> {
        Python::with_gil(|py| -> PyResult<()> {
            for eo in edge_pairs.bind(py).try_iter()? {
                let ep = eo?.extract::<(V, V)>()?;
                self.add_edge(ep, edgetype);
            }
            Ok(())
        })
    }

    fn remove_vertex(&mut self, v: V) {
        self.g.remove_vertex(v)
    }

    fn remove_edge(&mut self, e: (usize, usize)) {
        self.g.remove_edge(e.0, e.1)
    }

    fn add_to_phase(&mut self, v: usize, phase: PyObject) -> PyResult<()> {
        let phase1 = Python::with_gil(|py| pyfraction_to_phase(phase.bind(py)))?;
        self.g.add_to_phase(v, phase1);
        Ok(())
    }

    fn num_inputs(&self) -> usize {
        self.g.inputs().len()
    }

    fn num_outputs(&self) -> usize {
        self.g.outputs().len()
    }

    fn set_position(&mut self, vertex: V, q: f64, r: f64) {
        self.g.set_coord(vertex, Coord::new(r, q));
    }

    fn neighbors(&self, vertex: V) -> Vec<V> {
        Vec::from_iter(self.g.neighbors(vertex))
    }

    fn vertex_degree(&self, v: V) -> usize {
        self.g.degree(v)
    }

    fn edge_s(&self, edge: E) -> V {
        edge.0
    }

    fn edge_t(&self, edge: E) -> V {
        edge.1
    }

    fn vertex_set(&self) -> HashSet<V> {
        HashSet::from_iter(self.g.vertices())
    }

    fn edge_set(&self) -> HashSet<E> {
        HashSet::from_iter(self.g.edges().map(|(s, t, _)| (s, t)))
    }

    fn stats(&self) -> String {
        format!("TODO: stats here")
    }

    #[pyo3(signature = (adjoint=false, backend=None))]
    fn copy(&self, adjoint: bool, backend: Option<&str>) -> PyResult<VecGraph> {
        if backend != None && backend != Some("quizx_vec") {
            Err(PyNotImplementedError::new_err(
                "Copy to other backends not implemented on backend: quizx_vec",
            ))
        } else {
            let mut g = self.clone();
            if adjoint {
                g.adjoint();
            }
            Ok(g)
        }
    }

    fn adjoint(&mut self) {
        self.g.adjoint()
    }

    fn map_qubits(&self, qubit_map: HashMap<i32, (f64, f64)>) -> PyResult<()> {
        let _qubit_map = qubit_map;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn compose(&self, other: PyObject) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Compose with other backends not implemented on backend: quizx_vec",
        ))
    }

    fn tensor(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn __iadd__(&self, other: PyObject) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn __add__(&self, other: PyObject) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn __mul__(&self, other: PyObject) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn __matmul__(&self, other: PyObject) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn merge(&self, other: PyObject) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn subgraph_from_vertices(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn apply_state(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn apply_effect(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn to_tensor(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn to_matrix(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn to_dict(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn to_json(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn to_graphml(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn to_tikz(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    // fn from_json(cls, js:Union[str,Dict[str,Any]]) -> VecGraph: ...
    // fn from_tikz(cls, tikz: str, warn_overlap:bool= True, fuse_overlap:bool = True, ignore_nonzx:bool = False) -> VecGraph: ...

    fn is_id(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn pack_circuit_rows(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn qubit_count(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn auto_detect_io(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn normalize(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn translate(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn add_edge_table(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn set_phase_master(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn update_phase_index(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn fuse_phases(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn phase_negate(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn vertex_from_phase_index(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn remove_isolated_vertices(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn vdata_dict(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn set_vdata_dict(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn is_well_formed(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn get_auto_simplify(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn set_auto_simplify(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
    }

    fn is_phase_gadget(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx_vec",
        ))
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

fn pyfraction_to_phase(obj: &Bound<'_, PyAny>) -> PyResult<Phase> {
    let num = obj
        .getattr(intern!(obj.py(), "numerator"))?
        .extract::<i64>()?;
    let denom = obj
        .getattr(intern!(obj.py(), "denominator"))?
        .extract::<i64>()?;
    Ok(Phase::new(Rational64::new(num, denom)))
}
