use num::complex::Complex64;
use num::One;
use num::Rational64;
use num::Zero;
use pyo3::exceptions::*;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyType};
use quizx::graph::*;
use quizx::phase::*;
use quizx::scalar::Scalar4;
use quizx::vec_graph::Graph;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::scalar::PyScalar4;
use crate::util::*;

type E = (V, V);

/// Wrapper for quizx::vec_graph::Graph
#[pyclass(name = "VecGraph")]
pub struct PyVecGraph {
    pub g: ::quizx::vec_graph::Graph,
    // extra data PyZX expects in BaseGraph
    track_phases: bool,
    phase_index: Py<PyDict>,
    max_phase_index: isize,
    phase_master: Option<Py<PyAny>>,
    phase_mult: Py<PyDict>,
}

impl PyVecGraph {
    pub fn from_graph<'py>(py: Python<'py>, g: Graph) -> PyVecGraph {
        PyVecGraph {
            g,
            track_phases: false,
            phase_index: PyDict::new(py).unbind(),
            max_phase_index: -1,
            phase_master: None,
            phase_mult: PyDict::new(py).unbind(),
        }
    }
}

#[pymethods]
impl PyVecGraph {
    #[new]
    pub fn new<'py>(py: Python<'py>) -> PyVecGraph {
        PyVecGraph::from_graph(py, Graph::new())
    }

    #[classattr]
    fn backend() -> String {
        "quizx-vec".to_string()
    }

    /// Returns the graph scalar.
    ///
    /// Warning: this returns a *copy* of the scalar, so for changes to have
    /// effect, you must set the scalar back afters, e.g. in Python:
    ///     s = g.scalar
    ///     s.add_phase(Fraction(1,2))
    ///     g.set_scalar(s)
    ///
    /// This is different from the behaviour of other backends, but it seems to be
    /// unavoidable due to Rust's ownership limitations.
    ///
    #[getter]
    fn scalar<'py>(&mut self, py: Python<'py>) -> PyResult<PyObject> {
        to_pyzx_scalar(py, self.g.scalar())
    }

    /// Sets the graph scalar.
    #[setter]
    fn set_scalar<'py>(&mut self, py: Python<'py>, scalar: PyObject) -> PyResult<()> {
        *self.g.scalar_mut() = from_pyzx_scalar(py, scalar)?;
        Ok(())
    }

    /// Get a wrapped quizx native scalar
    ///
    /// This can be used to avoid losing precision converting to/from PyZX scalars
    #[getter]
    fn scalar4(&self) -> PyScalar4 {
        self.g.scalar().clone().into()
    }

    /// Set a wrapped quizx native scalar
    #[setter]
    fn set_scalar4(&mut self, s: PyScalar4) {
        *self.g.scalar_mut() = s.into();
    }

    fn clone<'py>(&self, py: Python<'py>) -> PyResult<PyVecGraph> {
        Ok(PyVecGraph {
            g: self.g.clone(),
            track_phases: false,
            phase_index: PyDict::new(py).unbind(),
            max_phase_index: -1,
            phase_master: None,
            phase_mult: PyDict::new(py).unbind(),
        })
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
        Vec::from_iter((0..amount).map(|_| self.g.add_vertex(VType::B)))
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
            3 => EType::Wio,
            _ => EType::N,
        };
        self.g.add_edge_smart(edge_pair.0, edge_pair.1, et)
    }

    fn remove_vertices(&mut self, vertices: &Bound<'_, PyAny>) -> PyResult<()> {
        for vo in vertices.try_iter()? {
            let v = vo?.extract::<V>()?;
            self.g.remove_vertex(v);
        }
        Ok(())
    }

    fn remove_edges(&mut self, edges: &Bound<'_, PyAny>) -> PyResult<()> {
        for eo in edges.try_iter()? {
            let e = eo?.extract::<E>()?;
            self.g.remove_edge(e.0, e.1);
        }
        Ok(())
    }

    fn num_vertices(&self) -> usize {
        self.g.num_vertices()
    }

    #[pyo3(signature = (s=None, t=None, et=None))]
    fn num_edges(&self, s: Option<V>, t: Option<V>, et: Option<u8>) -> usize {
        // For now, we ignore the filtering parameters and just return total count
        let _s = s;
        let _t = t;
        let _et = et;
        if s.is_some() || t.is_some() || et.is_some() {
            println!(
                "warning: num_edges filtering parameters are not implemented for quizx-vec backend"
            );
        }
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

    fn phase(&self, py: Python<'_>, v: usize) -> PyResult<PyObject> {
        let (phase, vars) = self.g.phase_and_vars(v);
        to_fraction_like(py, phase, vars)
    }

    fn set_phase(&mut self, v: usize, phase: Rational64) {
        self.g.set_phase(v, Phase::new(phase));
    }

    fn qubit(&self, v: usize) -> f64 {
        self.g.qubit(v)
    }

    fn set_qubit(&mut self, v: usize, q: f64) {
        self.g.set_qubit(v, q);
    }

    fn row(&self, v: usize) -> f64 {
        self.g.row(v)
    }

    fn set_row(&mut self, v: usize, r: f64) {
        self.g.set_row(v, r);
    }

    fn clear_vdata(&self, vertex: V) -> PyResult<()> {
        let _vertex = vertex;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn vdata_keys(&self, vertex: V) -> PyResult<Vec<String>> {
        let _vertex = vertex;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    #[pyo3(signature = (vertex, key, default=None))]
    fn vdata(
        &self,
        vertex: V,
        key: String,
        default: Option<PyObject>,
    ) -> PyResult<Option<PyObject>> {
        let _vertex = vertex;
        let _key = key;
        println!("warning: vdata is not implemented for quizx-vec backend");
        Ok(default)
    }

    fn set_vdata(&self, vertex: V, key: String, val: PyObject) -> PyResult<()> {
        let _vertex = vertex;
        let _key = key;
        let _val = val;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
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
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn is_hybrid(&self) -> bool {
        false
    }

    fn multigraph(&self) -> bool {
        false
    }

    fn phases(&self, py: Python<'_>) -> PyResult<HashMap<V, PyObject>> {
        let mut m = HashMap::default();
        for v in self.g.vertices() {
            m.insert(v, self.phase(py, v)?);
        }
        Ok(m)
    }

    fn types(&self) -> HashMap<V, u8> {
        let mut m = HashMap::default();
        for v in self.g.vertices() {
            m.insert(v, self.vertex_type(v));
        }
        m
    }

    fn qubits(&self) -> HashMap<V, f64> {
        let mut m = HashMap::default();
        for v in self.g.vertices() {
            m.insert(v, self.qubit(v));
        }
        m
    }

    fn rows(&self) -> HashMap<V, f64> {
        let mut m = HashMap::default();
        for v in self.g.vertices() {
            m.insert(v, self.row(v));
        }
        m
    }

    fn depth(&self) -> f64 {
        self.g.depth()
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
        phase: Option<Rational64>,
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
            Some(p) => Phase::new(p),
            None => Phase::zero(),
        };

        Ok(self.g.add_vertex_with_data(VData {
            ty: ty1,
            phase: phase1,
            qubit,
            row,
            ..Default::default()
        }))
    }

    #[pyo3(signature = (edge_pairs, edgetype=1))]
    fn add_edges(&mut self, edge_pairs: &Bound<'_, PyAny>, edgetype: u8) -> PyResult<()> {
        for eo in edge_pairs.try_iter()? {
            let ep = eo?.extract::<(V, V)>()?;
            self.add_edge(ep, edgetype);
        }
        Ok(())
    }

    fn remove_vertex(&mut self, v: V) {
        self.g.remove_vertex(v)
    }

    fn remove_edge(&mut self, e: (usize, usize)) {
        self.g.remove_edge(e.0, e.1)
    }

    fn add_to_phase(&mut self, v: usize, phase: Rational64) -> PyResult<()> {
        let phase1 = Phase::new(phase);
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

    fn vindex(&self) -> usize {
        println!("warning: vindex is not fully implemented for quizx-vec backend");
        self.g.num_vertices()
    }

    fn clear_edata(&self, _edge: E) -> PyResult<()> {
        println!("warning: clear_edata is not implemented for quizx-vec backend");
        Ok(())
    }

    fn edata_keys(&self, _edge: E) -> PyResult<Vec<String>> {
        println!("warning: edata_keys is not implemented for quizx-vec backend");
        Ok(vec![])
    }

    #[pyo3(signature = (edge, key, default=None))]
    fn edata(
        &self,
        py: Python<'_>,
        edge: E,
        key: String,
        default: Option<PyObject>,
    ) -> PyResult<PyObject> {
        let _edge = edge;
        let _key = key;
        println!("warning: edata is not implemented for quizx-vec backend");
        Ok(default.unwrap_or_else(|| py.None()))
    }

    fn set_edata(&self, _edge: E, _key: String, _val: PyObject) -> PyResult<()> {
        println!("warning: set_edata is not implemented for quizx-vec backend");
        Ok(())
    }

    fn stats(&self) -> String {
        let mut degrees: HashMap<usize, usize> = HashMap::new();
        for v in self.g.vertices() {
            let d = self.g.degree(v);
            degrees.insert(d, 1 + *degrees.get(&d).unwrap_or(&0));
        }

        let mut s = format!(
            "VecGraph ({} vertices, {} edges)\nDegrees:",
            self.g.num_vertices(),
            self.g.num_edges()
        );
        for e in degrees {
            s += &format!("  {}: {}\n", e.0, e.1);
        }
        s
    }

    #[pyo3(signature = (adjoint=false, backend=None))]
    fn copy<'py>(
        &self,
        py: Python<'py>,
        adjoint: bool,
        backend: Option<&str>,
    ) -> PyResult<PyVecGraph> {
        if backend.is_some() && backend != Some("quizx-vec") {
            Err(PyNotImplementedError::new_err(
                "Copy to other backends not implemented on backend: quizx-vec",
            ))
        } else {
            Ok(PyVecGraph::from_graph(py, self.g.copy(adjoint)))
        }
    }

    fn adjoint(&mut self) {
        self.g.adjoint()
    }

    fn map_qubits(&self, qubit_map: HashMap<i32, (f64, f64)>) -> PyResult<()> {
        let _qubit_map = qubit_map;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn compose(&mut self, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let other1 = other
            .downcast::<PyVecGraph>()
            .map_err(|_| {
                PyNotImplementedError::new_err(
                    "Operations with mixed backends not implemented on backend: quizx-vec",
                )
            })?
            .borrow();

        // TODO: GraphLike::plug seems to not be working right
        // TODO: fix coords
        self.g.plug(&other1.g);
        Ok(())
    }

    fn tensor(&mut self, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let other1 = other
            .downcast::<PyVecGraph>()
            .map_err(|_| {
                PyNotImplementedError::new_err(
                    "Operations with mixed backends not implemented on backend: quizx-vec",
                )
            })?
            .borrow();

        let mp = self.g.append_graph(&other1.g);

        for i in other1.g.inputs() {
            self.g.inputs_mut().push(mp[i]);
        }

        for o in other1.g.outputs() {
            self.g.outputs_mut().push(mp[o]);
        }

        // TODO: fix coords
        Ok(())
    }

    fn __iadd__(&mut self, other: &Bound<'_, PyAny>) -> PyResult<()> {
        // TODO: should return a reference to Self. I don't know how to do this w PyO3
        self.compose(other)
    }

    fn __add__<'py>(&self, py: Python<'py>, other: &Bound<'_, PyAny>) -> PyResult<PyVecGraph> {
        let mut g = self.clone(py)?;
        g.compose(other)?;
        Ok(g)
    }

    fn __mul__<'py>(&self, py: Python<'py>, other: &Bound<'_, PyAny>) -> PyResult<PyVecGraph> {
        let mut other1 = other
            .downcast::<PyVecGraph>()
            .map_err(|_| {
                PyNotImplementedError::new_err(
                    "Operations with mixed backends not implemented on backend: quizx-vec",
                )
            })?
            .borrow()
            .clone(py)?;
        other1.g.plug(&self.g);
        Ok(other1)
    }

    fn __matmul__<'py>(&self, py: Python<'py>, other: &Bound<'_, PyAny>) -> PyResult<PyVecGraph> {
        let mut other1 = other
            .downcast::<PyVecGraph>()
            .map_err(|_| {
                PyNotImplementedError::new_err(
                    "Operations with mixed backends not implemented on backend: quizx-vec",
                )
            })?
            .borrow()
            .clone(py)?;
        other1.g.plug(&self.g);
        Ok(other1)
    }

    fn merge(&self, other: &Bound<'_, PyAny>) -> PyResult<(Vec<V>, Vec<E>)> {
        println!("warning: merge is not implemented for quizx-vec backend");
        let _other = other;
        Ok((vec![], vec![]))
    }

    fn subgraph_from_vertices<'py>(&self, py: Python<'py>, verts: Vec<V>) -> PyVecGraph {
        PyVecGraph::from_graph(py, self.g.subgraph_from_vertices(verts))
    }

    fn apply_state(&mut self, state: String) {
        let e = Vec::from_iter(state.chars().map(|c| match c {
            '0' => BasisElem::Z0,
            '1' => BasisElem::Z1,
            '+' => BasisElem::X0,
            '-' => BasisElem::X1,
            _ => BasisElem::SKIP,
        }));
        self.g.plug_inputs(&e);
    }

    fn apply_effect(&mut self, state: String) {
        let e = Vec::from_iter(state.chars().map(|c| match c {
            '0' => BasisElem::Z0,
            '1' => BasisElem::Z1,
            '+' => BasisElem::X0,
            '-' => BasisElem::X1,
            _ => BasisElem::SKIP,
        }));
        self.g.plug_outputs(&e);
    }

    #[pyo3(signature = (preserve_scalar=true))]
    fn to_tensor(&self, py: Python<'_>, preserve_scalar: bool) -> PyResult<PyObject> {
        let m = PyModule::import(py, "pyzx.tensor")?;
        let f = m.getattr("tensorfy")?;
        Ok(f.call((self.clone(py)?, preserve_scalar), None)?.unbind())
    }

    #[pyo3(signature = (preserve_scalar=true))]
    fn to_matrix(&self, py: Python<'_>, preserve_scalar: bool) -> PyResult<PyObject> {
        let m = PyModule::import(py, "pyzx.tensor")?;
        let tensorfy = m.getattr("tensorfy")?;
        let tensor_to_matrix = m.getattr("tensor_to_matrix")?;

        let tensor = tensorfy.call((self.clone(py)?, preserve_scalar), None)?;
        Ok(tensor_to_matrix
            .call((tensor, self.num_inputs(), self.num_outputs()), None)?
            .unbind())
    }

    #[pyo3(signature = (include_scalar=true))]
    fn to_dict(&self, include_scalar: bool) -> PyResult<PyObject> {
        let _include_scalar = include_scalar;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    #[pyo3(signature = (include_scalar=true))]
    fn to_json(&self, include_scalar: bool) -> PyResult<String> {
        let _include_scalar = include_scalar;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn to_graphml(&self) -> PyResult<String> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    #[pyo3(signature = (draw_scalar=false))]
    fn to_tikz(&self, draw_scalar: bool) -> PyResult<String> {
        let _draw_scalar = draw_scalar;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    // fn from_json(cls, js:Union[str,Dict[str,Any]]) -> VecGraph: ...
    // fn from_tikz(cls, tikz: str, warn_overlap:bool= True, fuse_overlap:bool = True, ignore_nonzx:bool = False) -> VecGraph: ...

    fn is_id(&self) -> bool {
        println!("warning: is_id is not fully implemented for quizx-vec backend");
        let inputs = self.inputs();
        let outputs = self.outputs();

        if inputs.len() != outputs.len()
            || self.num_vertices() != 2 * inputs.len()
            || self.num_edges(None, None, None) != inputs.len()
        {
            return false;
        }

        for i in 0..inputs.len() {
            if !self.connected(inputs[i], outputs[i]) {
                return false;
            }
        }
        true
    }

    fn pack_circuit_rows(&self) -> PyResult<()> {
        println!("warning: pack_circuit_rows is not implemented for quizx-vec backend");
        Ok(())
    }

    fn qubit_count(&self) -> usize {
        println!("warning: qubit_count is not fully implemented for quizx-vec backend");
        self.num_inputs()
    }

    fn auto_detect_io(&self) -> PyResult<()> {
        println!("warning: auto_detect_io is not implemented for quizx-vec backend");
        Ok(())
    }

    fn normalize(&self) -> PyResult<()> {
        println!("warning: normalize is not implemented for quizx-vec backend");
        Ok(())
    }

    fn translate(&mut self, x: f64, y: f64) -> PyResult<()> {
        println!("warning: translate is not fully implemented for quizx-vec backend");
        let vertices: Vec<V> = self.g.vertices().collect();
        for v in vertices {
            self.g.set_row(v, self.g.row(v) + x);
            self.g.set_qubit(v, self.g.qubit(v) + y);
        }
        Ok(())
    }

    fn add_edge_table(&self, _etab: PyObject) -> PyResult<()> {
        println!("warning: add_edge_table is not implemented for quizx-vec backend");
        Ok(())
    }

    fn update_phase_index(&self, old: V, new: V) -> PyResult<()> {
        println!("warning: update_phase_index is not implemented for quizx-vec backend");
        let _old = old;
        let _new = new;
        Ok(())
    }

    fn fuse_phases(&self, p1: V, p2: V) -> PyResult<()> {
        println!("warning: fuse_phases is not implemented for quizx-vec backend");
        let _p1 = p1;
        let _p2 = p2;
        Ok(())
    }

    fn phase_negate(&self, v: V) -> PyResult<()> {
        println!("warning: phase_negate is not implemented for quizx-vec backend");
        let _v = v;
        Ok(())
    }

    fn vertex_from_phase_index(&self, i: isize) -> PyResult<V> {
        let _i = i;
        Err(PyNotImplementedError::new_err(
            "vertex_from_phase_index not implemented on backend: quizx-vec",
        ))
    }

    fn remove_isolated_vertices(&self) -> PyResult<()> {
        println!("warning: remove_isolated_vertices is not implemented for quizx-vec backend");
        Ok(())
    }

    fn vdata_dict<'py>(&self, py: Python<'py>, _v: PyObject) -> PyResult<Py<PyDict>> {
        println!("warning: set_vdata_dict is not implemented for quizx-vec backend");
        Ok(PyDict::new(py).unbind())
    }

    fn set_vdata_dict(&self, _v: PyObject, _d: PyObject) -> PyResult<()> {
        println!("warning: set_vdata_dict is not implemented for quizx-vec backend");
        Ok(())
    }

    fn edata_dict<'py>(&self, py: Python<'py>, _e: PyObject) -> PyResult<Py<PyDict>> {
        println!("warning: set_edata_dict is not implemented for quizx-vec backend");
        Ok(PyDict::new(py).unbind())
    }

    fn set_edata_dict(&self, _e: PyObject, _d: PyObject) -> PyResult<()> {
        println!("warning: set_edata_dict is not implemented for quizx-vec backend");
        Ok(())
    }

    fn is_well_formed(&self) -> bool {
        println!("warning: is_well_formed is not fully implemented for quizx-vec backend");
        // Basic check: all boundary vertices should have degree 1
        for v in self.g.vertices() {
            let vtype = self.vertex_type(v);
            if vtype == 0 && self.g.degree(v) != 1 {
                // VertexType::BOUNDARY = 0
                return false;
            }
        }
        true
    }

    fn get_auto_simplify(&self) -> bool {
        true
    }

    fn set_auto_simplify(&self) {}

    fn is_phase_gadget(&self, py: Python<'_>, v: V) -> bool {
        println!("warning: is_phase_gadget is not fully implemented for quizx-vec backend");
        let vtype = self.vertex_type(v);
        // Basic check: must be Z or X spider with phase 0 and degree >= 2
        if (vtype == 1 || vtype == 2) && self.g.degree(v) >= 2 {
            // Z or X spider
            // Check if phase is 0 (would need to extract the phase properly)
            let phase_result = self.phase(py, v);
            if let Ok(_phase_obj) = phase_result {
                // Would need to check if phase is actually 0, but for now just return false
                return false;
            }
        }
        false
    }

    // Properties for BaseGraph compatibility fields
    #[getter]
    fn track_phases(&self) -> bool {
        self.track_phases
    }

    #[setter]
    fn set_track_phases(&mut self, value: bool) {
        self.track_phases = value;
    }

    #[getter]
    fn phase_index<'py>(&self, py: Python<'py>) -> PyResult<Py<PyDict>> {
        Ok(self.phase_index.clone_ref(py))
    }

    #[setter]
    fn set_phase_index(&mut self, value: Py<PyDict>) {
        self.phase_index = value;
    }

    #[getter]
    fn max_phase_index(&self) -> isize {
        self.max_phase_index
    }

    #[setter]
    fn set_max_phase_index(&mut self, value: isize) {
        self.max_phase_index = value;
    }

    #[getter]
    fn phase_master<'py>(&self, py: Python<'py>) -> Option<Py<PyAny>> {
        self.phase_master.as_ref().map(|pm| pm.clone_ref(py))
    }

    #[setter]
    // These methods are no-ops for the quizx-vec backend, as it doesn't support
    // Poly or custom merge_vdata

    fn set_phase_master(&mut self, value: Option<Py<PyAny>>) {
        println!("warning: set_phase_master is not fully implemented for quizx-vec backend");
        self.phase_master = value;
    }

    #[getter]
    fn phase_mult<'py>(&self, py: Python<'py>) -> PyResult<Py<PyDict>> {
        Ok(self.phase_mult.clone_ref(py))
    }

    #[setter]
    fn set_phase_mult(&mut self, value: Py<PyDict>) {
        self.phase_mult = value;
    }

    // These methods are no-ops for the quizx-vec backend, as it doesn't support
    // Poly or custom merge_vdata

    #[getter]
    fn variable_types<'py>(&self, py: Python<'py>) -> PyResult<Py<PyDict>> {
        Ok(PyDict::new(py).unbind())
    }

    #[setter]
    fn set_variable_types(&mut self, _d: Py<PyDict>) {
        // No-op for quizx-vec, as it does not track variable types
    }

    #[getter]
    fn merge_vdata(&self) -> Option<PyObject> {
        None
    }

    #[setter]
    fn set_merge_vdata(&mut self, _d: Option<PyObject>) {
        // No-op for quizx-vec, as it does not use merge_vdata
    }

    // Class methods
    #[classmethod]
    fn from_json<'py>(
        _cls: &Bound<'py, PyType>,
        _py: Python<'py>,
        js: PyObject,
    ) -> PyResult<PyVecGraph> {
        let _js = js;
        Err(PyNotImplementedError::new_err(
            "from_json not implemented on backend: quizx-vec",
        ))
    }

    #[classmethod]
    #[pyo3(signature = (tikz, warn_overlap=true, fuse_overlap=true, ignore_nonzx=false))]
    fn from_tikz<'py>(
        _cls: &Bound<'py, PyType>,
        _py: Python<'py>,
        tikz: String,
        warn_overlap: bool,
        fuse_overlap: bool,
        ignore_nonzx: bool,
    ) -> PyResult<PyVecGraph> {
        let _tikz = tikz;
        let _warn_overlap = warn_overlap;
        let _fuse_overlap = fuse_overlap;
        let _ignore_nonzx = ignore_nonzx;
        Err(PyNotImplementedError::new_err(
            "from_tikz not implemented on backend: quizx-vec",
        ))
    }
}

/// Returns the scalar as a tuple of four integers.
pub fn from_pyzx_scalar<'py>(py: Python<'py>, pyzx_scalar: PyObject) -> PyResult<Scalar4> {
    let is_zero = pyzx_scalar.getattr(py, "is_zero")?.extract::<bool>(py)?;
    if is_zero {
        return Ok(Scalar4::zero());
    }

    let mut s = Scalar4::one();

    let phase = from_fraction_like(py, pyzx_scalar.getattr(py, "phase")?);
    s.mul_phase(phase);

    let power2 = pyzx_scalar
        .getattr(py, "power2")?
        .extract::<i32>(py)
        .unwrap_or_default();
    s.mul_sqrt2_pow(power2);

    pyzx_scalar
        .getattr(py, "phasenodes")?
        .extract::<Vec<PyObject>>(py)?
        .into_iter()
        .for_each(|f| {
            let s1 = Scalar4::one_plus_phase(from_fraction_like(py, f));
            s += s1;
        });

    let floatfactor: Scalar4 = pyzx_scalar
        .getattr(py, "floatfactor")?
        .extract::<Complex64>(py)?
        .into();

    if !floatfactor.is_one() {
        s *= floatfactor;
    }

    Ok(s)
}

pub fn to_pyzx_scalar<'py>(py: Python<'py>, s: &Scalar4) -> PyResult<PyObject> {
    let m = PyModule::import(py, "pyzx.graph.scalar")?;
    let scalar_class = m.getattr("Scalar")?;
    let scalar = scalar_class.call1(())?;

    if let Some((phase, pow)) = s.exact_phase_and_sqrt2_pow() {
        scalar.setattr("phase", Rational64::from(phase).into_pyobject(py)?)?;
        scalar.setattr("power2", pow.into_pyobject(py)?)?;
    } else {
        scalar.setattr("floatfactor", s.complex_value().into_pyobject(py)?)?;
    }

    Ok(scalar.unbind())
}
