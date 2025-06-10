use ::quizx::graph::*;
use ::quizx::phase::*;
use num::Rational64;
use num::Zero;
use pyo3::exceptions::*;
use pyo3::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::scalar::Scalar;
use crate::util::phase_and_vars_to_py;

type E = (V, V);

/// Wrapper for quizx::vec_graph::Graph
#[pyclass]
pub struct VecGraph {
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
    fn get_scalar(&mut self) -> Scalar {
        (*self.g.scalar()).into()
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

    fn phase(&self, py: Python<'_>, v: usize) -> PyResult<PyObject> {
        let (phase, vars) = self.g.phase_and_vars(v);
        phase_and_vars_to_py(py, phase, vars)
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

    fn vdata(&self, vertex: V, key: String, default: PyObject) -> PyResult<()> {
        let _vertex = vertex;
        let _key = key;
        let _default = default;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
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
    fn copy(&self, adjoint: bool, backend: Option<&str>) -> PyResult<VecGraph> {
        if backend.is_some() && backend != Some("quizx-vec") {
            Err(PyNotImplementedError::new_err(
                "Copy to other backends not implemented on backend: quizx-vec",
            ))
        } else {
            Ok(VecGraph {
                g: self.g.copy(adjoint),
            })
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
            .downcast::<VecGraph>()
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
            .downcast::<VecGraph>()
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

    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<VecGraph> {
        let mut g = self.clone();
        g.compose(other)?;
        Ok(g)
    }

    fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<VecGraph> {
        let mut other1 = other
            .downcast::<VecGraph>()
            .map_err(|_| {
                PyNotImplementedError::new_err(
                    "Operations with mixed backends not implemented on backend: quizx-vec",
                )
            })?
            .borrow()
            .clone();
        other1.g.plug(&self.g);
        Ok(other1)
    }

    fn __matmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<VecGraph> {
        let mut other1 = other
            .downcast::<VecGraph>()
            .map_err(|_| {
                PyNotImplementedError::new_err(
                    "Operations with mixed backends not implemented on backend: quizx-vec",
                )
            })?
            .borrow()
            .clone();
        other1.g.plug(&self.g);
        Ok(other1)
    }

    fn merge(&self, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let _other = other;
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn subgraph_from_vertices(&self, verts: Vec<V>) -> VecGraph {
        VecGraph {
            g: self.g.subgraph_from_vertices(verts),
        }
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
        Ok(f.call((self.clone(), preserve_scalar), None)?.unbind())
    }

    #[pyo3(signature = (preserve_scalar=true))]
    fn to_matrix(&self, py: Python<'_>, preserve_scalar: bool) -> PyResult<PyObject> {
        let m = PyModule::import(py, "pyzx.tensor")?;
        let tensorfy = m.getattr("tensorfy")?;
        let tensor_to_matrix = m.getattr("tensor_to_matrix")?;

        let tensor = tensorfy.call((self.clone(), preserve_scalar), None)?;
        Ok(tensor_to_matrix
            .call((tensor, self.num_inputs(), self.num_outputs()), None)?
            .unbind())
    }

    fn to_dict(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn to_json(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn to_graphml(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn to_tikz(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    // fn from_json(cls, js:Union[str,Dict[str,Any]]) -> VecGraph: ...
    // fn from_tikz(cls, tikz: str, warn_overlap:bool= True, fuse_overlap:bool = True, ignore_nonzx:bool = False) -> VecGraph: ...

    fn is_id(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn pack_circuit_rows(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn qubit_count(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn auto_detect_io(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn normalize(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn translate(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn add_edge_table(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn set_phase_master(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn update_phase_index(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn fuse_phases(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn phase_negate(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn vertex_from_phase_index(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn remove_isolated_vertices(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn vdata_dict(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn set_vdata_dict(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn is_well_formed(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }

    fn get_auto_simplify(&self) -> bool {
        true
    }

    fn set_auto_simplify(&self) {}

    fn is_phase_gadget(&self) -> PyResult<()> {
        Err(PyNotImplementedError::new_err(
            "Not implemented on backend: quizx-vec",
        ))
    }
}
