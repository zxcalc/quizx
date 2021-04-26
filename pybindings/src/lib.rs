use quizx::graph::*;
use quizx::extract::ToCircuit;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use num::Rational;

#[pymodule]
fn libquizx(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dummy, m)?)?;
    m.add_function(wrap_pyfunction!(interior_clifford_simp, m)?)?;
    m.add_function(wrap_pyfunction!(extract_circuit, m)?)?;
    m.add_class::<VecGraph>()?;
    m.add_class::<Circuit>()?;
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
fn extract_circuit(g: &mut VecGraph) -> Circuit {
    Circuit { c: g.g.into_circuit().unwrap() }
}

/// A (mostly) opaque wrapper for quizx circuits
#[pyclass]
struct Circuit { pub c: quizx::circuit::Circuit }

#[pymethods]
impl Circuit {
    #[staticmethod]
    fn from_qasm(qasm: String) -> Circuit {
        Circuit { c: quizx::circuit::Circuit::from_qasm(&qasm).unwrap() }
    }

    #[staticmethod]
    fn load(file: String) -> Circuit {
        Circuit { c: quizx::circuit::Circuit::from_file(&file).unwrap() }
    }

    fn to_qasm(&self) -> String { self.c.to_qasm() }
    fn to_graph(&self) -> VecGraph { VecGraph { g: self.c.to_graph() } }

    fn num_gates(&self) -> usize { self.c.num_gates() }
    fn stats(&self) -> String { self.c.stats().to_string() }
    fn stats_array(&self) -> [usize; 7] { self.c.stats().into_array() }
    fn num_twoq(&self) -> usize { self.c.stats().twoq }
    fn num_oneq(&self) -> usize { self.c.stats().oneq }
}

/// Wrapper for quizx::vec_graph::Graph
#[pyclass]
struct VecGraph { pub g: quizx::vec_graph::Graph }

#[pymethods]
impl VecGraph {
    #[new]
    fn new() -> VecGraph {
        VecGraph { g: quizx::vec_graph::Graph::new() }
    }

    fn vindex(&self) -> usize { self.g.vindex() }
    fn neighbor_at(&self, v: usize, n: usize) -> usize { self.g.neighbor_at(v, n) }
    fn num_vertices(&self) -> usize { self.g.num_vertices() }
    fn num_edges(&self) -> usize { self.g.num_edges() }
    fn add_vertex(&mut self,
                  ty_num: u8,
                  qubit: i32,
                  row: i32,
                  phase: (isize, isize)) -> usize
    {
        let ty = match ty_num {
            1 => VType::Z,
            2 => VType::X,
            3 => VType::H,
            _ => VType::B
        };
        let phase = Rational::new(phase.0, phase.1);
        self.g.add_vertex_with_data(VData {
            ty, phase, qubit, row,
        })
    }

    fn contains_vertex(&self, v: usize) -> bool { self.g.contains_vertex(v) }

    fn add_edge(&mut self, e: (usize, usize), et_num: u8) {
        let et = match et_num { 2 => EType::H, _ => EType::N };
        self.g.add_edge_with_type(e.0, e.1, et)
    }

    fn add_edge_smart(&mut self, e: (usize, usize), et_num: u8) {
        let et = match et_num { 1 => EType::H, _ => EType::N };
        self.g.add_edge_smart(e.0, e.1, et)
    }

    fn remove_vertex(&mut self, v: usize) { self.g.remove_vertex(v) }
    fn remove_edge(&mut self, e: (usize, usize)) { self.g.remove_edge(e.0, e.1) }
    fn degree(&self, v: usize) -> usize { self.g.degree(v) }
    fn connected(&self, s: usize, t: usize) -> bool { self.g.connected(s,t) }

    fn vertex_type(&self, v: usize) -> u8 {
        match self.g.vertex_type(v) {
            VType::Z => 1,
            VType::X => 2,
            VType::H => 3,
            VType::B => 0,
        }
    }

    fn set_vertex_type(&mut self, v: usize, ty_num: u8) {
        let ty = match ty_num {
            1 => VType::Z,
            2 => VType::X,
            3 => VType::H,
            _ => VType::B,
        };
        self.g.set_vertex_type(v, ty);
    }

    fn edge_type(&self, e: (usize, usize)) -> u8 {
        match self.g.edge_type_opt(e.0,e.1) {
            Some(EType::N) => 1,
            Some(EType::H) => 2,
            None => 0,
        }
    }

    fn set_edge_type(&mut self, e: (usize, usize), et_num: u8) {
        self.g.set_edge_type(e.0, e.1,
          if et_num == 2 { EType::H } else { EType::N });
    }

    fn phase(&self, v: usize) -> (isize, isize) {
        let p = self.g.phase(v);
        (*p.numer(), *p.denom())
    }

    fn set_phase(&mut self, v: usize, phase: (isize, isize)) {
        self.g.set_phase(v, Rational::new(phase.0, phase.1));
    }

    fn add_to_phase(&mut self, v: usize, phase: (isize, isize)) {
        self.g.add_to_phase(v, Rational::new(phase.0, phase.1));
    }

    fn qubit(&mut self, v: usize) -> i32 { self.g.qubit(v) }
    fn set_qubit(&mut self, v: usize, q: i32) { self.g.set_qubit(v, q); }
    fn row(&mut self, v: usize) -> i32 { self.g.row(v) }
    fn set_row(&mut self, v: usize, r: i32) { self.g.set_row(v, r); }
}
