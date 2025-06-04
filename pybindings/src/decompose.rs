// use crate::scalar::Scalar;
use crate::vec_graph::VecGraph;

use pyo3::prelude::*;

#[pyclass]
pub struct Decomposer {
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
    fn decomp_all(&mut self) {
        self.d.decompose();
    }
    fn decomp_until_depth(&mut self, depth: i64) {
        self.d.decomp_until_depth(depth);
    }
    fn decomp_parallel(&mut self) {
        self.d.decompose_parallel();
    }
    fn get_nterms(&self) -> usize {
        self.d.nterms
    }
}
