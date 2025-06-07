// use std::string;

// use crate::scalar::Scalar;
use crate::vec_graph::VecGraph;
use crate::Scalar;
// use num::integer;
use pyo3::prelude::*;
use quizx::decompose::Driver;

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

    pub fn scalar(&self) -> Scalar {
        self.d.scalar().into()
    }

    fn with_simp(&mut self, b: bool, clifford_only: bool) {
        if b && clifford_only {
            self.d.with_simp(::quizx::decompose::SimpFunc::CliffordSimp);
        } else if b {
            self.d.with_simp(::quizx::decompose::SimpFunc::FullSimp);
        } else {
            self.d.with_simp(::quizx::decompose::SimpFunc::NoSimp);
        }
    }

    fn with_split_graph_components(&mut self, b: bool) {
        self.d.with_split_graphs_components(b);
    }

    fn save(&mut self, b: bool) {
        self.d.with_save(b);
    }

    fn with_driver(&mut self, driver_type: &str, random_t: bool) {
        match driver_type {
            "BssTOnly" => {
                self.d.with_driver(Driver::BssTOnly(random_t));
            }
            "BssWithCats" => {
                self.d.with_driver(Driver::BssWithCats(random_t));
            }
            _ => {
                println!("Driver Not Supported!");
            }
        };
    }

    fn max_terms(&self) -> f64 {
        self.d.max_terms()
    }

    fn decompose(&mut self) {
        self.d.decompose();
    }

    fn decompose_parallel(&mut self) {
        self.d.decompose_parallel();
    }

    fn decompose_until_depth(&mut self, depth: i64) {
        self.d.decomp_until_depth(depth);
    }

    fn get_nterms(&self) -> usize {
        self.d.nterms
    }
}
