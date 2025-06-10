// use std::string;
// use crate::scalar::Scalar;
use crate::vec_graph::VecGraph;
use crate::Scalar;
// use num::integer;
use pyo3::prelude::*;
use quizx::decompose::Driver;

#[pyclass]
#[derive(Clone, Debug)]
pub enum SimpFunc {
    FullSimp,
    CliffordSimp,
    NoSimp,
}

impl From<SimpFunc> for ::quizx::decompose::SimpFunc {
    fn from(s: SimpFunc) -> Self {
        match s {
            SimpFunc::FullSimp => ::quizx::decompose::SimpFunc::FullSimp,
            SimpFunc::CliffordSimp => ::quizx::decompose::SimpFunc::CliffordSimp,
            SimpFunc::NoSimp => ::quizx::decompose::SimpFunc::NoSimp,
        }
    }
}

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
    #[pyo3(signature = (g, *, use_driver=None, save=None, simp=None))]
    fn new(
        g: &VecGraph,
        use_driver: Option<&str>,
        save: Option<bool>,
        simp: Option<SimpFunc>,
    ) -> Decomposer {
        let mut d = ::quizx::decompose::Decomposer::new(&g.g);

        if let Some(driver_str) = use_driver {
            match driver_str {
                "BssTOnly" => {
                    d.with_driver(Driver::BssTOnly(false));
                }
                "BssWithCats" => {
                    d.with_driver(Driver::BssWithCats(false));
                }
                _ => {}
            }
        }

        if let Some(save) = save {
            d.with_save(save);
        }

        if let Some(simp) = simp {
            d.with_simp(simp.into());
        }

        Decomposer { d }
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

    fn with_simp(&mut self, simp: SimpFunc) {
        self.d.with_simp(simp.into());
    }

    fn with_full_simp(&mut self) {
        self.d.with_full_simp();
    }

    fn with_clifford_simp(&mut self) {
        self.d.with_clifford_simp();
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

    fn decompose_until_depth(&mut self, depth: i64) {
        self.d.decomp_until_depth(depth);
    }

    fn get_nterms(&self) -> usize {
        self.d.nterms
    }

    fn __repr__(&self) -> String {
        format!(
            "Decomposer(nterms={}, done_size={})",
            self.d.nterms,
            self.d.done.len()
        )
    }
}
