// use std::string;

// use crate::scalar::Scalar;
use crate::vec_graph::VecGraph;
use crate::Scalar;
// use num::integer;
use pyo3::prelude::*;
use quizx::vec_graph::Graph;

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
    #[pyo3(signature = (g, *, use_cats=None, save=None, simp=None, random_t=None))]
    fn new(
        g: &VecGraph,
        use_cats: Option<bool>,
        save: Option<bool>,
        simp: Option<SimpFunc>,
        random_t: Option<bool>,
    ) -> Decomposer {
        let mut d = ::quizx::decompose::Decomposer::new(&g.g);

        if let Some(use_cats) = use_cats {
            d.use_cats(use_cats);
        }

        if let Some(save) = save {
            d.save(save);
        }

        if let Some(simp) = simp {
            d.with_simp(simp.into());
        }

        if let Some(random_t) = random_t {
            d.random_t(random_t);
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

    fn random_t(&mut self, b: bool) {
        self.d.random_t(b);
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

    fn with_simp(&mut self, simp: SimpFunc) {
        self.d.with_simp(simp.into());
    }

    fn with_full_simp(&mut self) {
        self.d.with_full_simp();
    }

    fn with_clifford_simp(&mut self) {
        self.d.with_clifford_simp();
    }

    fn random_t(&mut self, b: bool) {
        self.d.random_t(b);
    }

    fn max_terms(&self) -> f64 {
        self.d.max_terms()
    }

    fn decompose(&mut self) {
        self.d.decompose();

    fn decomp_top(&mut self) {
        self.d.decomp_top();
    }
    fn decomp_all(&mut self) {
        self.d.decomp_all();
    }


    fn decompose_until_depth(&mut self, depth: i64) {
        self.d.decomp_until_depth(depth);
    }

    fn use_cats(&mut self, b: bool) {
        self.d.use_cats(b);
    }
    fn get_nterms(&self) -> usize {
        self.d.nterms
    }

    fn split(&mut self) -> Vec<Decomposer> {
        self.d
            .clone()
            .split()
            .into_iter()
            .map(|d| Decomposer { d })
            .collect()
    }

    #[staticmethod]
    fn merge(ds: &Bound<pyo3::types::PyList>) -> PyResult<Decomposer> {
        let mut rust_decomposers = Vec::new();

        for item in ds.iter() {
            let decomposer: PyRef<Decomposer> = item.extract()?;
            rust_decomposers.push(decomposer.d.clone());
        }

        Ok(Decomposer {
            d: ::quizx::decompose::Decomposer::merge(rust_decomposers),
        })
    }

    fn pop_graph(&mut self) -> VecGraph {
        VecGraph {
            g: self.d.pop_graph(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Decomposer(nterms={}, scalar={}, stack_size={}, done_size={})",
            self.d.nterms,
            self.d.scalar,
            self.d.stack.len(),
            self.d.done.len()
        )
    }
}
