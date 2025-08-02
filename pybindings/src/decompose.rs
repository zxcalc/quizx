use quizx::decompose::{Decomposer, SimpFunc};

use crate::vec_graph::PyVecGraph;
use crate::PyScalar;
use pyo3::prelude::*;

#[pyclass(name = "SimpFunc")]
#[derive(Clone, Debug)]
pub enum PySimpFunc {
    FullSimp,
    CliffordSimp,
    NoSimp,
}

impl From<PySimpFunc> for SimpFunc {
    fn from(s: PySimpFunc) -> Self {
        match s {
            PySimpFunc::FullSimp => SimpFunc::FullSimp,
            PySimpFunc::CliffordSimp => SimpFunc::CliffordSimp,
            PySimpFunc::NoSimp => SimpFunc::NoSimp,
        }
    }
}

#[pyclass(name = "Decomposer")]
pub struct PyDecomposer {
    d: Decomposer<::quizx::vec_graph::Graph>,
}

#[pymethods]
impl PyDecomposer {
    #[staticmethod]
    fn empty() -> PyDecomposer {
        PyDecomposer {
            d: Decomposer::empty(),
        }
    }

    #[new]
    #[pyo3(signature = (g, *, save=None, simp=None))]
    fn new(g: &PyVecGraph, save: Option<bool>, simp: Option<PySimpFunc>) -> PyDecomposer {
        let mut d = Decomposer::new(&g.g);

        if let Some(save) = save {
            d.with_save(save);
        }

        if let Some(simp) = simp {
            d.with_simp(simp.into());
        }

        PyDecomposer { d }
    }

    fn done<'py>(&self, py: Python<'py>) -> PyResult<Vec<PyVecGraph>> {
        let mut gs = vec![];
        for g in &self.d.done {
            gs.push(PyVecGraph::from_graph(py, g.clone()));
        }
        Ok(gs)
    }

    fn with_simp(&mut self, simp: PySimpFunc) {
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

    fn with_save(&mut self, b: bool) {
        self.d.with_save(b);
    }

    fn max_terms(&self) -> f64 {
        self.d.max_terms()
    }

    #[pyo3(signature = (driver_type = "BssWithCats", random_t = false, sherlock_tries = Vec::new()))]
    fn decompose(&mut self, driver_type: &str, random_t: bool, sherlock_tries: Vec<usize>) {
        match driver_type {
            "BssTOnly" => {
                self.d
                    .decompose(&quizx::decompose::BssTOnlyDriver { random_t });
            }
            "BssWithCats" => {
                self.d
                    .decompose(&quizx::decompose::BssWithCatsDriver { random_t });
            }
            "DynamicT" => {
                self.d.decompose(&quizx::decompose::DynamicTDriver);
            }
            "Sherlock" => {
                self.d.decompose(&quizx::decompose::SherlockDriver {
                    tries: sherlock_tries,
                });
            }
            "SpiderCutting" => {
                self.d.decompose(&quizx::decompose::SpiderCuttingDriver);
            }
            _ => {
                println!("Driver Not Supported!");
            }
        };
    }

    #[pyo3(signature = (/, *, allow_threads=true, driver_type, random_t = false, sherlock_tries = Vec::new()))]
    fn decompose_parallel(
        &mut self,
        allow_threads: bool,
        driver_type: &str,
        random_t: bool,
        sherlock_tries: Vec<usize>,
    ) {
        if allow_threads {
            // Release the GIL for potentially long-running parallel computation
            pyo3::Python::with_gil(|py| {
                py.allow_threads(|| {
                    match driver_type {
                        "BssTOnly" => {
                            self.d
                                .decompose_parallel(&quizx::decompose::BssTOnlyDriver { random_t });
                        }
                        "BssWithCats" => {
                            self.d
                                .decompose_parallel(&quizx::decompose::BssWithCatsDriver {
                                    random_t,
                                });
                        }
                        "DynamicT" => {
                            self.d.decompose_parallel(&quizx::decompose::DynamicTDriver);
                        }
                        "Sherlock" => {
                            self.d
                                .decompose_parallel(&quizx::decompose::SherlockDriver {
                                    tries: sherlock_tries,
                                });
                        }
                        "SpiderCutting" => {
                            self.d
                                .decompose_parallel(&quizx::decompose::SpiderCuttingDriver);
                        }
                        _ => {
                            println!("Driver Not Supported!");
                        }
                    };
                });
            });
        }
    }

    #[pyo3(signature = (depth, driver_type, random_t = false, sherlock_tries = Vec::new()))]
    fn decompose_until_depth(
        &mut self,
        depth: i64,
        driver_type: &str,
        random_t: bool,
        sherlock_tries: Vec<usize>,
    ) {
        match driver_type {
            "BssTOnly" => {
                self.d
                    .decompose_until_depth(depth, &quizx::decompose::BssTOnlyDriver { random_t });
            }
            "BssWithCats" => {
                self.d.decompose_until_depth(
                    depth,
                    &quizx::decompose::BssWithCatsDriver { random_t },
                );
            }
            "DynamicT" => {
                self.d
                    .decompose_until_depth(depth, &quizx::decompose::DynamicTDriver);
            }
            "Sherlock" => {
                self.d.decompose_until_depth(
                    depth,
                    &quizx::decompose::SherlockDriver {
                        tries: sherlock_tries,
                    },
                );
            }
            "SpiderCutting" => {
                self.d
                    .decompose_until_depth(depth, &quizx::decompose::SpiderCuttingDriver);
            }
            _ => {
                println!("Driver Not Supported!");
            }
        };
    }

    fn get_nterms(&self) -> usize {
        self.d.nterms
    }

    fn get_scalar(&self) -> PyScalar {
        self.d.scalar().into()
    }

    fn __repr__(&self) -> String {
        format!(
            "Decomposer(nterms={}, done_size={})",
            self.d.nterms,
            self.d.done.len()
        )
    }
}
