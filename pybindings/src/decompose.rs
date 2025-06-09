use crate::scalar::Scalar;
use crate::vec_graph::VecGraph;

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

    #[getter]
    fn get_scalar(&self) -> Scalar {
        self.d.scalar.into()
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

    fn decomp_top(&mut self) {
        self.d.decomp_top();
    }

    fn decomp_all(&mut self) {
        self.d.decomp_all();
    }

    fn decomp_until_depth(&mut self, depth: usize) {
        self.d.decomp_until_depth(depth);
    }

    #[pyo3(signature = (depth, /, *, allow_threads=true))]
    fn decomp_parallel(&mut self, depth: usize, allow_threads: bool) {
        if allow_threads {
            // Release the GIL for potentially long-running parallel computation
            pyo3::Python::with_gil(|py| {
                py.allow_threads(|| {
                    self.d = self.d.clone().decomp_parallel(depth);
                });
            });
        } else {
            self.d = self.d.clone().decomp_parallel(depth);
        }
    }

    fn decomp_ts(
        &mut self,
        depth: usize,
        g: &VecGraph,
        ts: &Bound<pyo3::types::PyList>,
    ) -> PyResult<()> {
        let res: Vec<usize> = ts.extract()?;
        let rs_ts: &[usize] = &res;
        self.d.decomp_ts(depth, g.g.clone(), rs_ts);
        Ok(())
    }

    #[staticmethod]
    fn first_ts(g: &VecGraph) -> Vec<usize> {
        ::quizx::decompose::Decomposer::<Graph>::first_ts(&g.g)
    }

    #[staticmethod]
    #[pyo3(signature = (g, seed=None))]
    fn random_ts(g: &VecGraph, seed: Option<u64>) -> Vec<usize> {
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        ::quizx::decompose::Decomposer::<Graph>::random_ts(&g.g, &mut rng)
    }

    #[staticmethod]
    fn cat_ts(g: &VecGraph) -> Vec<usize> {
        ::quizx::decompose::Decomposer::<Graph>::cat_ts(&g.g)
    }

    fn use_cats(&mut self, b: bool) {
        self.d.use_cats(b);
    }

    #[getter]
    fn nterms(&self) -> usize {
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
