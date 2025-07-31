use pyo3::{exceptions::PyValueError, prelude::*, types::PyList, IntoPyObjectExt};

use crate::vec_graph::PyVecGraph;
use quizx::graph::V;
use quizx::rankwidth::annealer::RankwidthAnnealer;
use quizx::rankwidth::decomp_tree::{DecompNode, DecompTree};

#[pyclass(name = "DecompTree")]
#[derive(Clone, Debug, Default)]
pub struct PyDecompTree(DecompTree);

#[pymethods]
impl PyDecompTree {
    fn to_list<'py>(&self, py: Python<'py>) -> PyResult<PyObject> {
        if self.0.nodes.is_empty() {
            return Ok(PyList::empty(py).into_any().unbind());
        }

        let n = self.0.nodes[0].nhd()[0];
        let list = PyList::empty(py);
        list.append(self.write_list(py, 0, n)?)?;
        list.append(self.write_list(py, n, 0)?)?;
        Ok(list.into_any().unbind())
    }

    #[staticmethod]
    fn from_list<'py>(py: Python<'py>, list: PyObject) -> PyResult<Self> {
        let list = list.downcast_bound::<PyList>(py)?;
        if list.len() != 2 {
            return Err(PyValueError::new_err(
                "DecompTree must be constructed from nested lists of length 2",
            ));
        }

        let mut tree = PyDecompTree(DecompTree::new());
        let left = tree.read_list(py, list.get_item(0)?.unbind(), 0)?;
        let right = tree.read_list(py, list.get_item(1)?.unbind(), 0)?;
        tree.0.nodes[left].nhd_mut()[0] = right;
        tree.0.nodes[right].nhd_mut()[0] = left;

        Ok(tree)
    }

    /// Create a new empty decomposition tree
    #[new]
    pub fn new() -> Self {
        PyDecompTree(DecompTree::new())
    }

    /// Add a leaf node to the tree
    fn add_leaf(&mut self, parent: usize, vertex: V) -> usize {
        self.0.add_leaf(parent, vertex)
    }

    /// Add an interior node to the tree
    fn add_interior(&mut self, neighbors: [usize; 3]) -> usize {
        self.0.add_interior(neighbors)
    }

    /// Get all edges in the tree
    fn edges(&self, py: Python<'_>) -> PyResult<PyObject> {
        let edges = self.0.edges();
        let py_edges: Vec<(usize, usize)> = edges;
        py_edges.into_py_any(py)
    }

    /// Get the number of edges in the tree
    fn num_edges(&self) -> usize {
        self.0.num_edges()
    }

    /// Find a path between two nodes
    fn path(&self, py: Python<'_>, n1: usize, n2: usize) -> PyResult<PyObject> {
        let path = self.0.path(n1, n2);
        path.into_py_any(py)
    }

    /// Get the vertex partition defined by an edge
    fn partition(&self, py: Python<'_>, edge: (usize, usize)) -> PyResult<PyObject> {
        let (p1, p2) = self.0.partition(edge);
        let result = (p1, p2);
        result.into_py_any(py)
    }

    /// Compute ranks for all edges given a graph
    fn compute_ranks(&mut self, graph: &PyVecGraph) {
        self.0.compute_ranks(&graph.g)
    }

    /// Compute the rankwidth of the tree given a graph
    fn rankwidth(&mut self, graph: &PyVecGraph) -> usize {
        self.0.rankwidth(&graph.g)
    }

    /// Compute the rankwidth score (sum of squares of ranks)
    fn rankwidth_score(&mut self, graph: &PyVecGraph) -> usize {
        self.0.rankwidth_score(&graph.g)
    }

    /// Set the rank of an edge
    fn set_rank(&mut self, edge: (usize, usize), rank: usize) {
        self.0.set_rank(edge, rank)
    }

    /// Clear the rank of an edge
    fn clear_rank(&mut self, edge: (usize, usize)) {
        self.0.clear_rank(edge)
    }

    /// Get the rank of an edge
    fn rank(&mut self, edge: (usize, usize)) -> Option<usize> {
        self.0.rank(edge)
    }

    /// Check if the tree gives a valid partition for a graph
    fn is_valid_for_graph(&self, graph: &PyVecGraph) -> bool {
        self.0.is_valid_for_graph(&graph.g)
    }

    /// Generate a random decomposition tree for a graph
    #[staticmethod]
    #[pyo3(signature = (graph, seed=None))]
    fn random_decomp(graph: &PyVecGraph, seed: Option<u64>) -> Self {
        use rand::{rngs::StdRng, SeedableRng};
        let mut rng = if let Some(s) = seed {
            StdRng::seed_from_u64(s)
        } else {
            StdRng::from_os_rng()
        };
        PyDecompTree(DecompTree::random_decomp(&graph.g, &mut rng))
    }

    /// Perform a random swap of leaves
    #[pyo3(signature = (seed=None))]
    fn swap_random_leaves(&mut self, seed: Option<u64>) {
        use rand::{rngs::StdRng, SeedableRng};
        let mut rng = if let Some(s) = seed {
            StdRng::seed_from_u64(s)
        } else {
            StdRng::from_os_rng()
        };
        self.0.swap_random_leaves(&mut rng)
    }

    /// Perform a random subtree move
    #[pyo3(signature = (seed=None))]
    fn move_random_subtree(&mut self, seed: Option<u64>) {
        use rand::{rngs::StdRng, SeedableRng};
        let mut rng = if let Some(s) = seed {
            StdRng::seed_from_u64(s)
        } else {
            StdRng::from_os_rng()
        };
        self.0.move_random_subtree(&mut rng)
    }

    /// Perform a random local swap
    #[pyo3(signature = (seed=None))]
    fn random_local_swap(&mut self, seed: Option<u64>) {
        use rand::{rngs::StdRng, SeedableRng};
        let mut rng = if let Some(s) = seed {
            StdRng::seed_from_u64(s)
        } else {
            StdRng::from_os_rng()
        };
        self.0.random_local_swap(&mut rng)
    }

    /// Get the number of nodes in the tree
    #[getter]
    fn num_nodes(&self) -> usize {
        self.0.nodes.len()
    }

    /// Clone the decomposition tree
    fn clone(&self) -> Self {
        PyDecompTree(self.0.clone())
    }

    /// String representation for debugging
    fn __repr__(&self) -> String {
        format!(
            "DecompTree(nodes={}, leaves={}, interior={})",
            self.0.nodes.len(),
            self.0.leaves.len(),
            self.0.interior.len()
        )
    }
}

impl PyDecompTree {
    fn write_list<'py>(&self, py: Python<'py>, node: usize, parent: usize) -> PyResult<PyObject> {
        match self.0.nodes[node] {
            DecompNode::Leaf(_, v) => Ok(v.into_py_any(py)?),
            DecompNode::Interior(nhd) => {
                let list = PyList::empty(py);
                for n in nhd {
                    if parent != n {
                        list.append(self.write_list(py, n, node)?)?;
                    }
                }
                Ok(list.into_any().unbind())
            }
        }
    }

    fn read_list<'py>(&mut self, py: Python<'py>, obj: PyObject, parent: usize) -> PyResult<usize> {
        if let Ok(list) = obj.downcast_bound::<PyList>(py) {
            if list.len() != 2 {
                return Err(PyValueError::new_err(
                    "DecompTree must be constructed from nested lists of length 2",
                ));
            }

            let n = self.0.add_interior([parent, 0, 0]);
            let left = self.read_list(py, list.get_item(0)?.unbind(), n)?;
            let right = self.read_list(py, list.get_item(1)?.unbind(), n)?;
            self.0.nodes[n].nhd_mut()[1] = left;
            self.0.nodes[n].nhd_mut()[2] = right;
            Ok(n)
        } else {
            let v = obj.extract::<V>(py)?;
            Ok(self.0.add_leaf(parent, v))
        }
    }
}

#[pyclass(name = "RankwidthAnnealer")]
pub struct PyRankwidthAnnealer(RankwidthAnnealer<rand::rngs::StdRng, quizx::vec_graph::Graph>);

#[pymethods]
impl PyRankwidthAnnealer {
    /// Create a new annealer with a random initial decomposition
    #[new]
    #[pyo3(signature = (graph, seed=None))]
    fn new(graph: &PyVecGraph, seed: Option<u64>) -> Self {
        use rand::{rngs::StdRng, SeedableRng};
        let rng = if let Some(s) = seed {
            StdRng::seed_from_u64(s)
        } else {
            StdRng::from_os_rng()
        };
        PyRankwidthAnnealer(RankwidthAnnealer::new(graph.g.clone(), rng))
    }

    /// Set the initial decomposition tree
    fn set_init_decomp(&mut self, init_decomp: PyDecompTree) {
        self.0.set_init_decomp(init_decomp.0);
    }

    /// Set the initial temperature for simulated annealing
    fn set_init_temp(&mut self, init_temp: f64) {
        self.0.set_init_temp(init_temp);
    }

    /// Set the minimum temperature for simulated annealing
    fn set_min_temp(&mut self, min_temp: f64) {
        self.0.set_min_temp(min_temp);
    }

    /// Set the cooling rate for simulated annealing
    fn set_cooling_rate(&mut self, cooling_rate: f64) {
        self.0.set_cooling_rate(cooling_rate);
    }

    /// Set whether to use adaptive cooling
    fn set_adaptive_cooling(&mut self, adaptive_cooling: bool) {
        self.0.set_adaptive_cooling(adaptive_cooling);
    }

    /// Set the number of iterations to run
    fn set_iterations(&mut self, iterations: usize) {
        self.0.set_iterations(iterations);
    }

    /// Get the initial decomposition tree
    fn init_decomp(&self) -> PyDecompTree {
        PyDecompTree(self.0.init_decomp().clone())
    }

    /// Run the simulated annealing algorithm and return the best decomposition found
    fn run(&mut self) -> PyDecompTree {
        let result = self.0.run();
        PyDecompTree(result)
    }

    /// String representation for debugging
    fn __repr__(&self) -> String {
        "RankwidthAnnealer()".to_string()
    }
}
