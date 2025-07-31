use pyo3::{exceptions::PyValueError, prelude::*, types::PyList, IntoPyObjectExt};

use quizx::graph::V;
use quizx::rankwidth::decomp_tree::{DecompNode, DecompTree};

#[pyclass(name = "DecompTree")]
#[derive(Clone, Debug)]
pub struct PyDecompTree(DecompTree);

#[pymethods]
impl PyDecompTree {
    fn to_list<'py>(&self, py: Python<'py>) -> PyResult<PyObject> {
        if self.0.nodes.is_empty() {
            return Ok(PyList::empty(py).into_any().unbind());
        }

        let n = self.0.nodes[0].nhd()[0];
        let list = PyList::empty(py);
        list.append(self.to_list_helper(py, 0, n)?)?;
        list.append(self.to_list_helper(py, n, 0)?)?;
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
        let left = tree.from_list_helper(py, list.get_item(0)?.unbind(), 0)?;
        let right = tree.from_list_helper(py, list.get_item(1)?.unbind(), 0)?;
        tree.0.nodes[left].nhd_mut()[0] = right;
        tree.0.nodes[right].nhd_mut()[0] = left;

        Ok(tree)
    }
}

impl PyDecompTree {
    fn to_list_helper<'py>(
        &self,
        py: Python<'py>,
        node: usize,
        parent: usize,
    ) -> PyResult<PyObject> {
        match self.0.nodes[node] {
            DecompNode::Leaf(_, v) => Ok(v.into_py_any(py)?),
            DecompNode::Interior(nhd) => {
                let list = PyList::empty(py);
                for n in nhd {
                    if parent != n {
                        list.append(self.to_list_helper(py, n, node)?)?;
                    }
                }
                Ok(list.into_any().unbind())
            }
        }
    }

    fn from_list_helper<'py>(
        &mut self,
        py: Python<'py>,
        obj: PyObject,
        parent: usize,
    ) -> PyResult<usize> {
        if let Ok(list) = obj.downcast_bound::<PyList>(py) {
            if list.len() != 2 {
                return Err(PyValueError::new_err(
                    "DecompTree must be constructed from nested lists of length 2",
                ));
            }

            let n = self.0.add_interior([parent, 0, 0]);
            let left = self.from_list_helper(py, list.get_item(0)?.unbind(), n)?;
            let right = self.from_list_helper(py, list.get_item(1)?.unbind(), n)?;
            self.0.nodes[n].nhd_mut()[1] = left;
            self.0.nodes[n].nhd_mut()[2] = right;
            Ok(n)
        } else {
            let v = obj.extract::<V>(py)?;
            Ok(self.0.add_leaf(parent, v))
        }
    }
}
