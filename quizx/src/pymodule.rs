// use crate::graph::GraphLike;
// use crate::vec_graph::*;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pymodule]
fn libquizx(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dummy, m)?)?;
    m.add_class::<TestClass>()?;
    Ok(())
}

#[pyfunction]
fn dummy(a: i64) -> String {
    format!("FOO! {}", a)
}

#[pyclass]
struct TestClass {
    pub x: i64
}

#[pymethods]
impl TestClass {
    #[new]
    fn new() -> TestClass { TestClass { x: 12 } }
}
