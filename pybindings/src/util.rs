use num::{One, Rational64, Zero};
use pyo3::{prelude::*, types::PyList};
use quizx::{params::Parity, phase::Phase};

pub fn phase_and_vars_to_py(py: Python<'_>, phase: Phase, vars: Parity) -> PyResult<PyObject> {
    let p;
    if vars.is_empty() {
        p = phase.to_rational().into_pyobject(py)?;
    } else {
        let m = PyModule::import(py, "pyzx.symbolic")?;
        let poly = m.getattr("Poly")?;
        let term = m.getattr("Term")?;
        let var = m.getattr("Var")?;
        let mut ts = vec![];

        if !phase.is_zero() {
            // add a constant term for the phase if it != 0
            ts.push((phase.to_rational(), term.call1((PyList::empty(py),))?));
        }

        for v in vars.iter() {
            // add a linear term to the Poly for each var
            let py_v = (var.call1((format!("b{v}"), true))?, 1);
            ts.push((Rational64::one(), term.call1((PyList::new(py, [py_v])?,))?));
        }

        p = poly.call1((ts,))?;
    }

    Ok(p.unbind())
}
