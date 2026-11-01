use num::{One, Rational64, Zero};
use pyo3::{prelude::*, types::PyList};
use quizx::{params::Parity, phase::Phase};

pub fn from_fraction_like<'py>(py: Python<'py>, f: PyObject) -> Rational64 {
    // FractionLike can be int, rational, or Poly. If it is Poly, set to zero
    // TODO: return Parity for the Poly case, rather than zero
    if let Ok(p) = f.extract::<Rational64>(py) {
        p
    } else {
        f.extract::<i64>(py).unwrap_or_default().into()
    }
}

pub fn to_fraction_like<'py>(py: Python<'py>, phase: Phase, vars: Parity) -> PyResult<PyObject> {
    let p = if vars.is_empty() {
        phase.to_rational().into_pyobject(py)?
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

        poly.call1((ts,))?
    };

    Ok(p.unbind())
}
