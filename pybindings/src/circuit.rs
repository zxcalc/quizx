use ::quizx::circuit::Circuit;
use ::quizx::gate::GType::*;
use pyo3::prelude::*;
use pyo3::types::IntoPyDict;

pub fn to_pyzx_circuit(py: Python<'_>, c: Circuit) -> PyResult<PyObject> {
    let m = PyModule::import(py, "pyzx.circuit")?;
    let c1 = m.getattr("Circuit")?.call((c.num_qubits(),), None)?;

    for g in c.gates {
        match g.t {
            XPhase => {
                let kwargs = [("phase", g.phase.to_rational())].into_py_dict(py)?;
                c1.call_method("add_gate", ("XPhase", g.qs[0]), Some(&kwargs))?;
            }
            NOT => {
                c1.call_method("add_gate", ("NOT", g.qs[0]), None)?;
            }
            ZPhase => {
                let kwargs = [("phase", g.phase.to_rational())].into_py_dict(py)?;
                c1.call_method("add_gate", ("ZPhase", g.qs[0]), Some(&kwargs))?;
            }
            Z => {
                c1.call_method("add_gate", ("Z", g.qs[0]), None)?;
            }
            S => {
                c1.call_method("add_gate", ("S", g.qs[0]), None)?;
            }
            T => {
                c1.call_method("add_gate", ("T", g.qs[0]), None)?;
            }
            Sdg => {
                let kwargs = [("adjoint", true)].into_py_dict(py)?;
                c1.call_method("add_gate", ("S", g.qs[0]), Some(&kwargs))?;
            }
            Tdg => {
                let kwargs = [("adjoint", true)].into_py_dict(py)?;
                c1.call_method("add_gate", ("T", g.qs[0]), Some(&kwargs))?;
            }
            CNOT => {
                c1.call_method("add_gate", ("CNOT", g.qs[0], g.qs[1]), None)?;
            }
            CZ => {
                c1.call_method("add_gate", ("CZ", g.qs[0], g.qs[1]), None)?;
            }
            ParityPhase => {
                c1.call_method("add_gate", ("ParityPhase", g.qs[0], g.qs[1]), None)?;
            }
            XCX => {
                c1.call_method("add_gate", ("XCX", g.qs[0], g.qs[1]), None)?;
            }
            SWAP => {
                c1.call_method("add_gate", ("SWAP", g.qs[0], g.qs[1]), None)?;
            }
            HAD => {
                c1.call_method("add_gate", ("HAD", g.qs[0]), None)?;
            }
            TOFF => {
                c1.call_method("add_gate", ("TOF", g.qs[0], g.qs[1], g.qs[2]), None)?;
            }
            CCZ => {
                c1.call_method("add_gate", ("CCZ", g.qs[0], g.qs[1], g.qs[2]), None)?;
            }
            InitAncilla => {
                c1.call_method("add_gate", ("InitAncilla", g.qs[0]), None)?;
            }
            PostSelect => {
                c1.call_method("add_gate", ("PostSelect", g.qs[0]), None)?;
            }
            UnknownGate => {}
        }
    }

    Ok(c1.unbind())
}
