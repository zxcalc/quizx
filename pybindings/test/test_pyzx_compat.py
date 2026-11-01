"""Test the bindings against the PyZX API.

These tests are meant to be run against different versions of the pyzx library
to ensure compatibility. See the 'test-pybindings' job in `ci.yml`.
"""

from fractions import Fraction

import numpy as np
import pytest
import pyzx as zx
from pyzx.graph.scalar import Scalar

import quizx

QASM = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
h q[0];
cx q[0], q[1];
t q[1];
"""


def test_tensor_and_circuit_conversion():
    graph = quizx.qasm(QASM)
    expected = zx.Circuit.from_qasm(QASM).to_matrix()
    np.testing.assert_allclose(graph.to_matrix(), expected, atol=1e-10)
    quizx.full_simp(graph)
    extracted = quizx.extract_circuit(graph)
    assert isinstance(extracted, zx.Circuit)
    np.testing.assert_allclose(extracted.to_matrix(), expected, atol=1e-10)


def test_simplification_and_clone():
    graph = quizx.qasm(QASM)
    expected = graph.to_matrix()
    restored = graph.clone()
    quizx.full_simp(restored)
    np.testing.assert_allclose(restored.to_matrix(), expected, atol=1e-10)


def test_scalar_conversion():
    graph = quizx.VecGraph()
    scalar = Scalar()
    scalar.add_phase(Fraction(1, 4))
    scalar.add_power(2)
    scalar.add_node(Fraction(1, 2))
    graph.scalar = scalar
    np.testing.assert_allclose(graph.scalar.to_number(), scalar.to_number(), atol=1e-10)


def test_zero_scalar_conversion():
    graph = quizx.VecGraph()
    scalar = Scalar()
    scalar.is_zero = True
    graph.scalar = scalar
    assert graph.scalar.to_number() == 0


def test_edge_queries():
    graph = quizx.VecGraph()
    v, w, u = graph.add_vertices(3)
    graph.add_edge((v, w), zx.EdgeType.SIMPLE)
    graph.add_edge((w, u), zx.EdgeType.HADAMARD)
    assert graph.num_edges() == 2
    assert graph.num_edges(v, v) == 0
    assert graph.num_edges(v, w, zx.EdgeType.SIMPLE) == 1
    assert graph.num_edges(v, w, zx.EdgeType.HADAMARD) == 0
    assert graph.num_edges(w) == 2
    assert graph.num_edges(w, et=zx.EdgeType.HADAMARD) == 1
    assert graph.edge(w, v, zx.EdgeType.SIMPLE) == (v, w)


def test_tensor_scalar():
    graph = quizx.qasm(QASM)
    expected = graph.to_matrix()
    scalar = graph.scalar
    scalar.add_power(2)
    graph.scalar = scalar
    np.testing.assert_allclose(graph.to_matrix(), 2 * expected, atol=1e-10)
    np.testing.assert_allclose(
        graph.to_matrix(preserve_scalar=False), expected / 2**0.5, atol=1e-10
    )


def test_indexed_vertices():
    graph = quizx.VecGraph()
    assert graph.add_vertex(zx.VertexType.Z, phase=Fraction(1, 4), index=5) == 5
    assert graph.vertices() == [5]
    assert graph.phase(5) == Fraction(1, 4)
    assert graph.add_vertex(index=0) == 0
    with pytest.raises(ValueError):
        graph.add_vertex(index=5)
    graph.remove_vertex(5)
    assert graph.add_vertex(index=5) == 5


def test_tensor_strategy():
    graph = quizx.qasm(QASM)
    np.testing.assert_allclose(
        graph.to_matrix(strategy="naive"), graph.to_matrix(), atol=1e-10
    )
    with pytest.raises(ValueError, match="tensor strategies"):
        graph.to_tensor(strategy="unknown")


def test_adjoint_tensor_conversion():
    graph = quizx.qasm(QASM)
    original = graph.to_matrix()
    adjoint = graph.adjoint()
    np.testing.assert_allclose(adjoint.to_matrix(), original.conj().T, atol=1e-10)
    np.testing.assert_allclose(graph.to_matrix(), original, atol=1e-10)
