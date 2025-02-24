from quizx.graph import VecGraph
from pyzx import EdgeType


def test_add_edge():
    g = VecGraph()
    v1, v2, v3 = g.add_vertices(3)
    e1 = (v1, v2)
    g.add_edge(e1)
    assert g.edge_type(e1) == EdgeType.SIMPLE
    assert g.num_edges() == 1
    e2 = (v1, v3)
    g.add_edge(e2, EdgeType.HADAMARD)
    assert g.edge_type(e2) == EdgeType.HADAMARD
    assert g.num_edges() == 2
