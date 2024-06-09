from typing import List, Optional

from . import _quizx
from .graph import VecGraph
from .scalar import from_pyzx_scalar, to_pyzx_scalar
from pyzx.graph.scalar import Scalar


class Decomposer(object):
    _d: _quizx.Decomposer

    def __init__(self, graph: Optional[VecGraph] = None):
        if graph is None:
            self._d = _quizx.Decomposer.empty()
        else:
            self._d = _quizx.Decomposer(graph.get_raw_graph())

    def graphs(self) -> List[VecGraph]:
        return [VecGraph(g) for g in self._d.graphs()]

    def apply_optimizations(self, b: bool):
        self._d.apply_optimizations(b)

    def max_terms(self):
        return self._d.max_terms()

    def decomp_top(self):
        self._d.decomp_top()

    def decomp_all(self):
        self._d.decomp_all()

    def decomp_until_depth(self, depth: int):
        self._d.decomp_until_depth(depth)

    def use_cats(self, b: bool):
        self._d.use_cats(b)

    def get_nterms(self):
        return self._d.get_nterms()

    @property
    def scalar(self) -> Scalar:
        return to_pyzx_scalar(self._d.scalar)

    @scalar.setter
    def scalar(self, s: Scalar):
        self._d.scalar = from_pyzx_scalar(s)

    def is_ground(self, vertex):
        return False
