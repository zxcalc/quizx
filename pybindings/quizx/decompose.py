from typing import List, Optional

import libquizx

from .graph import VecGraph


class Decomposer(object):
    def __init__(self, graph: Optional[VecGraph] = None):
        if graph is None:
            self._d = libquizx.Decomposer.empty()
        self._d = libquizx.Decomposer(graph.get_raw_graph())

    def graphs(self) -> List[VecGraph]:
        return [VecGraph(g) for g in self._d.graphs()]

    def apply_optimizations(self, b: bool): self._d.apply_optimizations(b)
    def max_terms(self): return self._d.max_terms()
    def decomp_top(self): self._d.decomp_top()
    def decomp_all(self): self._d.decomp_all()
    def decomp_until_depth(self, depth: int): self._d.decomp_until_depth(depth)
