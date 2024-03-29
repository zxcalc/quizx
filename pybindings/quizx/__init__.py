from . import _quizx
from .graph import VecGraph
from .circuit import Circuit
from . import simplify
from .decompose import Decomposer

__all__ = ["VecGraph", "Circuit", "simplify", "Decomposer"]


def extract_circuit(g):
    c = Circuit()
    c._c = _quizx.extract_circuit(g._g)
    return c
