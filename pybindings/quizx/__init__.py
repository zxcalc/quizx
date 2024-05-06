from . import _quizx, simplify
from .graph import VecGraph
from .circuit import Circuit
from .decompose import Decomposer
from ._quizx import Scalar

__all__ = ["VecGraph", "Circuit", "simplify", "Decomposer", "Scalar"]


def extract_circuit(g):
    c = Circuit()
    c._c = _quizx.extract_circuit(g._g)
    return c
