import libquizx
from .graph import VecGraph
from .circuit import Circuit
from . import simplify
from .decompose import Decomposer

def extract_circuit(g):
    c = Circuit()
    c._c = libquizx.extract_circuit(g._g)
    return c

