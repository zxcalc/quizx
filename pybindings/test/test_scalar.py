from pyzx.graph.scalar import Scalar as PyzxScalar
from quizx._quizx import Scalar as QuizxScalar
from quizx.scalar import from_pyzx_scalar, to_pyzx_scalar


def test_scalars():
    s = PyzxScalar()
    s.add_phase(0.25)

    q = from_pyzx_scalar(s)
    assert isinstance(q, QuizxScalar)

    s2 = to_pyzx_scalar(q)
    assert isinstance(s2, PyzxScalar)

    assert s2.phase == s.phase
