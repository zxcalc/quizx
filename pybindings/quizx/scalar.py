"""Utilities for interfacing `quizx::Scalar` with `pyzx.Scalar`."""

from ._quizx import Scalar as QuizxScalar
from pyzx.graph.scalar import Scalar as PyzxScalar


def from_pyzx_scalar(s: PyzxScalar) -> QuizxScalar:
    """Convert a `pyzx.Scalar` to a `quizx::Scalar`."""
    return QuizxScalar.from_json(s.to_json())


def to_pyzx_scalar(s: QuizxScalar) -> PyzxScalar:
    """Convert a `pyzx.Scalar` to a `quizx::Scalar`."""
    return PyzxScalar().from_json(s.to_json())
