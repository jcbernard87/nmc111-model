"""Floating-point helpers that reproduce how the original Fortran programs evaluate expressions.

Faithful mode needs these to match the original output byte for byte (deviation D-4):

- `f32`: a default-kind (single-precision) literal or variable, rounded to float32.
- `powi`: an integer power, as gfortran evaluates it (libgcc ``__powidf2``: repeated squaring,
  with a reciprocal for negative exponents).
- `powr`: a real power, as gfortran evaluates it (libm ``pow``).
"""
from __future__ import annotations

import numpy as np


def f32(x):
    """Round to single precision; returns a float (or a float64 array)."""
    if np.ndim(x) == 0:
        return float(np.float32(x))
    return np.asarray(x, dtype=np.float32).astype(np.float64)


def powi(x, m: int):
    """x**m for an integer m, evaluated like libgcc's __powidf2 (vectorized over x)."""
    x = np.asarray(x, dtype=np.float64)
    n = abs(int(m))
    y = x.copy() if n % 2 else np.ones_like(x)
    n >>= 1
    while n:
        x = x * x
        if n % 2:
            y = y * x
        n >>= 1
    return 1.0 / y if m < 0 else y


def powr(x, e):
    """x**e for a real exponent e, through libm pow (avoids numpy's sqrt/square fast paths)."""
    x = np.asarray(x, dtype=np.float64)
    return np.power(x, np.full(x.shape, float(e)))


def mul32(*factors):
    """A product evaluated entirely in single precision, left to right (e.g. ``4.0/3.0*PI``)."""
    acc = np.float32(factors[0])
    for f in factors[1:]:
        acc = np.float32(acc * np.float32(f))
    return float(acc)
