"""Cycling protocols (docs/protocol.md).

A protocol string such as ``'cc C=1 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.2; cv V=4.2 Imin=0.05'``
is parsed into a list of Step objects. Both models use it (through nmc_model.driver), and the
Fortran and C++ programs implement the same grammar. `p` is any parameter object with the fields
C_rate, V_min, V_max, steps and cycles.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


KINDS = ("cc", "cv", "rest")
_KEYS = {"cc": {"c", "t", "vmin", "vmax"}, "cv": {"v", "t", "imin"}, "rest": {"t"}}


@dataclass(frozen=True)
class Step:
    kind: str
    C: float = 0.0                  # C-rate for cc (positive = discharge)
    V: float = 0.0                  # held voltage for cv [V]
    t: Optional[float] = None       # maximum duration [s]
    Vmin: float = 0.0
    Vmax: float = 0.0
    Imin: Optional[float] = None    # cv: end when |I| <= Imin (C-rate)


def parse(text: str, p) -> list[Step]:
    """Parse a protocol string; an empty string gives the classic discharge at p.C_rate."""
    text = (text or "").strip()
    if not text:
        return [Step("cc", C=p.C_rate, Vmin=p.V_min, Vmax=p.V_max)]
    steps = []
    for n, part in enumerate(filter(None, (s.strip() for s in text.split(";"))), start=1):
        words = part.split()
        kind = words[0].lower()
        if kind not in KINDS:
            raise ValueError(f"step {n}: unknown step type {words[0]!r} (expected cc, cv or rest)")
        kv = {}
        for w in words[1:]:
            if "=" not in w:
                raise ValueError(f"step {n}: expected key=value, got {w!r}")
            k, v = w.split("=", 1)
            k = k.lower()
            if k not in _KEYS[kind]:
                raise ValueError(f"step {n}: {kind} does not take {k!r}")
            kv[k] = float(v.lower().replace("d", "e"))
        t = kv.get("t")
        if t is not None and t <= 0:
            raise ValueError(f"step {n}: t must be positive")
        if kind == "cc":
            if "c" not in kv:
                raise ValueError(f"step {n}: cc needs C=")
            steps.append(Step("cc", C=kv["c"], t=t, Vmin=kv.get("vmin", p.V_min), Vmax=kv.get("vmax", p.V_max)))
        elif kind == "cv":
            if "v" not in kv:
                raise ValueError(f"step {n}: cv needs V=")
            if t is None and "imin" not in kv:
                raise ValueError(f"step {n}: cv needs t= or Imin= to end")
            steps.append(Step("cv", V=kv["v"], t=t, Imin=kv.get("imin"), Vmin=p.V_min, Vmax=p.V_max))
        else:
            if t is None:
                raise ValueError(f"step {n}: rest needs t=")
            steps.append(Step("rest", t=t, Vmin=p.V_min, Vmax=p.V_max))
    return steps


def expand(p) -> list[Step]:
    return parse(p.steps, p) * max(1, p.cycles)
