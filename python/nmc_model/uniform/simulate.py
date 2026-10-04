"""Time integration and output (docs/model.md sections 7 and 9, docs/protocol.md).

Faithful mode reproduces the original constant-current discharge: one linearized solve
per step, output rows written before each step, and the original exit tests.
Corrected mode runs a protocol of cc / cv / rest steps, each time step solved with Newton.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

import bandsolver

from . import kinetics
from .model import Assembler, C, CS, P1, P2
from .params import Params, f32
from ..driver import limit_reason as driver_limit, run_protocol

HEADER = ("State", "Time", "Voltage", "Equivalence", "Anode_Eta", "anode_exchange_c", "Edge_c0")
UNITS = ("CDR", "hours", "Volts", "electron_equivs", "mV", "mA/cm2", "mol/cm3")
HEADER_EXTRA = ("Current", "Step", "Li_Nernst")   # corrected mode only
UNITS_EXTRA = ("mA/cm2", "#", "mV")


def _fmt_fixed(v: float) -> str:
    return f"{'NaN':>12}" if math.isnan(v) else f"{v:12.5f}"


def _fmt_sci(v: float) -> str:
    return f"{'NaN':>15}" if math.isnan(v) else f"{v:15.5E}"


def format_header(extended: bool = False) -> str:
    def row(cols):
        # Fortran A<w> right-justifies shorter strings and keeps the leftmost w characters of longer ones
        return (f"{cols[0][:5]:>5} " + " ".join(f"{c[:12]:>12}" for c in cols[1:3]) + " "
                + " ".join(f"{c[:15]:>15}" for c in cols[3:]))
    h, u = (HEADER + HEADER_EXTRA, UNITS + UNITS_EXTRA) if extended else (HEADER, UNITS)
    return row(h) + "\n" + row(u) + "\n"


def format_row(state: str, values) -> str:
    t, v, *rest = values
    cols = [f"{x:15d}" if isinstance(x, int) else _fmt_sci(x) for x in rest]
    return f"{state:>5} {_fmt_fixed(t)} {_fmt_fixed(v)} " + " ".join(cols) + "\n"


@dataclass
class Result:
    rows: list = field(default_factory=list)      # (state, t_h, V, equiv, eta_mV, i0_mA, c_edge[, I_mA, step])
    exit_reason: str = ""
    steps: int = 0
    final_state: Optional[np.ndarray] = None

    @property
    def array(self) -> np.ndarray:
        """Numeric columns: t [h], V, equivalents, eta [mV], i0 [mA/cm2], c_edge (, I [mA/cm2], step)."""
        return np.array([r[1:] for r in self.rows], dtype=float)

    def write(self, path) -> None:
        with open(path, "w") as fh:
            fh.write(format_header(extended=bool(self.rows) and len(self.rows[0]) > 7))
            for r in self.rows:
                fh.write(format_row(r[0], r[1:]))


def initial_state(p: Params) -> np.ndarray:
    c = np.empty((p.nj, 4))
    c[:, C] = p.c_bulk
    c[:, P1] = p.phi1_init
    c[:, P2] = p.phi2_init
    c[:, CS] = p.cs_init
    return c


def _output_row(p: Params, state: str, t: float, c: np.ndarray, mAhg: float):
    lit36 = f32(3.6) if p.mode == "faithful" else 3.6
    equiv = mAhg * p.M * lit36 / p.F
    i0_li = float(kinetics.li_exchange_current(p, c[0, C]))
    alpha = 0.5
    with np.errstate(invalid="ignore", divide="ignore"):
        if state == "C":
            eta = 0.5 * math.log(p.i_app / i0_li) / (alpha * p.F / (p.R * p.T))
        elif state == "D":
            eta = -(0.5 * math.log(p.i_app / i0_li)) / (alpha * p.F / (p.R * p.T))
        else:
            eta = 0.0
    return (state, t / float(3600), c[-1, P1] + eta, equiv, eta * 1.0e3, i0_li * 1.0e3, c[0, C])


def run(p: Params, *, backend: str = "fortran", pivot: Optional[str] = None, max_steps: Optional[int] = None) -> Result:
    """Run the model: the original discharge in faithful mode, the protocol in corrected mode."""
    if p.mode == "faithful":
        return _run_faithful(p, backend=backend, pivot=pivot or "legacy", max_steps=max_steps)
    return _run_protocol(p, backend=backend, pivot=pivot or "partial", max_steps=max_steps)


def _run_faithful(p: Params, *, backend: str, pivot: str, max_steps: Optional[int]) -> Result:
    """The original program's constant-current discharge, step for step."""
    asm = Assembler(p)
    c = initial_state(p)
    dt = p.dt
    t = 0.0
    mAhg = 0.0
    state = "D" if p.C_rate >= 0 else "C"
    last_write = 0                          # an integer in the original (deviation D-4)
    write_every = p.t_max / p.n_steps / 200
    dc = np.zeros_like(c)
    res = Result()
    n = p.n_steps if max_steps is None else max_steps

    for it in range(1, n + 1):
        if it == 1:
            res.rows.append(_output_row(p, state, t, c, mAhg))
        elif (t - last_write) / 3600 >= write_every:
            res.rows.append(_output_row(p, state, t, c, mAhg))
            last_write = int(t - dt)
        elif it >= p.n_steps:
            res.rows.append(_output_row(p, state, t, c, mAhg))
        elif c[-1, P1] >= 99.0 and state == "C":
            res.rows.append(_output_row(p, state, t, c, mAhg))
            res.exit_reason = "end_of_charge"
            break
        elif math.isnan(dc[0, C]):
            res.rows.append(_output_row(p, state, t, c, mAhg))
            res.exit_reason = "nan"
            break
        elif t >= 99.0 * 3600.0:
            res.rows.append(_output_row(p, state, t, c, mAhg))
            res.exit_reason = "max_time"
            break

        # Coulomb counting happens before the solve in the original
        if state == "D":
            mAhg = mAhg + 1000.0 * p.i_specific * dt / 3600.0
        elif state == "C":
            mAhg = mAhg - 1000.0 * p.i_specific * dt / 3600.0

        A, B, D, G, _ = asm.assemble(c, dt)
        try:
            dc = bandsolver.solve(A, B, D, G, pivot=pivot, backend=backend)
        except (bandsolver.NonFiniteError, bandsolver.SingularBlockError):
            dc = np.full_like(c, np.nan)
        c = c + dc
        res.steps = it

        if state == "R":
            dt = dt * 1.0001
        else:
            dt = p.t_max / float(p.n_steps)
        t = t + dt
    else:
        res.exit_reason = res.exit_reason or "max_steps"

    res.final_state = c
    return res


# ------------------------------------------------------------------ corrected mode: protocols

def li_eta(p: Params, c_edge: float, I: float) -> float:
    """Signed overpotential of the lithium counter electrode (symmetric Butler-Volmer, fixes D-12).

    Reported as a negative number on discharge (I > 0).
    """
    return -float(kinetics.li_foil(p, float(c_edge), I)[1])


def _c_foil(p: Params, x: np.ndarray) -> float:
    """Electrolyte concentration at the foil from a corrected-mode state (column 0 is u = ln(c/c_bulk))."""
    return float(p.c_bulk * math.exp(x[0, C]))


def foil_shift(p: Params, c: np.ndarray, I: float) -> float:
    """Potential of the lithium foil metal on the solver's scale [V] (corrected mode).

    The equations depend only on potential differences, so the solver fixes the gauge with
    phi2 = 0 at the foil face. The foil metal then sits at U_Li + eta_Li; subtracting it from
    phi1 and phi2 gives potentials referenced to the foil (0 V), exactly. (Imposing the foil
    potential as the boundary condition instead makes the node-0 block singular at rest.)
    """
    u_li, eta_li = kinetics.li_foil(p, _c_foil(p, c), I)
    return float(u_li + eta_li)


def foil_referenced(p: Params, c: np.ndarray, I: float) -> np.ndarray:
    """State with phi1 and phi2 referenced to the lithium foil (0 V)."""
    out = c.copy()
    s = foil_shift(p, c, I)
    out[:, P1] -= s
    out[:, P2] -= s
    return out


def cell_voltage(p: Params, c: np.ndarray, I: float) -> float:
    """Cell voltage against the lithium foil: phi1 at the current collector, foil-referenced (D-16)."""
    return float(c[-1, P1]) - foil_shift(p, c, I)


def _row(p: Params, t: float, c: np.ndarray, mAhg: float, I: float, step: int):
    state = "D" if I > 0 else ("C" if I < 0 else "R")
    c0 = _c_foil(p, c)
    u_li, eta_li = kinetics.li_foil(p, c0, I)
    i0 = float(kinetics.li_exchange_current(p, c0))
    return (state, t / 3600.0, cell_voltage(p, c, I), mAhg * p.M * 3.6 / p.F, -float(eta_li) * 1.0e3,
            i0 * 1.0e3, c0, I * 1.0e3, step, float(u_li) * 1.0e3)


class UniformStepper:
    """The uniform-particle model (corrected mode, log variables) as a stepper for nmc_model.driver."""

    def __init__(self, p: Params, *, backend: str, pivot: str):
        from .logmodel import LogUniformModel
        self.p, self.model, self.backend = p, LogUniformModel(p), backend

    def initial_state(self):
        return self.model.initial_state()

    def newton_step(self, c, h, I):
        return self.model.newton_step(c, h, I, backend=self.backend)

    def voltage(self, c, I):
        return cell_voltage(self.p, c, I)

    def row(self, t, c, mAhg, I, k):
        return _row(self.p, t, c, mAhg, I, k)

    @staticmethod
    def finite(c):
        return bool(np.all(np.isfinite(c)))

    def limit_reason(self, x):
        s = self.model.mesh.s
        th = self.model.cs(x)[s:] / self.model.cs_max
        return driver_limit(float(self.model.conc(x).min()), self.p.c_bulk, float(th.min()), float(th.max()))


def _run_protocol(p: Params, *, backend: str, pivot: str, max_steps: Optional[int]) -> Result:
    return run_protocol(UniformStepper(p, backend=backend, pivot=pivot), max_steps=max_steps, result=Result())
