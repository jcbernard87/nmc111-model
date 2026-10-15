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
HEADER_CRYSTAL = ("Th_surf_sep", "Th_mean_sep", "Th_surf_col", "Th_mean_col")   # crystal model only
UNITS_CRYSTAL = ("-", "-", "-", "-")


def _fmt_fixed(v: float) -> str:
    return f"{'NaN':>12}" if math.isnan(v) else f"{v:12.5f}"


def _fmt_sci(v: float) -> str:
    return f"{'NaN':>15}" if math.isnan(v) else f"{v:15.5E}"


def format_header(extended: bool = False, crystal: bool = False) -> str:
    def row(cols):
        # Fortran A<w> right-justifies shorter strings and keeps the leftmost w characters of longer ones
        return (f"{cols[0][:5]:>5} " + " ".join(f"{c[:12]:>12}" for c in cols[1:3]) + " "
                + " ".join(f"{c[:15]:>15}" for c in cols[3:]))
    h, u = (HEADER + HEADER_EXTRA, UNITS + UNITS_EXTRA) if extended else (HEADER, UNITS)
    if crystal:
        h, u = h + HEADER_CRYSTAL, u + UNITS_CRYSTAL
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
    columns: str = ""                             # the time-series columns to write ('' = all; nmc_model.output)
    profiles: list = field(default_factory=list)  # [(t [s], x [um], {variable: values})] (nmc_model.output)
    profile_vars: tuple = ()

    def header_units(self, n):
        h, u = HEADER + HEADER_EXTRA + HEADER_CRYSTAL, UNITS + UNITS_EXTRA + UNITS_CRYSTAL
        return h[:n], u[:n]

    def write_columns(self, fh) -> None:
        from .. import output
        n = len(self.rows[0]) if self.rows else 0
        header, units = self.header_units(n)
        chosen = output.choose(header[1:], self.columns, "output column")
        keep = [0, 1] + [header.index(c) for c in chosen if c != "Time"]   # the state and time always
        output.write_table(fh, header, units, self.rows, keep)

    @property
    def array(self) -> np.ndarray:
        """Numeric columns: t [h], V, equivalents, eta [mV], i0 [mA/cm2], c_edge (, I [mA/cm2], step)."""
        return np.array([r[1:] for r in self.rows], dtype=float)

    def write(self, path) -> None:
        with open(path, "w") as fh:
            if self.columns:
                self.write_columns(fh)
                return
            n = len(self.rows[0]) if self.rows else 0
            fh.write(format_header(extended=n > 7, crystal=n > 10))
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
        if p.particle_model == "crystal":
            from .crystal import CrystalModel
            self.model = CrystalModel(p)
        else:
            from .logmodel import LogUniformModel
            self.model = LogUniformModel(p)
        self.p, self.backend = p, backend
        from ..output import ProfileRecorder
        names = ["c", "phi1", "phi2", "theta"] + (["theta_surf"] if p.particle_model == "crystal" else [])
        self.recorder = ProfileRecorder(p, names, self.profile)

    def profile(self, c, I):
        """x [um] and the profile variables at every electrode node (nmc_model.output)."""
        m, p = self.model, self.p
        el = self.electrode(c)
        s = m.mesh.s
        shift = foil_shift(p, el, I)
        nan = np.full(p.nj, np.nan)
        phi1 = el[:, P1] - shift
        phi1[:s] = np.nan
        vals = {"c": p.c_bulk * np.exp(el[:, 0]), "phi1": phi1, "phi2": el[:, P2] - shift}
        if hasattr(c, "xc"):
            th = m.cs(c) / m.cs_max
            mean, surf = nan.copy(), nan.copy()
            mean[m.nodes], surf[m.nodes] = th @ (m.cm.V / m.cm.V.sum()), th[:, -1]
            vals["theta"], vals["theta_surf"] = mean, surf
        else:
            th = 1.0 / (1.0 + np.exp(-el[:, 3]))
            th[:s] = np.nan
            vals["theta"] = th
        return m.mesh.x * 1.0e4, vals

    def initial_state(self):
        return self.model.initial_state()

    def newton_step(self, c, h, I):
        return self.model.newton_step(c, h, I, backend=self.backend)

    # ---- restart files (nmc_model.restart) ----
    def restart_arrays(self, c):
        if hasattr(c, "xc"):
            return [("c", c.x), ("xc", c.xc)]
        if hasattr(c, "ca"):
            return [("c", c.c), ("xc", c.ca.reshape(c.ca.shape[0], -1))]
        return [("c", c)]

    def restart_shapes(self):
        return {name: a.shape for name, a in self.restart_arrays(self.model.initial_state())}

    def from_restart(self, arrays):
        st = self.model.initial_state()
        if hasattr(st, "xc"):
            return type(st)(arrays["c"].copy(), arrays["xc"].copy())
        if hasattr(st, "ca"):
            return type(st)(arrays["c"].copy(), arrays["xc"].reshape(st.ca.shape).copy())
        return arrays["c"].copy()

    @staticmethod
    def flat(c):
        """The state as one vector, for the adaptive driver's error estimate."""
        return np.concatenate((c.x.ravel(), c.xc.ravel())) if hasattr(c, "xc") else c.ravel()

    @staticmethod
    def electrode(c):
        return c.x if hasattr(c, "xc") else c

    def voltage(self, c, I):
        return cell_voltage(self.p, self.electrode(c), I)

    def row(self, t, c, mAhg, I, k):
        self.recorder.row(t, c, I)
        r = _row(self.p, t, self.electrode(c), mAhg, I, k)
        if hasattr(c, "xc"):                     # crystal model: theta at the first and last crystal
            m = self.model
            th = m.cs(c) / m.cs_max
            w = m.cm.V / m.cm.Vtot
            r = r + (float(th[0, -1]), float(th[0] @ w), float(th[-1, -1]), float(th[-1] @ w))
        return r

    @staticmethod
    def finite(c):
        if hasattr(c, "xc"):
            return bool(np.all(np.isfinite(c.x)) and np.all(np.isfinite(c.xc)))
        return bool(np.all(np.isfinite(c)))

    def limit_reason(self, x):
        if hasattr(x, "xc"):
            return self.model.limit_reason(x)
        s = self.model.mesh.s
        th = self.model.cs(x)[s:] / self.model.cs_max
        return driver_limit(float(self.model.conc(x).min()), self.p.c_bulk, float(th.min()), float(th.max()))


def run_with_restarts(stepper, model: str, result, max_steps=None):
    """run_protocol, from p.restart_read if set, writing p.restart_write if set (nmc_model.restart)."""
    from .. import restart
    p = stepper.p
    start = None
    if p.restart_read:
        r = restart.read(p.restart_read, model, stepper.restart_shapes())
        start = (stepper.from_restart(r["arrays"]), r["t"], r["mAhg"], r["I"])
    res = run_protocol(stepper, max_steps=max_steps, result=result, start=start)
    res.columns = p.columns
    stepper.recorder.finish(res, p.profile_file)
    if p.restart_write:
        t, mAhg, I = res.end
        restart.write(p.restart_write, model, t, mAhg, I, stepper.restart_arrays(res.final_state))
    return res


def _run_protocol(p: Params, *, backend: str, pivot: str, max_steps: Optional[int]) -> Result:
    return run_with_restarts(UniformStepper(p, backend=backend, pivot=pivot), p.particle_model, Result(), max_steps)
