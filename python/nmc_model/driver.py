"""Corrected-mode time driver shared by both models (docs/protocol.md).

A protocol of constant-current, constant-voltage and rest steps is run with backward Euler,
each step solved to convergence by the model's Newton step. The driver is independent of the
model: it works with any *stepper*, an object with

- ``p``: parameters with the protocol fields (``steps``, ``cycles``, ``C_rate``, ``V_min``,
  ``V_max``, ``write_interval``, ``dt``, ``i_1C``, ``mass_area``)
- ``initial_state()``
- ``newton_step(state, h, I) -> state``: one converged backward-Euler step, or ``SolverFailure``
- ``voltage(state, I) -> float``: the cell voltage against the lithium foil
- ``row(t, state, mAhg, I, k) -> tuple``: one output row
- ``finite(state) -> bool``
- ``limit_reason(state) -> str or None`` (optional): the physical limit the state has reached, if any,
  used as the exit reason when a step cannot be solved (see LIMIT_* below)

Sub-steps are halved when Newton fails, voltage cutoffs are located to within EVENT_DV, and a
constant-voltage step finds its current by bracketing (see cv_step).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Optional

from .protocol import Step, expand


class SolverFailure(RuntimeError):
    pass


EVENT_DV = 1.0e-4      # a cutoff crossing is located to within this voltage [V] ...
MIN_SUBSTEP = 1.0e-10     # [s] Newton failures halve the sub-step down to this, then give up ...
MAX_FAILURES = 200        # ... or give up after this many Newton failures within one time step
# Physical limits reported as the exit reason when a step cannot be solved
LIMIT_DEPLETED = 1.0e-3   # electrolyte below this fraction of c_bulk anywhere: 'electrolyte_depleted'
LIMIT_THETA = 1.0e-3      # particles within this of full (or empty): 'particles_full' / 'particles_empty'


def limit_reason(c_min: float, c_bulk: float, theta_min: float, theta_max: float):
    """The physical limit a state has reached, or None."""
    if c_min < LIMIT_DEPLETED * c_bulk:
        return "electrolyte_depleted"
    if theta_max > 1.0 - LIMIT_THETA:
        return "particles_full"
    if theta_min < LIMIT_THETA:
        return "particles_empty"
    return None
EVENT_MIN_DT = 1.0e-12  # ... or this sub-step length [s]
CV_TOL = 1.0e-9        # [V]


@dataclass
class ProtocolResult:
    rows: list = field(default_factory=list)
    exit_reason: str = ""
    steps: int = 0
    final_state: Any = None


def advance(stepper, state, dt: float, I: float, *, margin=None):
    """Advance by dt with backward Euler, halving the sub-step on Newton failure.

    `margin(state)` is the distance to the nearest voltage cutoff (negative once crossed). A
    sub-step that crosses by more than EVENT_DV is retried with half the length, so the step
    ends within EVENT_DV of the cutoff. Newton failures halve the sub-step, down to MIN_SUBSTEP and
    at most MAX_FAILURES times per step, then raise SolverFailure; after each success the sub-step doubles again (up to dt), so a
    failure does not leave the rest of the step crawling at a tiny sub-step.
    Returns (state, time advanced, stopped).
    """
    t_done, h, failures = 0.0, dt, 0
    while t_done < dt:
        h = min(h, dt - t_done)
        try:
            new = stepper.newton_step(state, h, I)
        except SolverFailure:
            failures += 1
            if h / 2 < MIN_SUBSTEP or failures >= MAX_FAILURES:
                raise
            h = h / 2
            continue
        if margin is not None:
            m = margin(new)
            if m < 0.0:
                if m < -EVENT_DV and h / 2 >= EVENT_MIN_DT:
                    h = h / 2
                    continue
                return new, t_done + h, True
        state, t_done = new, t_done + h
        h = 2.0 * h
    return state, t_done, False


def cv_step(stepper, state, h: float, V_set: float, I_guess: float):
    """One time step at constant voltage: find the current I with V(I) = V_set.

    f(I) = V(I) - V_set decreases with I. A current the cell cannot carry for the whole step
    (Newton fails) is treated as f = +inf when charging and -inf when discharging, which keeps
    f monotone. The root is bracketed (expanding from the previous current, through I = 0) and
    then found by the Illinois variant of regula falsi, with bisection when a bound is infinite.
    """
    p = stepper.p
    states = {}

    def f(I):
        try:
            new = stepper.newton_step(state, h, I)
        except SolverFailure:
            return math.inf if I < 0 else -math.inf
        states[I] = new
        return stepper.voltage(new, I) - V_set

    I = I_guess
    fI = f(I)
    if abs(fI) <= CV_TOL:
        return states[I], I
    # bracket: a < b with f(a) > 0 > f(b)
    grow = max(abs(I), 1.0e-2 * p.i_1C)
    a = b = None
    fa = fb = None
    for _ in range(60):
        if fI > 0:
            a, fa = I, fI
            if b is not None:
                break
            I = 0.0 if I < 0 else I + grow
        else:
            b, fb = I, fI
            if a is not None:
                break
            I = 0.0 if I > 0 else I - grow
        grow *= 2.0
        fI = f(I)
        if abs(fI) <= CV_TOL:
            return states[I], I
    if a is None or b is None:
        raise SolverFailure("constant-voltage current could not be bracketed")
    side = 0
    for _ in range(200):
        if math.isfinite(fa) and math.isfinite(fb):
            I = (a * fb - b * fa) / (fb - fa)
            if not (a < I < b):
                I = 0.5 * (a + b)
        else:
            I = 0.5 * (a + b)
        fI = f(I)
        if abs(fI) <= CV_TOL or (b - a) <= 1.0e-14 * p.i_1C:
            if I in states:
                return states[I], I
            raise SolverFailure("constant-voltage step: no feasible current at the set voltage")
        if fI > 0:
            a, fa = I, fI
            if side == 1 and math.isfinite(fb):
                fb *= 0.5
            side = 1
        else:
            b, fb = I, fI
            if side == -1 and math.isfinite(fa):
                fa *= 0.5
            side = -1
    raise SolverFailure("constant-voltage current iteration did not converge")


def run_protocol(stepper, *, max_steps: Optional[int] = None, result=None) -> ProtocolResult:
    """Run the protocol p.steps (repeated p.cycles times)."""
    p = stepper.p
    steps = expand(p)
    state = stepper.initial_state()
    dt = p.dt
    t, mAhg, n_done = 0.0, 0.0, 0
    res = result if result is not None else ProtocolResult()

    def first_current(st: Step) -> float:
        return st.C * p.i_1C if st.kind == "cc" else 0.0

    def finish(reason, k, I):
        res.rows.append(stepper.row(t, state, mAhg, I, k))
        res.exit_reason, res.steps, res.final_state = reason, n_done, state
        return res

    I = first_current(steps[0])
    res.rows.append(stepper.row(t, state, mAhg, I, 1))
    last_write = t
    reason = ""
    for k, st in enumerate(steps, start=1):
        t_step = 0.0
        if st.kind != "cv":
            I = first_current(st)
        while True:
            if max_steps is not None and n_done >= max_steps:
                res.exit_reason, res.steps, res.final_state = "max_steps", n_done, state
                return res
            h = dt if st.t is None else min(dt, st.t - t_step)
            try:
                if st.kind == "cv":
                    new, I = cv_step(stepper, state, h, st.V, I)
                    h_done = h
                    stopped = st.Imin is not None and abs(I) <= st.Imin * p.i_1C
                    why = "current_limit"
                else:
                    margin = None
                    if st.kind == "cc":
                        def margin(sn, I=I, st=st):
                            v = stepper.voltage(sn, I)
                            return min(v - st.Vmin, st.Vmax - v)
                    new, h_done, stopped = advance(stepper, state, h, I, margin=margin)
                    why = "cutoff_low" if stopped and stepper.voltage(new, I) <= st.Vmin else "cutoff_high"
            except SolverFailure:
                why = getattr(stepper, "limit_reason", lambda s: None)(state)
                return finish(why or "solver_fail", k, I)
            mAhg = mAhg + 1000.0 * (I / p.mass_area) * h_done / 3600.0
            state, t, t_step, n_done = new, t + h_done, t_step + h_done, n_done + 1
            if not stepper.finite(state):
                return finish("nan", k, I)
            if stopped or (st.t is not None and t_step >= st.t * (1.0 - 1.0e-12)):
                res.rows.append(stepper.row(t, state, mAhg, I, k))
                last_write = t
                reason = why if stopped else "duration"
                break
            if t - last_write >= p.write_interval:
                res.rows.append(stepper.row(t, state, mAhg, I, k))
                last_write = t
            if t >= 99.0 * 3600.0:
                return finish("max_time", k, I)
    res.exit_reason = reason if len(steps) == 1 else "end_of_protocol"
    res.steps, res.final_state = n_done, state
    return res
