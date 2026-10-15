"""Configurable output (#12): the time-series columns to write, and profiles against position.

Profiles go to one long table: two header rows (names, units), then one row per node per profile time with the
columns Time [h], x [um] and the chosen variables, each 15 characters wide (ES15.7). They are written on output
rows, so they never change the time steps: with every `profile_every`-th row, at the first row at or after each
of `profile_times`, and with the last row. A value that does not exist at a node (the solid potential and theta in
the separator, theta where a node has no particle) is NaN.
"""
from __future__ import annotations

import math

PROFILE_UNITS = {"c": "mol/cm3", "phi1": "V", "phi2": "V", "theta": "-", "theta_surf": "-"}


def names(text: str) -> list:
    return [w for w in text.replace(",", " ").split() if w]


def times(text: str) -> list:
    return sorted(float(w) for w in names(text))


def choose(available, wanted: str, what: str) -> list:
    """The wanted names (case-insensitive) in the available order; all of them if `wanted` is empty."""
    if not wanted.strip():
        return list(available)
    low = {a.lower(): a for a in available}
    out = []
    for w in names(wanted):
        if w.lower() not in low:
            raise ValueError(f"unknown {what} '{w}'")
        out.append(low[w.lower()])
    return [a for a in available if a in out]


def _sci(v: float, w: int = 15, d: int = 5) -> str:
    return f"{'NaN':>{w}}" if math.isnan(v) else f"{v:{w}.{d}E}"


def write_table(fh, header, units, rows, keep) -> None:
    """The time series with the columns `keep` (indices into header; 0, the state, always first). Widths and
    formats as the Fortran program's: the state A5, Time and Voltage F12.5, Step I15, the others ES15.5."""
    fixed = {"Time", "Voltage"}

    def head(cols):
        return " ".join(f"{cols[i][:5]:>5}" if i == 0 else (f"{cols[i][:12]:>12}" if header[i] in fixed
                        else f"{cols[i][:15]:>15}") for i in keep)
    fh.write(head(header) + "\n" + head(units) + "\n")
    for r in rows:
        out = []
        for i in keep:
            v = r[i]
            if i == 0:
                out.append(f"{v:>5}")
            elif header[i] in fixed:
                out.append(f"{'NaN':>12}" if math.isnan(v) else f"{v:12.5f}")
            elif isinstance(v, int):
                out.append(f"{v:15d}")
            else:
                out.append(_sci(v))
        fh.write(" ".join(out) + "\n")


def write_profiles(path, var_names, profiles) -> None:
    """profiles: [(t [s], x [um] array, {name: array})]."""
    cols = ["Time", "x"] + list(var_names)
    units = ["hours", "um"] + [PROFILE_UNITS[v] for v in var_names]
    with open(path, "w") as fh:
        fh.write(" ".join(f"{c[:15]:>15}" for c in cols) + "\n")
        fh.write(" ".join(f"{u[:15]:>15}" for u in units) + "\n")
        for t, x, vals in profiles:
            for j in range(len(x)):
                fh.write(" ".join([_sci(t / 3600.0, 15, 7), _sci(float(x[j]), 15, 7)]
                                  + [_sci(float(vals[v][j]), 15, 7) for v in var_names]) + "\n")


class ProfileRecorder:
    """The profile schedule of a run: with every p.profile_every-th output row, at the first row at or after each of
    p.profile_times, and with the last row. `profile(state, I)` gives (x [um], {variable: values})."""

    def __init__(self, p, available, profile):
        self.on = bool(p.profile_file)
        self.every, self.profile = p.profile_every, profile
        self.times = times(p.profile_times) if self.on else []
        self.vars = choose(available, p.profile_vars, "profile variable") if self.on else []
        self.profiles, self.n_rows, self.last = [], 0, -1

    def row(self, t, state, I):
        if not self.on:
            return
        i = self.n_rows
        self.n_rows += 1
        due = self.every > 0 and i % self.every == 0
        while self.times and t >= self.times[0] * (1.0 - 1.0e-12):
            self.times.pop(0)
            due = True
        if due:
            self._take(t, state, I)
            self.last = i

    def finish(self, res, path):
        """The last row's profile, if not already taken; write the file."""
        if not self.on:
            return
        if self.last != self.n_rows - 1:
            t, _, I = res.end
            self._take(t, res.final_state, I)
        res.profiles, res.profile_vars = self.profiles, tuple(self.vars)
        write_profiles(path, self.vars, self.profiles)

    def _take(self, t, state, I):
        x, vals = self.profile(state, I)
        self.profiles.append((t, x, {v: vals[v] for v in self.vars}))
