"""Restart files (#13): the state, time, capacity and current at the end of a run, to start another run from.

Plain text, read and written alike by the Python package and the Fortran and C++ programs:

    # nmc111-model restart 1
    model <particle model>
    t <time [s]>
    mAhg <capacity [mAh/g]>
    I <last current [A/cm2]>
    array <name> <rows> <cols>
    <rows*cols values, row-major, one per line, 17 significant digits>
    ...

17 significant digits round-trip every double exactly, so a run continued from a restart repeats the run it
continues. The arrays are the electrode state `c` (nj, 4) and, for the crystal and agglomerate models, `xc`
(crystal log-odds (nl, nj_crystal), or the agglomerates' unknowns (nl, nja*4)).
"""
from __future__ import annotations

import numpy as np

HEADER = "# nmc111-model restart 1"


def write(path, model: str, t: float, mAhg: float, I: float, arrays) -> None:
    lines = [HEADER, f"model {model}", f"t {t:.16e}", f"mAhg {mAhg:.16e}", f"I {I:.16e}"]
    for name, a in arrays:
        a = np.asarray(a, dtype=float)
        lines.append(f"array {name} {a.shape[0]} {a.shape[1]}")
        lines.extend(f"{v:.16e}" for v in a.ravel())
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def read(path, model: str, shapes: dict) -> dict:
    """The restart's t, mAhg, I and arrays; `shapes` maps each array name to the shape this input needs."""
    with open(path) as fh:
        tok = fh.read().split("\n")
    if not tok or tok[0].strip() != HEADER:
        raise ValueError(f"{path}: not a restart file of this model")
    out, arrays, i = {}, {}, 1
    while i < len(tok):
        words = tok[i].split()
        i += 1
        if not words:
            continue
        if words[0] == "array":
            name, n, m = words[1], int(words[2]), int(words[3])
            arrays[name] = np.array([float(v) for v in tok[i:i + n * m]]).reshape(n, m)
            i += n * m
        else:
            out[words[0]] = words[1]
    if out.get("model") != model:
        raise ValueError(f"restart file does not match the input: model {out.get('model')}, input {model}")
    for name, shape in shapes.items():
        if name not in arrays or arrays[name].shape != tuple(shape):
            got = arrays[name].shape if name in arrays else "missing"
            raise ValueError(f"restart file does not match the input: array {name} {got}, input {tuple(shape)}")
    return dict(t=float(out["t"]), mAhg=float(out["mAhg"]), I=float(out["I"]), arrays=arrays)
