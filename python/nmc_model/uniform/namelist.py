"""Reader for the shared input file (a Fortran namelist subset).

The Fortran program reads the same file with native NAMELIST input, and the C++
program uses an equivalent parser. Supported syntax:

    ! comment
    &group
      name = 1.0d-3, other = 22   ! comments after values
      mode = 'faithful'
    /

Values are integers, reals (e/E/d/D exponents), quoted strings, or .true./.false.
Names are case-insensitive. Arrays and repeat counts are not supported.
"""
from __future__ import annotations

import re
from dataclasses import fields

from .params import Params, f32

GROUPS = {
    "cell": ["L_cath_um", "L_sep", "nj", "sep_node", "eps", "eps_AM", "eps_sep", "tau_sep", "bruggeman"],
    "electrolyte": ["D", "t_plus", "c_bulk", "z_plus", "z_minus", "kappa_bg"],
    "active": ["sigma", "M", "rho", "Q_th", "R_p", "k_rxn", "alpha_a", "alpha_c", "k_Li", "c_Li_ref"],
    "constants": ["R", "T", "F"],
    "operation": ["C_rate", "phi1_init", "phi2_init", "cs_init", "t_max", "n_steps", "V_min", "V_max"],
    "numerics": ["fd_step", "newton_tol", "newton_max_iter", "mode"],
    "protocol": ["steps", "cycles"],
    "output": ["file", "write_interval"],
}

# Values the original program stored from single-precision literals (deviation D-4).
# In faithful mode, values read from the input file for these names are rounded to float32.
FAITHFUL_F32 = {"R", "c_bulk", "Q_th", "M", "rho", "phi1_init", "eps_sep", "eps_AM", "tau_sep", "C_rate"}

_TOKEN = re.compile(r"""\s*(?:'([^']*)'|"([^"]*)"|([^,\s/=]+))""")


def _strip_comment(line: str) -> str:
    out, quote = [], None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            break
        out.append(ch)
    return "".join(out)


def _value(text: str, quoted: bool):
    if quoted:
        return text
    low = text.lower()
    if low in (".true.", "t", ".t."):
        return True
    if low in (".false.", "f", ".f."):
        return False
    try:
        return int(text)
    except ValueError:
        return float(low.replace("d", "e"))


def parse(text: str) -> dict:
    """Parse namelist text into {group: {name: value}} with lower-case keys."""
    body = "\n".join(_strip_comment(l) for l in text.splitlines())
    out: dict = {}
    for m in re.finditer(r"&(\w+)(.*?)/", body, flags=re.S):
        group, content = m.group(1).lower(), m.group(2)
        entries = {}
        for part in re.finditer(r"(\w+)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s/]+)", content):
            name, raw = part.group(1).lower(), part.group(2)
            quoted = raw[:1] in "'\""
            entries[name] = _value(raw[1:-1] if quoted else raw, quoted)
        if group in out:
            raise ValueError(f"namelist group &{group} appears twice")
        out[group] = entries
    return out


def load(path) -> tuple[Params, dict]:
    """Read an input file; return (Params, extra) where extra holds non-model settings (output file)."""
    with open(path) as fh:
        return from_dict(parse(fh.read()))


def from_dict(data: dict) -> tuple[Params, dict]:
    """Parsed namelist groups -> (Params, extra)."""
    known = {g: {n.lower(): n for n in names} for g, names in GROUPS.items()}
    values, extra = {}, {}
    for group, entries in data.items():
        if group not in known:
            raise ValueError(f"unknown namelist group &{group}")
        for key, val in entries.items():
            if key not in known[group]:
                raise ValueError(f"unknown name {key!r} in &{group}")
            name = known[group][key]
            (extra if name == "file" else values)[name] = val
    mode = values.pop("mode", "faithful")
    types = {f.name: f.type for f in fields(Params)}
    for name, val in list(values.items()):
        if types[name] in ("int", int):
            values[name] = int(val)
        elif types[name] in ("str", str):
            values[name] = str(val)
        else:
            values[name] = float(val)
    if mode == "faithful":
        if "k_rxn" not in values or "sigma" not in values:
            raise ValueError("faithful runs need k_rxn and sigma (the original's fitted values are not distributed)")
        for name in FAITHFUL_F32 & values.keys():
            values[name] = f32(values[name])
        base = Params.faithful(k_rxn=values["k_rxn"], sigma=values["sigma"])
        p = base.with_(**values)
    else:
        p = Params(mode=mode, **values)
    return p, extra
