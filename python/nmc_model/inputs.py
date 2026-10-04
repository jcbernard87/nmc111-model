"""Read the shared input file for either particle model, and run it.

The Fortran and C++ programs read the same file. `&model particle_model` chooses the model:
'uniform' (the default) or 'agglomerate'. See input/*.nml and docs/parameters.md.
"""
from __future__ import annotations

from dataclasses import fields
from pathlib import Path

from .agglomerate import AggParams
from .uniform import namelist as _unl

# namelist name -> AggParams field, per group
AGG_NAMES = {
    "cell": {"l_cath_um": "thickness_um", "l_sep": "L_sep", "nj": "nj", "sep_node": "sep_node", "eps": "eps",
             "eps_sep": "eps_sep", "tau_sep": "tau_sep", "bruggeman": "bruggeman", "eps_am": None},
    "electrolyte": {"d": "D0", "t_plus": "t_plus", "c_bulk": "c_bulk", "kappa_bg": "kappa_bg", "z_plus": None,
                    "z_minus": None},
    "active": {"sigma": "sigma", "m": "M", "rho": "rho", "q_th": "Q_th", "k_rxn": "rxn_k", "alpha_a": "alpha_a",
               "alpha_c": "alpha_c", "k_li": "k_Li", "c_li_ref": "c_Li_ref", "r_p": None},
    "constants": {"r": "R", "t": "T", "f": "F"},
    "operation": {"c_rate": "C_rate", "phi1_init": "phi1_init", "cs_init": "cs_init", "t_max": "t_max",
                  "v_min": "V_min", "v_max": "V_max", "phi2_init": None, "n_steps": None},
    "numerics": {"fd_step": "fd_step", "newton_tol": "newton_tol", "newton_max_iter": "newton_max_iter",
                 "mode": "mode"},
    "protocol": {"steps": "steps", "cycles": "cycles"},
    "output": {"file": None, "write_interval": "write_interval"},
    "agglomerate": {"nja": "nja", "time_mod": "time_mod", "r_agg": "R_agg", "r_xtal": "R_xtal",
                    "eps_agg": "eps_agg", "d_agg": "D_agg", "tortuosity_e": "tortuosity",
                    "mass_loading": "mass_loading", "percent_active": "percent_active", "mol_vol": "mol_vol",
                    "c0_init": "c0_init", "tau_agg": "tau_agg", "sigma_agg": "sigma_agg", "dt_s": "dt_s"},
}


def _check_mesh(p, model):
    """The separator and the cathode each need interior nodes, an agglomerate at least four."""
    if p.sep_node < 3:
        raise ValueError("sep_node must be at least 3")
    if p.nj - p.sep_node < 3:
        raise ValueError("nj - sep_node must be at least 3")
    if model == "agglomerate" and p.nja < 4:
        raise ValueError("nja must be at least 4")


def load(path):
    """Read an input file. Returns (params, particle_model, output_file)."""
    text = Path(path).read_text()
    data = _unl.parse(text)
    model = str(data.pop("model", {}).get("particle_model", "uniform")).lower()
    out_file = data.get("output", {}).get("file", "Time_Voltage.txt")
    if model == "uniform":
        if "agglomerate" in data:
            raise ValueError("&agglomerate is only for particle_model = 'agglomerate'")
        p, extra = _unl.from_dict(data)
        _check_mesh(p, model)
        return p, model, extra.get("file", out_file)
    if model != "agglomerate":
        raise ValueError("particle_model must be 'uniform' or 'agglomerate'")
    values = {}
    for group, entries in data.items():
        if group not in AGG_NAMES:
            raise ValueError(f"unknown namelist group &{group}")
        for key, val in entries.items():
            if key not in AGG_NAMES[group]:
                raise ValueError(f"unknown name {key!r} in &{group}")
            name = AGG_NAMES[group][key]
            if name is not None:
                values[name] = val
    types = {f.name: f.type for f in fields(AggParams)}
    for name, val in list(values.items()):
        t = str(types[name])
        values[name] = int(val) if t == "int" else (str(val) if t == "str" else float(val))
    mode = values.pop("mode", "corrected")
    if mode not in ("faithful", "corrected"):
        raise ValueError(f"mode must be 'faithful' or 'corrected', got {mode!r}")
    if mode == "faithful":
        need = {"D_agg", "rxn_k", "tortuosity", "mass_loading"}
        if not need <= values.keys():
            raise ValueError("faithful agglomerate runs need D_agg, k_rxn, tortuosity_e and mass_loading "
                             "(the original's fitted values are not distributed)")
        p = AggParams.faithful(**values)
    else:
        p = AggParams(mode=mode, **values)
    _check_mesh(p, model)
    return p, model, out_file


def run_file(path, out=None):
    """Run the model an input file describes and write its output file. Returns the result."""
    p, model, out_file = load(path)
    if model == "uniform":
        from .uniform.simulate import run
        r = run(p)
    elif p.mode == "faithful":
        from .agglomerate import run
        r = run(p)
    else:
        from .agglomerate_corrected import run_corrected
        r = run_corrected(p)
    r.write(out or Path(path).parent / out_file)
    return r
