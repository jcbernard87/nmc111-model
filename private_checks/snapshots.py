"""Read the run-specific values the original run generator substituted into an archived snapshot.

Used only by the private checks: the snapshot files are part of the private oracle and are
never copied into this repository.
"""
import re
from pathlib import Path

from nmc_model.agglomerate import AggParams


def _active(text):
    """Source lines with Fortran comments removed."""
    for line in text.splitlines():
        code = line.split("!", 1)[0].strip()
        if code:
            yield code


def agg_run_values(snapshot: Path) -> dict:
    text = Path(snapshot).read_text()
    vals = {}
    for code in _active(text):
        m = re.match(r"parameter\(N=4,NJ=([\d.]+),N_c=4,NJ_c=([\d.]+)", code)
        if m:
            vals["nj"], vals["nja"] = int(float(m.group(1))), int(float(m.group(2)))
        m = re.match(r"double precision :: Time_mod = ([\d.]+)$", code)
        if m:
            vals["time_mod"] = float(m.group(1))
        for key, pat in (("C_rate", r"C_rate = ([-\d.eE+]+)$"), ("eps", r"eps = ([-\d.eE+]+)$"),
                         ("sigma", r"sigma = ([-\d.eE+]+)$")):
            m = re.match(pat, code)
            if m:
                vals[key] = float(m.group(1))
        m = re.match(r"xmax = (\d+)/10000\.0$", code)
        if m:
            vals["thickness_um"] = float(m.group(1))
        # the original's fitted values (never distributed with the repository)
        for key, pat in (("mass_loading", r"Massloading = ([\d.]+)$"), ("tortuosity", r"tortuosity = ([\d.]+)$"),
                         ("percent_active", r"percent_active = ([\d.]+)$"),
                         ("D_agg_mult", r"diff_agg = 1\.0d-10 \* ([\d.]+)$"),
                         ("k_exp", r"rxn_k = +10\.0\*\*\(-1 \* ([\d.]+)\)$")):
            m = re.match(pat, code)
            if m:
                vals[key] = float(m.group(1))
    missing = {"nj", "nja", "time_mod", "C_rate", "eps", "sigma", "thickness_um", "mass_loading", "tortuosity",
               "percent_active", "D_agg_mult", "k_exp"} - vals.keys()
    if missing:
        raise ValueError(f"{snapshot}: could not read {sorted(missing)}")
    return vals


def uniform_fit_values(oracle_dir: Path) -> dict:
    """The uniform-particle model's fitted values, read from the original source (never distributed)."""
    text = (Path(oracle_dir) / "source" / "crystal" / "questions_NMC111_Electrode_Crystal.f95").read_text()
    vals = {}
    for line in text.splitlines():
        code = line.split("!", 1)[0].strip()
        m = re.search(r"rxn_k = 1\.0d-8\*10\*\*\(([\d.]+)\)$", code)
        if m:
            vals["k_exp"] = float(m.group(1))
        m = re.match(r"sigma=([\d.]+)d([-\d]+)$", code)
        if m:
            vals["sigma"] = float(f"{m.group(1)}e{m.group(2)}")
    if set(vals) != {"k_exp", "sigma"}:
        raise ValueError(f"could not read the uniform model's fitted values: {vals}")
    return vals


def uniform_namelist(oracle_dir: Path, out_file: str = "Time_Voltage.txt") -> str:
    """Input file for the compiled programs reproducing the uniform-particle model's archived run."""
    from nmc_model.uniform.params import Params
    p = Params.faithful(C_rate=1.0, **uniform_fit_values(oracle_dir))
    return f"""&model
  particle_model = 'uniform'
/
&cell
  L_cath_um = 24, L_sep = 25.0d-4, nj = 43, sep_node = 22
  eps = 0.5, eps_AM = 0.4, eps_sep = 0.39, tau_sep = 4.8, bruggeman = -0.5
/
&electrolyte
  D = 2.0d-6, t_plus = 0.25, c_bulk = 1.0d-3
/
&active
  sigma = {p.sigma!r}, k_rxn = {p.k_rxn!r}, M = 96.46, rho = 4.6, Q_th = 0.150, R_p = 200.0d-7
  alpha_a = 0.5, alpha_c = 0.5
/
&operation
  C_rate = 1.0
  phi1_init = 4.1, phi2_init = 0.0, cs_init = 1.0d-5
  t_max = 36000.0, n_steps = 36000
/
&numerics
  mode = 'faithful'
/
&output
  file = '{out_file}'
/
"""


def agg_params(snapshot: Path) -> AggParams:
    return AggParams.faithful(**agg_run_values(snapshot))


def agg_namelist(snapshot: Path, out_file: str = "Time_Voltage.txt") -> str:
    """Input file for the compiled programs reproducing an archived agglomerate snapshot.

    D_agg and k_rxn are passed as the exact doubles the original computed in single precision.
    """
    v = agg_run_values(snapshot)
    p = agg_params(snapshot)
    return f"""&model
  particle_model = 'agglomerate'
/
&cell
  nj = {v['nj']}, L_cath_um = {v['thickness_um']}, eps = {v['eps']!r}
/
&active
  sigma = {v['sigma']!r}, k_rxn = {p.rxn_k!r}
/
&operation
  C_rate = {v['C_rate']!r}
/
&agglomerate
  nja = {v['nja']}, time_mod = {v['time_mod']!r}, D_agg = {p.D_agg!r}, tortuosity_e = {v['tortuosity']!r},
  mass_loading = {v['mass_loading']!r}, percent_active = {v['percent_active']!r}
/
&numerics
  mode = 'faithful'
/
&output
  file = '{out_file}'
/
"""
