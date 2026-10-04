"""The input files, the Python reader and the compiled programs' defaults agree."""
import re
import subprocess
from dataclasses import replace
from pathlib import Path

import pytest

from nmc_model.agglomerate import AggParams
from nmc_model.inputs import load, run_file
from nmc_model.uniform.params import Params

ROOT = Path(__file__).resolve().parents[1]
SHORT = "cc C=1 t=120"   # a short run keeps the comparisons quick


def test_input_files_hold_the_defaults():
    p, model, out = load(ROOT / "input" / "uniform.nml")
    assert model == "uniform" and out == "Time_Voltage.txt"
    assert p == Params(mode="corrected")
    p, model, _ = load(ROOT / "input" / "agglomerate.nml")
    assert model == "agglomerate"
    assert p == replace(AggParams(mode="corrected"), tau_agg=-1.0, sigma_agg=-1.0)


def _short(text, name):
    text, n = re.subn(r"steps = '[^']*'", f"steps = '{SHORT}'", text)
    assert n == 1
    return text.replace("file = 'Time_Voltage.txt'", f"file = '{name}'")


@pytest.mark.parametrize("model", ["uniform", "agglomerate"])
def test_programs_default_to_the_input_file(model, fortran_exe, cpp_exe, tmp_path):
    """A minimal input (model and mode only) gives the same run as the full input file, in every language."""
    full = _short((ROOT / "input" / f"{model}.nml").read_text(), "full.txt")
    minimal = (f"&model\n  particle_model = '{model}'\n/\n&numerics\n  mode = 'corrected'\n/\n"
               f"&protocol\n  steps = '{SHORT}'\n/\n&output\n  file = 'min.txt'\n/\n")
    (tmp_path / "full.nml").write_text(full)
    (tmp_path / "min.nml").write_text(minimal)
    outputs = []
    for exe in (fortran_exe, cpp_exe):
        for nml, name in (("full.nml", "full.txt"), ("min.nml", "min.txt")):
            subprocess.run([str(exe), nml], cwd=tmp_path, check=True, capture_output=True)
            outputs.append((tmp_path / name).read_bytes())
    run_file(tmp_path / "full.nml", out=tmp_path / "py.txt")
    outputs.append((tmp_path / "py.txt").read_bytes())
    if model == "agglomerate":
        assert all(o == outputs[0] for o in outputs)
    else:
        assert all(o == outputs[0] for o in outputs[:4])   # Python agrees to rtol 1e-5 (test_cross_language)


@pytest.mark.parametrize("model,missing", [("uniform", "k_rxn and sigma"), ("agglomerate", "D_agg")])
def test_faithful_needs_the_fitted_values(model, missing, fortran_exe, cpp_exe, tmp_path):
    nml = tmp_path / "f.nml"
    nml.write_text(f"&model\n  particle_model = '{model}'\n/\n&numerics\n  mode = 'faithful'\n/\n")
    for exe in (fortran_exe, cpp_exe):
        r = subprocess.run([str(exe), str(nml)], cwd=tmp_path, capture_output=True, text=True)
        assert r.returncode != 0 and missing in (r.stdout + r.stderr)
    with pytest.raises(ValueError):
        load(nml)


BAD = [
    ("&model\n  particle_model = 'agglomerat'\n/\n", "particle_model must be 'uniform' or 'agglomerate'"),
    ("&numerics\n  mode = 'correctd'\n/\n", "mode must be 'faithful' or 'corrected'"),
    ("&protocol\n  steps = 'cx C=1'\n/\n", "unknown step type"),
    ("&cell\n  nj = 30, sep_node = 2\n/\n", "sep_node must be at least 3"),
    ("&cell\n  nj = 24, sep_node = 22\n/\n", "nj - sep_node must be at least 3"),
    ("&model\n  particle_model = 'agglomerate'\n/\n&agglomerate\n  nja = 3\n/\n", "nja must be at least 4"),
    ("&initial\n  cs_init = 0.01\n/\n", "unknown namelist group &initial"),
    ("&cell\n  cs_init = 0.01\n/\n", ("cs_init", "&cell")),          # a name in the wrong group
]


@pytest.mark.parametrize("text,msg", BAD)
def test_programs_reject_bad_values(text, msg, fortran_exe, cpp_exe, tmp_path):
    """A bad option value, mesh size, group or misplaced name stops each program with a non-zero exit and that
    specific message (not just the key's name, which a generic 'unknown name' error would also contain) (#5).
    The Fortran program skipped unknown groups silently, and the C++ program accepted a name in any group."""
    if "mode = " not in text:
        text += "&numerics\n  mode = 'corrected'\n/\n"
    nml = tmp_path / "bad.nml"
    nml.write_text(text + "&output\n  file = 'out.txt'\n/\n")
    msgs = msg if isinstance(msg, tuple) else (msg,)
    for exe in (fortran_exe, cpp_exe):
        r = subprocess.run([str(exe), str(nml)], cwd=tmp_path, capture_output=True, text=True)
        assert r.returncode != 0 and all(m in r.stdout + r.stderr for m in msgs), (exe, r.stdout + r.stderr)
    with pytest.raises(ValueError) as err:
        p, model, _ = load(nml)
        from nmc_model.protocol import parse
        parse(p.steps, p)                           # the protocol is parsed when a run starts
    assert all(m in str(err.value) for m in msgs), str(err.value)
