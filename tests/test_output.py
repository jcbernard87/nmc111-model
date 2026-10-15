"""Configurable output (#12): chosen time-series columns, and profiles against x in a long table."""
import subprocess

import numpy as np
import pytest

from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import run_corrected
from nmc_model.uniform.params import Params
from nmc_model.uniform.simulate import run as run_uniform

STEPS = "cc C=2 t=600; rest t=300"
AGG = dict(nj=30, nja=10, dt_s=5.0)
COLUMNS = {"uniform": "Voltage, Current", "crystal": "Voltage, Th_mean_sep", "agglomerate": "Voltage, x_front"}


def _run(model, **kw):
    if model == "agglomerate":
        return run_corrected(AggParams(mode="corrected", **AGG, **kw))
    return run_uniform(Params(mode="corrected", particle_model=model, **kw))


def _groups(model, **output):
    g = {"model": {"particle_model": model}, "numerics": {"mode": "corrected"}, "protocol": {"steps": STEPS},
         "output": output}
    if model == "agglomerate":
        g["cell"] = {"nj": AGG["nj"]}
        g["agglomerate"] = {"nja": AGG["nja"], "dt_s": AGG["dt_s"]}
    return g


def _table(path):
    lines = path.read_text().splitlines()
    return lines[:2], np.array([[float(v) for v in ln.split()] for ln in lines[2:]])


@pytest.mark.parametrize("model", ["uniform", "crystal", "agglomerate"])
def test_profiles_and_columns_in_every_language(fortran_exe, cpp_exe, run_native, tmp_path, model):
    out = dict(profile_every=7, profile_times="100 450", columns=COLUMNS[model])
    files = {}
    for exe, tag in ((fortran_exe, "f"), (cpp_exe, "c")):
        pf = tmp_path / f"{tag}_prof.txt"
        files[tag] = (run_native(exe, _groups(model, profile_file=str(pf), **out), name=f"{tag}.txt"), pf)
    assert files["f"][0].read_bytes() == files["c"][0].read_bytes()       # the time series identical
    hf, af = _table(files["f"][1])
    hc, ac = _table(files["c"][1])
    assert hf == hc and af.shape == ac.shape                              # profiles to round-off (Li_Nernst)
    np.testing.assert_allclose(af, ac, rtol=1e-12, atol=1e-15, equal_nan=True)
    pp = tmp_path / "p_prof.txt"
    r = _run(model, steps=STEPS, profile_file=str(pp), **out)
    r.write(tmp_path / "p.txt")
    assert files["f"][0].read_text().splitlines()[:2] == (tmp_path / "p.txt").read_text().splitlines()[:2]
    hp, ap = _table(pp)
    assert hp == hf and ap.shape == af.shape
    np.testing.assert_allclose(ap, af, rtol=1e-6, atol=1e-12, equal_nan=True)


def test_default_output_is_unchanged(tmp_path):
    r = _run("uniform", steps="cc C=1 t=60")
    r.write(tmp_path / "o.txt")
    assert (tmp_path / "o.txt").read_text().splitlines()[0].split()[-3:] == ["Current", "Step", "Li_Nernst"]
    assert r.profiles == []


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
@pytest.mark.parametrize("kw,msg", [(dict(columns="Voltage, x_front"), "unknown output column 'x_front'"),
                                    (dict(profile_vars="c, theta_surf"), "unknown profile variable 'theta_surf'")])
def test_unknown_names_are_rejected(fortran_exe, cpp_exe, run_native, tmp_path, which, kw, msg):
    """x_front is the agglomerate model's column, theta_surf the crystal and agglomerate models' variable."""
    kw = dict(kw, profile_file=str(tmp_path / "pr.txt")) if "profile_vars" in kw else dict(kw)
    if which == "python":
        with pytest.raises(ValueError, match=msg):
            _run("uniform", steps="cc C=1 t=60", **kw).write(tmp_path / "o.txt")
        return
    exe = fortran_exe if which == "fortran" else cpp_exe
    with pytest.raises(subprocess.CalledProcessError) as e:
        run_native(exe, _groups("uniform", **kw))
    assert msg in e.value.stderr.decode(), e.value.stderr.decode()
