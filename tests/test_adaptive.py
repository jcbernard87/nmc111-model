"""Adaptive time stepping (dt_mode = 'adaptive', #8): backward Euler with the local error estimated from a
linear predictor and controlled against rtol/atol, in both particle models; steps land on the output times."""
import subprocess

import numpy as np
import pytest

from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import run_corrected
from nmc_model.uniform.params import Params
from nmc_model.uniform.simulate import run as run_uniform

from .test_cross_language import AGG, AGG_CYCLE, agg_groups, assert_same

S = "cc C=2 Vmin=3.0; rest t=600"
AG = dict(nj=30, nja=10)


def run_model(model, **kw):
    kw.setdefault("steps", S)
    if model == "uniform":
        return run_uniform(Params(mode="corrected", **kw))
    return run_corrected(AggParams(mode="corrected", **AG, **kw))


_REF = {}


def reference(model):
    """Fixed steps of 0.25 s; rows every 18 s (on the same times as adaptive rows)."""
    if model not in _REF:
        _REF[model] = (run_model("uniform", t_max=3600, n_steps=14400) if model == "uniform"
                       else run_model("agglomerate", dt_s=0.25))
    return _REF[model]


def discharge_error(a, ref):
    """Largest voltage difference [V] at the rows both runs write during the discharge (step 1)."""
    ta = {round(t * 3600, 6): v for t, v, k in a[:, [0, 1, 7]] if k == 1}
    tr = {round(t * 3600, 6): v for t, v, k in ref[:, [0, 1, 7]] if k == 1}
    common = sorted(set(ta) & set(tr))[:-1]                  # the cutoff rows are located, not on the grid
    assert len(common) > 50
    return max(abs(ta[t] - tr[t]) for t in common)


def test_fixed_is_the_default():
    assert Params().dt_mode == AggParams().dt_mode == "fixed"
    for cls in (Params, AggParams):
        with pytest.raises(ValueError, match="dt_mode must be one of"):
            cls(dt_mode="variable")
        with pytest.raises(ValueError, match="needs mode='corrected'"):
            cls(dt_mode="adaptive", mode="faithful")


@pytest.mark.parametrize("model", ["uniform", "agglomerate"])
def test_adaptive_converges_with_rtol_and_takes_fewer_steps(model):
    ref = reference(model)
    loose = run_model(model, dt_mode="adaptive", rtol=1e-3, atol=1e-3)
    tight = run_model(model, dt_mode="adaptive", rtol=1e-5, atol=1e-5)
    e_loose, e_tight = discharge_error(loose.array, ref.array), discharge_error(tight.array, ref.array)
    assert e_tight < 1e-5 and e_tight < e_loose / 5           # [V]
    assert loose.steps < tight.steps < ref.steps / 4
    assert loose.exit_reason == tight.exit_reason == ref.exit_reason == "end_of_protocol"
    gap = np.diff(tight.array[:, 0] * 3600)                   # rows write_interval apart, but cutoff and end
    assert np.sum(np.abs(gap - 18.0) > 1e-9) <= 2


@pytest.mark.parametrize("steps", [S, AGG_CYCLE], ids=["discharge", "cccv"])
def test_adaptive_agglomerate_identical_in_every_language(fortran_exe, cpp_exe, run_native, tmp_path, steps):
    """The agglomerate model is byte-identical in Python, Fortran and C++ (Li_Nernst to round-off), and so are
    its adaptive steps: the error norm is summed in the same order in all three."""
    g = agg_groups(steps=steps)
    g["numerics"] = dict(g["numerics"], dt_mode="adaptive")
    f = run_native(fortran_exe, g, name="f.txt")
    c = run_native(cpp_exe, g, name="c.txt")
    r = run_corrected(AggParams(mode="corrected", steps=steps, dt_mode="adaptive", **AGG))
    py = tmp_path / "p.txt"
    r.write(py)
    assert_same(f, c)
    assert_same(f, py)


@pytest.mark.parametrize("steps", [S, AGG_CYCLE.replace("4.3", "4.2")], ids=["discharge", "cccv"])
def test_adaptive_uniform_agrees_in_every_language(fortran_exe, cpp_exe, run_native, tmp_path, steps):
    """Fortran and C++ are identical; Python agrees to the printed precision, with the same rows."""
    g = {"model": {"particle_model": "uniform"}, "protocol": {"steps": steps},
         "numerics": {"mode": "corrected", "dt_mode": "adaptive"}}
    f = run_native(fortran_exe, g, name="f.txt")
    c = run_native(cpp_exe, g, name="c.txt")
    assert_same(f, c)
    fa = np.array([[float(x) for x in ln.split()[1:]] for ln in f.read_text().splitlines()[2:]])
    pa = run_model("uniform", steps=steps, dt_mode="adaptive").array
    assert fa.shape[0] == pa.shape[0]
    np.testing.assert_allclose(fa[:, :3], pa[:, :3], rtol=0, atol=1e-5)          # 5 printed decimals


@pytest.mark.parametrize("which", ["fortran", "cpp"])
@pytest.mark.parametrize("numerics, msg", [
    ({"mode": "corrected", "dt_mode": "variable"}, "dt_mode must be one of ('fixed', 'adaptive'), got 'variable'"),
    ({"mode": "faithful", "dt_mode": "adaptive"}, "dt_mode='adaptive' needs mode='corrected'"),
])
def test_native_programs_check_dt_mode(fortran_exe, cpp_exe, run_native, which, numerics, msg):
    exe = fortran_exe if which == "fortran" else cpp_exe
    with pytest.raises(subprocess.CalledProcessError) as e:
        run_native(exe, {"model": {"particle_model": "uniform"}, "numerics": numerics})
    assert msg in e.value.stderr.decode()
