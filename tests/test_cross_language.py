"""Fortran, C++ and Python implementations agree (corrected mode).

The corrected agglomerate model does the same operations in the same order in all three
languages, so the written files are identical (apart from round-off in the Li_Nernst column
during rests; see assert_same). The uniform model's Python Newton step uses
numpy expressions that round differently, so it agrees with the compiled programs to the
printed precision's last digit.
"""
import numpy as np
import pytest

from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import run_corrected
from nmc_model.uniform.params import Params
from nmc_model.uniform.simulate import run as run_uniform

AGG = dict(nj=30, nja=10, dt_s=5.0)
AGG_CYCLE = "cc C=1 Vmin=3.0; rest t=600; cc C=-1 Vmax=4.3; cv V=4.3 Imin=0.05; rest t=600"


LI_NERNST = 9   # token index of the Li_Nernst column


def assert_same(a, b):
    """Files identical, except that Li_Nernst = (RT/F) ln(c/c_ref) may differ in round-off.

    While c at the foil relaxes back to c_ref during a rest, that column is round-off
    (about 1e-10 mV) and its printed digits differ between compilers.
    """
    la, lb = a.read_text().splitlines(), b.read_text().splitlines()
    assert len(la) == len(lb) and la[:2] == lb[:2]
    for x, y in zip(la[2:], lb[2:]):
        tx, ty = x.split(), y.split()
        assert tx[:LI_NERNST] + tx[LI_NERNST + 1:] == ty[:LI_NERNST] + ty[LI_NERNST + 1:], (x, y)
        assert abs(float(tx[LI_NERNST]) - float(ty[LI_NERNST])) <= 1e-6 + 1e-4 * abs(float(ty[LI_NERNST]))


def agg_groups(C_rate=1.0, steps=""):
    g = {"model": {"particle_model": "agglomerate"},
         "cell": {"nj": AGG["nj"]},
         "agglomerate": {"nja": AGG["nja"], "dt_s": AGG["dt_s"]},
         "operation": {"C_rate": C_rate},
         "numerics": {"mode": "corrected"}}
    if steps:
        g["protocol"] = {"steps": steps}
    return g


def py_agg(tmp_path, name, C_rate=1.0, steps=""):
    r = run_corrected(AggParams(mode="corrected", C_rate=C_rate, steps=steps, **AGG))
    out = tmp_path / name
    r.write(out)
    return out


@pytest.mark.parametrize("C_rate", [1.0, 2.0])
def test_agglomerate_fortran_identical_to_python(fortran_exe, run_native, tmp_path, C_rate):
    f = run_native(fortran_exe, agg_groups(C_rate), name=f"f{C_rate}.txt")
    assert f.read_bytes() == py_agg(tmp_path, f"p{C_rate}.txt", C_rate).read_bytes()


def test_agglomerate_cycle_fortran_identical_to_python(fortran_exe, run_native, tmp_path):
    f = run_native(fortran_exe, agg_groups(steps=AGG_CYCLE), name="fcyc.txt")
    assert_same(f, py_agg(tmp_path, "pcyc.txt", steps=AGG_CYCLE))


def read(path):
    lines = path.read_text().splitlines()[2:]
    return [l.split()[0] for l in lines], np.array([[float(x) for x in l.split()[1:]] for l in lines])


@pytest.mark.parametrize("C_rate", [1.0, 2.0])
def test_uniform_fortran_matches_python(fortran_exe, run_native, tmp_path, C_rate):
    f = run_native(fortran_exe, {"model": {"particle_model": "uniform"}, "operation": {"C_rate": C_rate},
                                 "numerics": {"mode": "corrected"}}, name=f"u{C_rate}.txt")
    r = run_uniform(Params(mode="corrected", C_rate=C_rate))
    p = tmp_path / f"pu{C_rate}.txt"
    r.write(p)
    sf, vf = read(f)
    sp, vp = read(p)
    assert sf == sp and vf.shape == vp.shape
    np.testing.assert_allclose(vf, vp, rtol=1e-5, atol=1e-9)


@pytest.mark.parametrize("C_rate", [1.0, 2.0])
def test_agglomerate_cpp_identical_to_python(cpp_exe, run_native, tmp_path, C_rate):
    c = run_native(cpp_exe, agg_groups(C_rate), name=f"c{C_rate}.txt")
    assert c.read_bytes() == py_agg(tmp_path, f"pc{C_rate}.txt", C_rate).read_bytes()


def test_agglomerate_cycle_cpp_identical_to_fortran(cpp_exe, fortran_exe, run_native):
    f = run_native(fortran_exe, agg_groups(steps=AGG_CYCLE), name="fcyc2.txt")
    c = run_native(cpp_exe, agg_groups(steps=AGG_CYCLE), name="ccyc2.txt")
    assert_same(c, f)


@pytest.mark.parametrize("C_rate", [1.0, 2.0])
def test_uniform_cpp_identical_to_fortran(cpp_exe, fortran_exe, run_native, C_rate):
    g = {"model": {"particle_model": "uniform"}, "operation": {"C_rate": C_rate}, "numerics": {"mode": "corrected"}}
    f = run_native(fortran_exe, g, name=f"uf{C_rate}.txt")
    c = run_native(cpp_exe, g, name=f"uc{C_rate}.txt")
    assert c.read_bytes() == f.read_bytes()


def test_drained_agglomerates_reach_the_cutoff_in_every_language(fortran_exe, cpp_exe, run_native, tmp_path):
    """Slow agglomerate transport at 3C drains the salt in the agglomerate cores. Every language still runs
    smoothly to the voltage cutoff, quickly, with the same output (docs/model.md section 8)."""
    import time
    g = agg_groups(3.0)
    g["agglomerate"] = dict(g["agglomerate"], tau_agg=2240.0, dt_s=10.0)
    outs = []
    for exe, name in ((fortran_exe, "fd.txt"), (cpp_exe, "cd.txt")):
        t0 = time.time()
        outs.append(run_native(exe, g, name=name))
        assert time.time() - t0 < 30.0
    r = run_corrected(AggParams(mode="corrected", C_rate=3.0, tau_agg=2240.0, **dict(AGG, dt_s=10.0)))
    assert r.exit_reason == "cutoff_low"
    p = tmp_path / "pd.txt"
    r.write(p)
    assert_same(outs[0], p)
    assert_same(outs[1], p)


# ------------------------------------------------------------------ driver: cutoffs and failures
def _last_time_voltage(path):
    """(time [h], voltage [V]) of the last row (two header lines; columns State, Time, Voltage, ...)."""
    tok = path.read_text().splitlines()[-1].split()
    return float(tok[1]), float(tok[2])


def _uniform(which, steps, fortran_exe, cpp_exe, run_native, tmp_path, name):
    if which == "python":
        r = run_uniform(Params(mode="corrected", steps=steps))
        out = tmp_path / f"{name}.txt"
        r.write(out)
        return out, r.exit_reason
    exe = fortran_exe if which == "fortran" else cpp_exe
    if exe is None:
        pytest.skip(f"{which} not built")
    groups = {"numerics": {"mode": "corrected"}, "protocol": {"steps": steps}}
    return run_native(exe, groups, name=f"{which}_{name}.txt"), None


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
def test_discharge_ignores_its_upper_bound(fortran_exe, cpp_exe, run_native, tmp_path, which):
    """A discharge ends at Vmin only: an upper bound below the starting voltage does not stop it
    (the drivers applied both bounds to every cc step and stopped it at once as cutoff_high)."""
    out, reason = _uniform(which, "cc C=1 Vmax=3.5", fortran_exe, cpp_exe, run_native, tmp_path, "bound")
    assert reason in (None, "cutoff_low")
    assert _last_time_voltage(out)[1] == pytest.approx(2.5, abs=1e-4)


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
def test_a_physical_limit_keeps_the_progress_of_its_last_step(fortran_exe, cpp_exe, run_native, tmp_path, which):
    """When a step cannot be completed (here the particles fill at 5C), the exit row is the last converged
    sub-step, inside the time step (the drivers reported the state at its start, a whole number of dt)."""
    out, reason = _uniform(which, "cc C=5 Vmin=0.5", fortran_exe, cpp_exe, run_native, tmp_path, "limit")
    assert reason in (None, "particles_full")
    dt = Params(mode="corrected").dt
    t_h = _last_time_voltage(out)[0]                          # printed in hours to 5 decimals
    grid_h = round(t_h * 3600.0 / dt) * dt / 3600.0
    assert abs(t_h - grid_h) > 2e-5


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
def test_cv_hold_after_a_discharge_proceeds_in_sub_steps(fortran_exe, cpp_exe, run_native, tmp_path, which):
    """A 4.2 V hold right after a 2C discharge cannot be held for whole time steps at first; it proceeds in
    sub-steps and ends on the current limit (without sub-stepping it stopped with particles_full, #3)."""
    out, reason = _uniform(which, "cc C=2 Vmin=2.5; cv V=4.2 Imin=0.05", fortran_exe, cpp_exe, run_native,
                           tmp_path, "cvsub")
    assert reason in (None, "end_of_protocol")
    tok = out.read_text().splitlines()[-1].split()
    assert float(tok[2]) == pytest.approx(4.2, abs=1e-6)                               # voltage held
    assert abs(float(tok[7])) <= 0.05 * Params(mode="corrected").i_1C * 1e3 * (1 + 1e-9)  # current at the limit
