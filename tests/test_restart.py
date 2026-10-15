"""Restart files (#13): a run continued from a restart repeats the run it continues."""
import numpy as np
import pytest

from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import run_corrected
from nmc_model.uniform.params import Params
from nmc_model.uniform.simulate import run as run_uniform

S1, S2 = "cc C=2 t=900", "rest t=300; cc C=-1 Vmax=4.2; cv V=4.2 Imin=0.2"
S2_CV = "cv V=3.8 Imin=0.5"                          # a constant-voltage first step uses the stored current
AGG = dict(nj=30, nja=10, dt_s=5.0)


def _run(model, **kw):
    if model == "agglomerate":
        return run_corrected(AggParams(mode="corrected", **AGG, **kw))
    return run_uniform(Params(mode="corrected", particle_model=model, **kw))


@pytest.mark.parametrize("model", ["uniform", "crystal", "agglomerate"])
@pytest.mark.parametrize("s2", [S2, S2_CV], ids=["cc_rest_cv", "cv_first"])
def test_restart_continues_the_run_exactly(tmp_path, model, s2):
    whole = _run(model, steps=f"{S1}; {s2}").array
    f = tmp_path / "r.txt"
    _run(model, steps=S1, restart_write=str(f))
    part = _run(model, steps=s2, restart_read=str(f)).array
    n1 = int(np.sum(whole[:, 7] == 1))                # S1's rows, its end row included
    a, b = whole[n1:], part[1:]                       # a restarted run first repeats the end row
    assert a.shape == b.shape
    cols = [c for c in range(a.shape[1]) if c != 7]   # every column but the step index
    np.testing.assert_array_equal(a[:, cols], b[:, cols])


def test_restart_must_match_the_input(tmp_path):
    f = tmp_path / "r.txt"
    _run("uniform", steps="cc C=1 t=60", restart_write=str(f))
    with pytest.raises(ValueError, match="does not match the input: model"):
        _run("crystal", steps="rest t=60", restart_read=str(f))
    with pytest.raises(ValueError, match="does not match the input: array c"):
        _run("uniform", nj=51, steps="rest t=60", restart_read=str(f))
    with pytest.raises(ValueError, match="mode='corrected'"):
        Params(restart_write=str(f))


def _groups(model, steps, **output):
    g = {"model": {"particle_model": model}, "numerics": {"mode": "corrected"}, "protocol": {"steps": steps}}
    if model == "agglomerate":
        g["cell"] = {"nj": AGG["nj"]}
        g["agglomerate"] = {"nja": AGG["nja"], "dt_s": AGG["dt_s"]}
    if output:
        g["output"] = output
    return g


@pytest.mark.parametrize("model", ["uniform", "crystal", "agglomerate"])
def test_native_restart_continues_the_run_exactly(fortran_exe, cpp_exe, run_native, tmp_path, model):
    """Fortran and C++: S1 then S2 from a restart equals S1; S2 in one run (byte for byte, Li_Nernst aside);
    Python continues from the Fortran program's restart file to the printed precision."""
    for exe, tag in ((fortran_exe, "f"), (cpp_exe, "c")):
        whole = run_native(exe, _groups(model, f"{S1}; {S2}"), name=f"{tag}w.txt").read_text().splitlines()
        rf = tmp_path / f"{tag}r.txt"
        run_native(exe, _groups(model, S1, restart_write=str(rf)), name=f"{tag}1.txt")
        part = run_native(exe, _groups(model, S2, restart_read=str(rf)), name=f"{tag}2.txt").read_text().splitlines()
        rows_w = [r.split() for r in whole[2:]]
        n1 = sum(1 for r in rows_w if r[8] == "1")
        drop = lambda r: r[:8] + r[10:]                 # the step index and Li_Nernst (round-off at rest)
        assert [drop(r) for r in rows_w[n1:]] == [drop(r.split()) for r in part[3:]], tag
    py = _run(model, steps=S2, restart_read=str(tmp_path / "fr.txt")).array
    f2 = np.array([[float(x) for x in r.split()[1:]] for r in (tmp_path / "f2.txt").read_text().splitlines()[2:]])
    np.testing.assert_allclose(py[:, :3], f2[:, :3], rtol=0, atol=1e-5)


@pytest.mark.parametrize("which", ["fortran", "cpp"])
def test_native_restart_must_match_the_input(fortran_exe, cpp_exe, run_native, tmp_path, which):
    import subprocess
    exe = fortran_exe if which == "fortran" else cpp_exe
    f = tmp_path / "r.txt"
    run_native(exe, _groups("uniform", "cc C=1 t=60", restart_write=str(f)), name="a.txt")
    for model, msg in (("crystal", "does not match the input: model"),):
        with pytest.raises(subprocess.CalledProcessError) as e:
            run_native(exe, _groups(model, "rest t=60", restart_read=str(f)), name="b.txt")
        assert msg in e.value.stderr.decode(), e.value.stderr.decode()
