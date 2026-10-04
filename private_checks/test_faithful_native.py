"""Compiled programs in faithful mode vs the original programs (private oracle).

The agglomerate references are the archived snapshots rebuilt with gfortran -O0
-finit-local-zero (runs/agg_zero); the uniform-particle reference is the archived output.
Full agglomerate runs take a few minutes each, so set NMC_FULL_RUNS=1 to run all of them;
otherwise only the shortest snapshot is checked.
"""
import os
import subprocess
from pathlib import Path

import pytest

from .snapshots import agg_namelist, uniform_namelist

ROOT = Path(__file__).resolve().parents[1]


def _exe(env_name, default):
    p = Path(os.environ.get(env_name, ROOT / default))
    if not p.is_absolute():
        p = (Path.cwd() / p).resolve()
    return p if p.is_file() and os.access(p, os.X_OK) else None


def _agg_dirs():
    root = os.environ.get("NMC_ORACLE_DIR")
    if not root or not Path(root).is_dir():
        return []
    dirs = sorted(p for p in (Path(root) / "runs" / "agg_zero").iterdir() if p.is_dir())
    if os.environ.get("NMC_FULL_RUNS") != "1":
        dirs = [min(dirs, key=lambda d: (d / "Time_Voltage.txt").stat().st_size)]
    return dirs


EXES = [("NMC_FORTRAN_EXE", "build/fortran/nmc_f"), ("NMC_CPP_EXE", "build/cpp/nmc_cpp")]


@pytest.mark.parametrize("exe_var,default", EXES)
@pytest.mark.parametrize("run_dir", _agg_dirs(), ids=lambda d: d.name[:40])
def test_agglomerate_byte_identical(run_dir, exe_var, default, tmp_path):
    exe = _exe(exe_var, default)
    if exe is None:
        pytest.skip(f"{default} not built")
    (tmp_path / "in.nml").write_text(agg_namelist(run_dir / "NMC111_agg_Values.f95", "tv.txt"))
    subprocess.run([str(exe), "in.nml"], cwd=tmp_path, check=True, capture_output=True)
    assert (tmp_path / "tv.txt").read_bytes() == (run_dir / "Time_Voltage.txt").read_bytes()


@pytest.mark.parametrize("exe_var,default", EXES)
def test_uniform_byte_identical(oracle_dir, exe_var, default, tmp_path):
    exe = _exe(exe_var, default)
    if exe is None:
        pytest.skip(f"{default} not built")
    (tmp_path / "in.nml").write_text(uniform_namelist(oracle_dir, "tv.txt"))
    subprocess.run([str(exe), "in.nml"], cwd=tmp_path, check=True, capture_output=True)
    ref = oracle_dir / "archive_outputs" / "crystal" / "blackbox_Time_Voltage.txt"
    assert (tmp_path / "tv.txt").read_bytes() == ref.read_bytes()
