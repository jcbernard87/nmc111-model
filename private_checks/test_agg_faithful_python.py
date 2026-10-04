"""Faithful Python agglomerate model vs the original program (private oracle).

The reference outputs are the original snapshots rebuilt on this machine with
`gfortran -O0 -finit-local-zero` (deviation D-19), in `runs/agg_zero/`. The Python port must
reproduce every line it writes. Full runs take hours in Python, so this checks the first steps
of every snapshot; the compiled ports are checked over full runs.

Requires bandsolver built without FMA contraction (see docs/validation.md).
"""
from pathlib import Path


from nmc_model.agglomerate import run

from .snapshots import agg_params

STEPS = 1500


def _runs(oracle_dir: Path):
    return sorted(p for p in (oracle_dir / "runs" / "agg_zero").iterdir() if p.is_dir())


def pytest_generate_tests(metafunc):
    import os
    if "run_dir" in metafunc.fixturenames:
        root = os.environ.get("NMC_ORACLE_DIR")
        dirs = _runs(Path(root)) if root and Path(root).is_dir() else []
        metafunc.parametrize("run_dir", dirs, ids=[d.name[:40] for d in dirs])


def test_prefix_byte_identical(run_dir, tmp_path):
    ref = (run_dir / "Time_Voltage.txt").read_text().splitlines(True)
    r = run(agg_params(run_dir / "NMC111_agg_Values.f95"), max_steps=STEPS)
    out = tmp_path / "tv.txt"
    r.write(out)
    mine = out.read_text().splitlines(True)
    n = min(len(mine), len(ref))
    assert n > 2
    bad = next((i for i in range(n) if mine[i] != ref[i]), None)
    assert bad is None, f"line {bad + 1} differs:\n ref  {ref[bad]} mine {mine[bad]}"
