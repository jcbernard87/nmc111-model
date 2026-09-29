import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _exe(env_name, default):
    p = Path(os.environ.get(env_name, ROOT / default))
    if not p.is_absolute():
        p = (Path.cwd() / p).resolve()
    return p if p.is_file() and os.access(p, os.X_OK) else None


@pytest.fixture
def fortran_exe():
    exe = _exe("NMC_FORTRAN_EXE", "build/fortran/nmc_f")
    if exe is None:
        pytest.skip("Fortran program not built (cmake -S . -B build && cmake --build build)")
    return exe


@pytest.fixture
def cpp_exe():
    exe = _exe("NMC_CPP_EXE", "build/cpp/nmc_cpp")
    if exe is None:
        pytest.skip("C++ program not built (cmake -S . -B build && cmake --build build)")
    return exe


@pytest.fixture
def run_native(tmp_path):
    """Run a compiled program on a namelist given as {group: {name: value}}; return the output path."""

    def _fmt(v):
        if isinstance(v, str):
            return f"'{v}'"
        if isinstance(v, bool):
            return ".true." if v else ".false."
        return repr(v)

    def _run(exe, groups, name="out.txt"):
        groups = {k: dict(v) for k, v in groups.items()}
        groups.setdefault("output", {})["file"] = name
        text = "".join(f"&{g}\n" + "".join(f"  {k} = {_fmt(v)}\n" for k, v in kv.items()) + "/\n"
                       for g, kv in groups.items())
        nml = tmp_path / f"in_{name}.nml"
        nml.write_text(text)
        subprocess.run([str(exe), str(nml)], cwd=tmp_path, check=True, capture_output=True)
        return tmp_path / name

    return _run
