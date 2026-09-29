"""Faithful Python uniform-particle model vs the original program's archived output (private oracle)."""
from nmc_model.uniform.params import Params
from nmc_model.uniform.simulate import run

from .snapshots import uniform_fit_values


def test_byte_identical_time_voltage(oracle_dir, tmp_path):
    r = run(Params.faithful(C_rate=1.0, **uniform_fit_values(oracle_dir)))
    assert r.exit_reason == "nan"          # deviation D-3: the original ends on a NaN
    out = tmp_path / "Time_Voltage.txt"
    r.write(out)
    ref = oracle_dir / "archive_outputs" / "crystal" / "blackbox_Time_Voltage.txt"
    assert out.read_bytes() == ref.read_bytes()
