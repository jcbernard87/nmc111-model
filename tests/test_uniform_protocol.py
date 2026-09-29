"""Cycling protocols (docs/protocol.md), uniform model, corrected mode."""
import numpy as np
import pytest

from nmc_model.uniform import kinetics
from nmc_model.uniform.params import Params
from nmc_model.protocol import Step, parse
from nmc_model.uniform.simulate import run

from .uniform_reference_residual import to_physical

CYCLE = "cc C=2 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=600"


def test_parse():
    p = Params(mode="corrected", V_min=2.5, V_max=4.2)
    steps = parse("cc C=1; rest t=60; CV v=4.0 imin=0.05; cc C=-0.5 Vmax=4.1 t=100", p)
    assert steps == [
        Step("cc", C=1.0, Vmin=2.5, Vmax=4.2),
        Step("rest", t=60.0, Vmin=2.5, Vmax=4.2),
        Step("cv", V=4.0, Imin=0.05, Vmin=2.5, Vmax=4.2),
        Step("cc", C=-0.5, t=100.0, Vmin=2.5, Vmax=4.1),
    ]
    assert parse("", p) == [Step("cc", C=p.C_rate, Vmin=2.5, Vmax=4.2)]


@pytest.mark.parametrize("bad", ["cx C=1", "cc", "cc C=1 V=3", "rest", "cv V=4", "cc C=1 t=-1", "cc C"])
def test_parse_rejects(bad):
    with pytest.raises(ValueError):
        parse(bad, Params(mode="corrected"))


@pytest.fixture(scope="module")
def cycle():
    p = Params(mode="corrected", steps=CYCLE)
    return p, run(p)


def test_cycle_completes(cycle):
    p, r = cycle
    assert r.exit_reason == "end_of_protocol"
    a = r.array
    step = a[:, 7]
    ends = [a[step == k][-1] for k in range(1, 6)]
    assert ends[0][1] == pytest.approx(2.5, abs=1e-4)          # discharge stops at Vmin
    assert ends[1][6] == 0.0                                   # rest carries no current
    assert ends[2][1] == pytest.approx(4.0, abs=1e-4)          # charge stops at Vmax
    assert abs(ends[3][6]) <= 0.05 * p.i_1C * 1e3              # CV ends on the current limit
    assert ends[3][1] == pytest.approx(4.0, abs=1e-8)


def test_cycle_conserves_lithium(cycle):
    """The lithium left in the particles after the cycle equals the net charge passed."""
    from nmc_model.uniform.model import make_mesh
    p, r = cycle
    m = make_mesh(p)
    s = m.s
    cs = to_physical(p, r.final_state)[:, 3]
    solid = (p.vf_AM * (cs[s + 1:-1] - p.cs_init) * m.dx[s + 1:-1]).sum()   # mol/cm2
    net = r.array[-1, 2] * p.mass_area / p.M               # electron equivalents -> mol/cm2
    assert solid == pytest.approx(net, rel=1e-6)


def test_rest_relaxes_to_ocp(cycle):
    p, r = cycle
    c = to_physical(p, r.final_state)
    u = kinetics.ocp(p, c[p.sep_node:-1, 0], c[p.sep_node:-1, 3])
    np.testing.assert_allclose(c[p.sep_node:-1, 1] - c[p.sep_node:-1, 2], u, atol=2e-3)
