"""Corrected agglomerate model (docs/model.md sections 6 and 8): conservation, the condensed Newton step,
limits (exhausted electrolyte, full or empty crystals), convergence."""
import math
from dataclasses import replace

import numpy as np
import pytest

from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import AggState, CorrectedModel, run_corrected

U, P1, P2, S = 0, 1, 2, 3          # the state holds u = ln(c/c_bulk) and s = logit(theta)


def params(**kw):
    base = dict(mode="corrected", nj=30, nja=10, dt_s=5.0)
    base.update(kw)
    return AggParams(**base)


def salt(m, p, st):
    """Salt in the macro-pores and in the agglomerate pores [mol/cm2]."""
    me = m.mesh
    s = me.s
    ce, ca = m.conc(st)
    macro = (p.eps_sep * ce[1:s] * me.dx[1:s]).sum() + (p.eps * ce[s + 1:-1] * me.dx[s + 1:-1]).sum()
    v_one = 4.0 / 3.0 * math.pi * p.R_agg ** 3
    micro = ((p.eps_agg * ca * m.dV).sum(axis=1) / v_one * m.v_agg * me.dx[m.nodes]).sum()
    return macro + micro


def lithium_in_crystals(m, p, st):
    v_one = 4.0 / 3.0 * math.pi * p.R_agg ** 3
    per_agg = ((m.cs(st) - p.cs_init) * m.dV).sum(axis=1) * (1.0 - p.eps_agg) / v_one * m.v_agg
    return (per_agg * m.mesh.dx[m.nodes]).sum()


@pytest.fixture(scope="module")
def discharge():
    p = params()
    return p, run_corrected(p)


def test_discharge_ends_at_cutoff(discharge):
    p, r = discharge
    assert r.exit_reason == "cutoff_low"
    v = r.array[:, 1]
    assert v[-1] <= p.V_min and v[-2] > p.V_min
    assert abs(v[-1] - p.V_min) < 1e-4


def test_lithium_conserved(discharge):
    """Lithium in the crystals equals the charge passed (1C throughout), to round-off."""
    p, r = discharge
    m = CorrectedModel(p)
    passed = r.array[-1, 0] * 3600.0 * p.i_1C / p.F
    assert lithium_in_crystals(m, p, r.final_state) / passed == pytest.approx(1.0, abs=1e-12)


def test_salt_conserved(discharge):
    """Salt in the macro-pores plus the agglomerate pores is constant (Li+ in at the foil = Li+ into the crystals)."""
    p, r = discharge
    m = CorrectedModel(p)
    s0 = salt(m, p, m.initial_state())
    assert abs(salt(m, p, r.final_state) - s0) / s0 < 1e-13


def test_condensed_step_solves_the_full_system():
    """After a condensed Newton step, both scales' uncondensed residuals vanish.

    Each residual column is compared with its size at a randomly perturbed state.
    """
    p = params()
    m = CorrectedModel(p)
    old = m.initial_state()
    st = m.newton_step(old, p.dt, p.i_1C)
    rng = np.random.default_rng(7)
    x = AggState(st.c + 1e-3 * rng.standard_normal(st.c.shape), st.ca + 1e-3 * rng.standard_normal(st.ca.shape))
    Re, Ra = m.residuals(st, old, p.dt, p.i_1C)
    Rex, Rax = m.residuals(x, old, p.dt, p.i_1C)
    for R, Rx, cols in ((Re, Rex, (U, P1, P2)), (Ra, Rax, (U, P1, P2, S))):
        for col in cols:
            assert np.abs(R[..., col]).max() <= 1e-9 * np.abs(Rx[..., col]).max()


def test_newton_converges_quadratically():
    """From a perturbed iterate, the scaled update falls quadratically: the condensed Jacobian is exact."""
    p = params(newton_tol=1e-14)
    m = CorrectedModel(p)
    old = m.initial_state()
    sol = m.newton_step(old, p.dt, p.i_1C)
    rng = np.random.default_rng(3)
    x = AggState(sol.c.copy(), sol.ca.copy())
    x.c[:, P1:P2 + 1] += 1e-3 * rng.standard_normal((p.nj, 2))
    x.ca[..., P1:P2 + 1] += 1e-3 * rng.standard_normal(x.ca[..., P1:P2 + 1].shape)
    x.ca[..., U] += 1e-2 * rng.standard_normal(x.ca[..., U].shape)
    hist = []
    m.newton_step(old, p.dt, p.i_1C, start=x, history=hist)
    assert len(hist) <= 6
    for a, b in zip(hist, hist[1:]):
        if a > 1e-7:
            assert b < 50.0 * a ** 2


def test_charge_neutral_agglomerates(discharge):
    """Each agglomerate takes in as much electronic as it gives ionic current: i1,in + i2,in = 0."""
    p, r = discharge
    m = CorrectedModel(p)
    q, _ = m.surface_flux(r.final_state.ca)
    assert np.abs(q[:, 1] + q[:, 2]).max() <= 1e-7 * np.abs(q[:, 1]).max()   # Newton tolerance level


def test_rest_at_equilibrium_stays_put():
    p = params()
    m = CorrectedModel(p)
    st = m.initial_state()
    u = float(m.kin.ocp(0.0, st.ca[0, 0, S]))
    st.c[:, P1] = u
    st.ca[..., P1] = u
    new = m.newton_step(st, 10.0, 0.0)
    assert np.abs(new.c - st.c).max() < 1e-12 and np.abs(new.ca - st.ca).max() < 1e-12


def test_fast_agglomerate_transport_makes_agglomerates_uniform():
    """With very fast pore transport and conduction inside the agglomerates, their crystals lithiate evenly."""
    slow = params(dt_s=10.0)
    fast = replace(slow, tau_agg=1e-4, sigma_agg=1e3)
    spread = []
    for p in (slow, fast):
        st = run_corrected(p, max_steps=120).final_state
        cs = CorrectedModel(p).cs(st)[:, 1:-1]
        spread.append(float(((cs.max(axis=1) - cs.min(axis=1)) / cs.mean(axis=1)).max()))
    assert spread[1] < 1e-3 and spread[1] < spread[0]


def test_protocol_cycle_returns_the_lithium():
    p = params(steps="cc C=1 Vmin=3.0; rest t=600; cc C=-1 Vmax=4.3; cv V=4.3 Imin=0.05; rest t=600")
    r = run_corrected(p)
    assert r.exit_reason == "end_of_protocol"
    m = CorrectedModel(p)
    li = lithium_in_crystals(m, p, r.final_state)
    net = r.array[-1, 2] * p.mass_area / p.M        # electron equivalents x moles of active material
    assert li == pytest.approx(net, rel=1e-9, abs=1e-15)


def _v600(**kw):
    """Cell voltage after 600 s at 2C."""
    base = dict(nj=30, nja=10, dt_s=10.0, C_rate=2.0, steps="cc C=2 t=600")
    base.update(kw)
    return run_corrected(params(**base)).array[-1, 1]


def _orders(values):
    e = np.abs(np.diff(values))
    return [math.log2(e[i] / e[i + 1]) for i in range(len(e) - 1)]


def test_time_convergence_first_order():
    """Backward Euler: halving the step halves the error."""
    orders = _orders([_v600(dt_s=d) for d in (20.0, 10.0, 5.0, 2.5)])
    assert all(0.85 < o < 1.15 for o in orders), orders


def test_radial_mesh_convergence_second_order():
    """Halving the agglomerate's radial spacing (4, 8, 16, 32 volumes) quarters the error."""
    orders = _orders([_v600(nja=n) for n in (6, 10, 18, 34)])
    assert all(1.9 < o < 2.1 for o in orders), orders


def test_electrode_mesh_convergence_second_order():
    """Halving the cathode spacing (8, 16, 32, 64 volumes) quarters the error."""
    orders = _orders([_v600(nj=22 + 1 + n) for n in (8, 16, 32, 64)])
    assert all(1.9 < o < 2.1 for o in orders), orders


# ------------------------------------------------------------------ physical limits (docs/model.md section 8)

@pytest.mark.parametrize("C_rate,tau_agg", [(1.0, 22400.0), (3.0, 2240.0)])
def test_exhausted_agglomerate_cores_run_to_the_cutoff(C_rate, tau_agg):
    """With slow transport inside the agglomerates their cores run out of salt; the discharge still ends
    smoothly at the voltage cutoff, with lithium and salt conserved."""
    p = params(dt_s=10.0, C_rate=C_rate, tau_agg=tau_agg)
    r = run_corrected(p)
    assert r.exit_reason == "cutoff_low"
    m = CorrectedModel(p)
    st = r.final_state
    _, ca = m.conc(st)
    assert ca.min() < 1e-6 * p.c_bulk                     # the cores are drained
    passed = r.array[-1, 0] * 3600.0 * C_rate * p.i_1C / p.F
    assert lithium_in_crystals(m, p, st) / passed == pytest.approx(1.0, abs=1e-12)
    s0 = salt(m, p, m.initial_state())
    assert abs(salt(m, p, st) - s0) / s0 < 1e-12


def test_charge_empties_the_crystals_smoothly():
    p = params(dt_s=10.0, steps="cc C=1 Vmin=3.0; rest t=600; cc C=-1 Vmax=4.4; cv V=4.4 Imin=0.02; rest t=600")
    r = run_corrected(p)
    assert r.exit_reason == "end_of_protocol"
    th = CorrectedModel(p).cs(r.final_state) / CorrectedModel(p).cs_max
    assert th.max() < 1e-3                                  # nearly empty, reached without clipping


def test_sigmoid_keeps_precision_at_the_limits():
    """theta = sigmoid(s) and 1 - theta = sigmoid(-s) keep full relative precision far into the tails
    (a 0.5 (1 + tanh(s/2)) form rounds theta to exactly 0 below s ~ -37)."""
    from nmc_model.logcore import sigmoid
    for s in (-40.0, -200.0, -700.0):
        assert sigmoid(s) == pytest.approx(math.exp(s), rel=1e-12)
        assert sigmoid(-s) == 1.0
    assert sigmoid(0.0) == 0.5
