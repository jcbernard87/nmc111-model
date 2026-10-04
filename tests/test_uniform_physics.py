"""Physics and numerics checks. They need no stored data: everything is computed here.

Corrected mode must pass every check. Faithful mode is checked to still show the documented
defects (docs/deviations.md), so that a change which silently alters faithful behaviour fails.
"""
import numpy as np
import pytest

import bandsolver
from nmc_model.uniform import kinetics
from nmc_model.uniform.model import Assembler, make_mesh
from nmc_model.uniform.params import Params
from nmc_model.uniform.logmodel import LogUniformModel
from nmc_model.uniform.simulate import initial_state, run

from .uniform_reference_residual import residual, sg_residual, to_physical


# Faithful mode needs the original's fitted k and sigma, which are not distributed. Its defects are
# checked here with round values; k = 1e-7 keeps its one linearized solve per step stable (the
# generic corrected-mode k diverges at once in faithful mode, which has no Newton iteration, D-7).
FAITHFUL_K = 1.0e-7


def corrected(**kw):
    return Params(mode="corrected", **kw)


def inventories(p, c):
    """Salt in the electrolyte and lithium in the solid, per unit area [mol/cm2] (a state in c-form)."""
    m = make_mesh(p)
    s, dx = m.s, m.dx
    salt = (p.eps_sep * c[1:s, 0] * dx[1:s]).sum() + (p.eps * c[s + 1:-1, 0] * dx[s + 1:-1]).sum()
    solid = (p.vf_AM * c[s + 1:-1, 3] * dx[s + 1:-1]).sum()
    return salt, solid


@pytest.fixture(scope="module")
def mid_discharge():
    """Corrected-mode state (log variables) after 600 s at 1C, and the previous state."""
    p = corrected(C_rate=1.0)
    x_prev = run(p, max_steps=599).final_state
    x = run(p, max_steps=600).final_state
    return p, x_prev, x


# ------------------------------------------------------------------ residual and Jacobian

def test_residual_matches_independent_reference(mid_discharge):
    """The assembled corrected-mode residual equals an independent node-by-node implementation."""
    p, x_old, x = mid_discharge
    rng = np.random.default_rng(1)
    xp = x + 1e-4 * rng.standard_normal(x.shape)            # away from the solution
    R = LogUniformModel(p).el.residual_and_blocks(xp, x_old, 1.0, p.i_app)[0]
    R_ref = sg_residual(p, xp, x_old, 1.0, p.i_app)
    rel = np.abs(R - R_ref).max(axis=0) / np.abs(R_ref).max(axis=0)
    assert rel.max() < 1e-10, rel


def test_faithful_residual_shows_D2_and_D11(mid_discharge):
    """Faithful mode differs from the intended equations in the cation row (D-2) and current row (D-11)."""
    pc, x_old, x = mid_discharge
    c_old, c = to_physical(pc, x_old), to_physical(pc, x)
    p = Params.faithful(k_rxn=FAITHFUL_K, sigma=0.1, C_rate=1.0)
    asm = Assembler(p)
    _, _, _, G, _ = asm.assemble(c, 1.0)
    G = G - asm.time_terms(1.0) * (c - c_old)
    R = residual(p.with_(mode="corrected"), c, c_old, 1.0)
    err = np.abs(G + R).max(axis=0) / np.abs(R).max(axis=0)
    assert err[0] > 0.1 and err[2] > 0.1           # D-2 (cation flux) and D-11 (current)


def _fill(asm, c_old, dt=1.0):
    def fill(c):
        A, B, D, G, _ = asm.assemble(c, dt)
        return A, B, D, G - asm.time_terms(dt) * (c - c_old)
    return fill


def test_jacobian_corrected(mid_discharge):
    p, x_old, x = mid_discharge
    el = LogUniformModel(p).el

    def fill(y):
        R, A, B, D = el.residual_and_blocks(y, x_old, 1.0, p.i_app)
        return A, B, D, -R
    chk = bandsolver.check_jacobian(fill, x)
    assert chk.max_error < 1e-4, chk.worst()


def test_jacobian_faithful_shows_D1(mid_discharge):
    """The Li-foil solid-potential row has the wrong sign in faithful mode (D-1)."""
    pc, x_old, x = mid_discharge
    c_old, c = to_physical(pc, x_old), to_physical(pc, x)
    asm = Assembler(Params.faithful(k_rxn=FAITHFUL_K, sigma=0.1, C_rate=1.0))
    chk = bandsolver.check_jacobian(_fill(asm, c_old), c)
    name, worst = chk.worst()
    assert (name, worst.node, worst.row, worst.col) == ("B", 0, 1, 1)
    assert worst.user == pytest.approx(-worst.fd, rel=1e-6)


# ------------------------------------------------------------------ conservation

def test_conservation_corrected():
    p = corrected(C_rate=1.0)
    m = LogUniformModel(p)
    x0 = m.initial_state()
    x = run(p, max_steps=1800).final_state
    s0, so0 = inventories(p, to_physical(p, x0))
    s1, so1 = inventories(p, to_physical(p, x))
    q = p.i_app * 1800 / p.F
    assert abs(s1 - s0) / s0 < 1e-12                  # salt is conserved
    assert abs((so1 - so0) / q - 1) < 1e-12           # lithium into the solid equals It/F
    el = m.el
    Fv, _, _ = el.tr.fluxes(x[:-1], x[1:], el.g, el.gs)
    sep = slice(1, m.mesh.s - 1)
    np.testing.assert_allclose(Fv[sep, 2], p.i_app, rtol=1e-9)        # the separator carries the applied current
    N_minus = (Fv[sep, 0] * p.F - Fv[sep, 2]) / p.F                  # i2 = F (N+ - N-), background current negligible
    assert np.abs(N_minus).max() < 1e-6 * p.i_app / p.F              # and no anion flux (quasi-steady)


def test_conservation_faithful_shows_D2():
    p = Params.faithful(k_rxn=FAITHFUL_K, sigma=0.1, C_rate=1.0)
    c0 = initial_state(p)
    c = run(p, max_steps=1800).final_state
    s0, _ = inventories(p, c0)
    s1, _ = inventories(p, c)
    assert -4e-3 < (s1 - s0) / s0 < -3e-3              # 0.34 % of the salt is lost at start-up


# ------------------------------------------------------------------ equilibrium

def test_rest_at_equilibrium_stays_put():
    """No current, open-circuit potential everywhere: one step must change nothing."""
    p = corrected(C_rate=0.0)
    m = LogUniformModel(p)
    x0 = m.initial_state()
    x0[:, 1] = float(m.kin.ocp(0.0, x0[m.mesh.s, 3]))
    x = m.newton_step(x0, 10.0, 0.0)
    assert np.abs(x - x0).max() < 1e-12


# ------------------------------------------------------------------ convergence

def test_time_convergence_first_order():
    """Backward Euler: the Li-face concentration at t = 8 s converges at first order in dt."""
    vals = []
    for dt in (1.0, 0.5, 0.25, 0.125):
        p = corrected(C_rate=2.0, n_steps=int(36000 / dt), newton_tol=1e-13)
        vals.append(run(p, max_steps=int(8 / dt)).final_state[0, 0])
    d = np.diff(vals)
    orders = np.log2(np.abs(d[:-1] / d[1:]))
    assert np.all((orders > 0.9) & (orders < 1.1)), orders


def test_mesh_convergence_second_order():
    """The finite-volume discretization converges at second order in the mesh size."""
    vals = []
    for sep, nj in ((7, 26), (12, 51), (22, 101), (42, 201)):
        p = corrected(C_rate=2.0, sep_node=sep, nj=nj)
        vals.append(run(p, max_steps=900).final_state[-1, 1])
    d = np.diff(vals)
    orders = np.log2(np.abs(d[:-1] / d[1:]))
    assert np.all((orders > 1.8) & (orders < 2.6)) and orders[-1] < orders[0], orders   # approaching 2


# ------------------------------------------------------------------ end of discharge

@pytest.mark.parametrize("C_rate", [2.0, 1.0])
def test_corrected_ends_at_cutoff(C_rate):
    r = run(corrected(C_rate=C_rate))
    assert r.exit_reason == "cutoff_low"
    v = r.array[:, 1]
    assert v[-1] <= 2.5 and v[-2] > 2.5
    assert 0.54 < r.array[-1, 2] < 0.55       # electron equivalents: nearly all of x_max = 0.55


# ------------------------------------------------------------------ reference electrode and OCP

def test_potentials_are_gauge_invariant(mid_discharge):
    """The residual depends only on potential differences, so the foil reference is an exact shift."""
    from nmc_model.uniform.simulate import foil_referenced
    p, x_prev, x = mid_discharge
    rng = np.random.default_rng(1)
    y = x + 1.0e-3 * rng.standard_normal(x.shape)            # off the solution, so the residual is not round-off
    r0 = sg_residual(p, y, x_prev, p.dt, p.i_app)
    r1 = sg_residual(p, foil_referenced(p, y, p.i_app), x_prev, p.dt, p.i_app)
    r0[0, 2] = r1[0, 2] = 0.0                                 # the gauge row phi2(foil face) = 0 itself
    np.testing.assert_allclose(r1, r0, rtol=0, atol=1e-9 * np.abs(r0).max())


def test_cell_voltage_against_foil(mid_discharge):
    """V = phi1(collector) - (U_Li + eta_Li): phi2 at the foil face is -(U_Li + eta_Li) on the foil scale."""
    from nmc_model.uniform.simulate import cell_voltage, foil_referenced
    p, _, x = mid_discharge
    c0 = p.c_bulk * np.exp(x[0, 0])
    u_li, eta_li = kinetics.li_foil(p, c0, p.i_app)
    ref = foil_referenced(p, x, p.i_app)
    assert ref[0, 2] == pytest.approx(-(u_li + eta_li), abs=1e-15)
    assert cell_voltage(p, x, p.i_app) == pytest.approx(ref[-1, 1], abs=1e-15)
    assert eta_li > 0 and c0 > p.c_bulk and u_li > 0          # discharge: salt builds up at the foil


@pytest.mark.parametrize("model", ["uniform", "agglomerate"])
def test_ocp_electrolyte_term(model):
    """The OCP carries (RT/F) ln(c/c_bulk): doubling c raises U by (RT/F) ln 2 to 1e-12, in the corrected
    kinetics of both models (log variables) and in the c-form OCP that faithful mode uses (the original has the
    term too), and the two forms agree at the same state (#6). theta = 1/2 is included: the Redlich-Kister sum's
    k = 0 term used to evaluate 0/0 there (NaN) in corrected mode."""
    import math
    from nmc_model.uniform import kinetics
    from nmc_model.agglomerate import AggParams
    from nmc_model.agglomerate_corrected import CorrectedModel
    if model == "uniform":
        p = Params(mode="corrected")
        kin = LogUniformModel(p).kin
    else:
        p = AggParams(mode="corrected")
        kin = CorrectedModel(p).kin
    rtf = p.R * p.T / p.F
    th = np.array([0.3, 0.5, 0.7])
    s = np.log(th / (1.0 - th))
    shift = kin.ocp(np.full(3, math.log(2.0)), s) - kin.ocp(np.zeros(3), s)
    np.testing.assert_allclose(shift, rtf * math.log(2.0), rtol=1e-12, atol=0.0)
    if model == "uniform":
        cs = th * kinetics.cs_max(p)
        for q in (p, Params.faithful(k_rxn=1e-6, sigma=0.1)):
            # faithful mode keeps the original's 0/0 at theta = 1/2 exactly (D-20): compare away from it
            d = kinetics.ocp(q, 2.0 * q.c_bulk, cs) - kinetics.ocp(q, q.c_bulk, cs)
            if q.mode == "faithful":
                d = d[[0, 2]]
            np.testing.assert_allclose(d, q.R * q.T / q.F * math.log(2.0), rtol=1e-6, atol=0.0)
        np.testing.assert_allclose(kin.ocp(np.zeros(3), s), kinetics.ocp(p, p.c_bulk, cs), rtol=1e-12)
