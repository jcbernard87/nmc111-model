"""The crystal model (particle_model = 'crystal', #4): solid diffusion in the uniform model's particles, ported
from lfp-model with NMC111's kinetics. The same checks as lfp-model's tests/test_crystal.py."""
import numpy as np
import pytest

import bandsolver
from nmc_model.logcore import logit, sigmoid
from nmc_model.uniform.crystal import CrystalMesh, CrystalModel, CrystalState, diffusion_rows
from nmc_model.uniform.params import Params
from nmc_model.uniform.simulate import run


def _run_flux(k, nc, J=0.05, th0=0.2, dt=0.05, t_end=2.0):
    """Unit crystal (R = D = cs_max = 1), constant inward flux J; Newton in s each step."""
    cm = CrystalMesh(1.0, nc, k)
    xc = np.full((1, nc), logit(th0))
    for _ in range(round(t_end / dt)):
        old = xc.copy()
        for _ in range(30):
            R, B, Dn, An = diffusion_rows(cm, 1.0, 1.0, xc, old, dt)
            R[:, -1] -= cm.A_R * J               # the caller adds the outward surface flux A_R i/F; inward J: i/F = -J
            A3 = An[0][:, None, None]; B3 = B[0][:, None, None]; D3 = Dn[0][:, None, None]
            dx = bandsolver.solve(A3, B3, D3, -R[0][:, None])[:, 0]
            xc = xc + dx[None, :]
            if np.abs(dx).max() < 1e-13:
                break
    th = sigmoid(xc[0])
    exact = th0 + (k + 1) * J * t_end + J / 2 * (cm.r ** 2 - (k + 1) / (k + 3))
    return cm, th, exact


@pytest.mark.parametrize("k", [0, 1, 2])
def test_constant_flux_matches_analytic(k):
    cm, th, exact = _run_flux(k, 41)
    assert np.abs(th - exact).max() < 2e-4
    mean = (cm.V * th).sum() / cm.Vtot
    assert mean == pytest.approx(0.2 + (k + 1) * 0.05 * 2.0, abs=1e-12)    # lithium balance, exact


@pytest.mark.parametrize("k", [0, 2])
def test_crystal_mesh_second_order(k):
    errs = []
    for nc in (11, 21, 41):
        _, th, exact = _run_flux(k, nc)
        errs.append(abs(th[-1] - exact[-1]))    # surface value
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))   # h = R/10, R/20, R/40
    assert np.all(orders > 1.8), orders


def test_diffusion_rows_jacobian():
    rng = np.random.default_rng(3)
    cm = CrystalMesh(2e-5, 9, 2)
    xc = rng.normal(0.0, 2.0, (2, 9)); old = xc + rng.normal(0.0, 0.1, xc.shape)
    R, B, Dn, An = diffusion_rows(cm, 5e-12, 0.0478, xc, old, 1.0)
    h = 1e-6
    for j in range(9):
        e = np.zeros_like(xc); e[:, j] = h
        dR = (diffusion_rows(cm, 5e-12, 0.0478, xc + e, old, 1.0)[0]
              - diffusion_rows(cm, 5e-12, 0.0478, xc - e, old, 1.0)[0]) / (2 * h)
        want = np.zeros_like(xc)
        want[:, j] = B[:, j]
        if j > 0:
            want[:, j - 1] = Dn[:, j - 1]
        if j < 8:
            want[:, j + 1] = An[:, j + 1]
        np.testing.assert_allclose(dR, want, rtol=1e-6, atol=1e-6 * np.abs(want).max())


def crystal(**kw):
    return Params(mode="corrected", particle_model="crystal", **kw)


def test_crystal_area_by_shape():
    for shape, k in (("slab", 0), ("cylinder", 1), ("sphere", 2)):
        p = crystal(crystal_shape=shape)
        assert CrystalModel(p).a == pytest.approx((k + 1) * p.vf_AM / p.R_p, rel=1e-15)


def test_condensed_step_solves_full_system():
    """At the condensed step's solution, a Newton update of the full (uncondensed) system is round-off,
    and the condensed iteration converges quadratically."""
    p = crystal(C_rate=1.0)
    m = CrystalModel(p)
    st = m.initial_state()
    for _ in range(20):
        st = m.newton_step(st, 1.0, p.i_app)
    hist = []
    new = m.newton_step(st, 1.0, p.i_app, history=hist)
    assert len(hist) <= 5 and hist[-1] < 1e-12
    J = m.jacobian_dense(new, st, 1.0, p.i_app)
    Re, Rc = m.residuals(new, st, 1.0, p.i_app)
    d = np.linalg.solve(J, -np.concatenate([Re.ravel(), Rc.ravel()]))
    assert np.abs(d).max() < 1e-12


def test_crystal_conservation():
    p = crystal(C_rate=1.0)
    m = CrystalModel(p)
    st = m.initial_state()
    x0 = st.copy()
    for _ in range(300):
        st = m.newton_step(st, 1.0, p.i_app)
    dx = m.mesh.dx[m.nodes]
    li = lambda s_: (m.cs_mean(s_) * p.vf_AM * dx).sum()
    q = p.i_app * 300 / p.F
    assert abs((li(st) - li(x0)) / q - 1.0) < 1e-11           # a discharge (I > 0) lithiates: + It/F
    salt = lambda s_: (m.el.eps_node * m.mesh.dx * m.conc(s_)).sum()
    assert abs(salt(st) / salt(x0) - 1.0) < 1e-12


def test_condensed_jacobian_matches_fd():
    """Full (uncondensed) Jacobian of (Re, Rc) vs central differences at a mid-discharge state."""
    p = crystal(C_rate=1.0, nj=26, sep_node=7, nj_crystal=6)
    m = CrystalModel(p)
    st = m.initial_state()
    for _ in range(200):
        st = m.newton_step(st, 1.0, p.i_app)
    old = st.copy()
    rng = np.random.default_rng(2)
    st = CrystalState(st.x + 1e-5 * rng.standard_normal(st.x.shape), st.xc + 1e-3 * rng.standard_normal(st.xc.shape))
    J = m.jacobian_dense(st, old, 1.0, p.i_app)
    flat = lambda s_: np.concatenate([a.ravel() for a in m.residuals(s_, old, 1.0, p.i_app)])
    z0 = np.concatenate([st.x.ravel(), st.xc.ravel()])
    ne = st.x.size
    Jfd = np.zeros_like(J)
    for q in range(z0.size):
        h = 1e-7 * max(1.0, abs(z0[q]))
        zp, zm = z0.copy(), z0.copy(); zp[q] += h; zm[q] -= h
        unpack = lambda z: CrystalState(z[:ne].reshape(st.x.shape), z[ne:].reshape(st.xc.shape))
        Jfd[:, q] = (flat(unpack(zp)) - flat(unpack(zm))) / (2 * h)
    scale = np.abs(Jfd).max(axis=1, keepdims=True) + 1e-300
    assert (np.abs(J - Jfd) / scale).max() < 1e-5


def test_uniform_limit():
    """Fast solid diffusion: the crystal model's voltage approaches the uniform-particle model's, with an error
    that falls exactly as 1/D_c (max|dV| D_c = 1.16e-14 V cm2/s measured at D_c = 1e-10, 1e-9 and 1e-8)."""
    a_u = run(Params(mode="corrected", C_rate=1.0), max_steps=1800).array
    err = {}
    for D_c in (1e-9, 1e-6):
        r_c = run(crystal(C_rate=1.0, D_c=D_c), max_steps=1800)
        assert isinstance(r_c.final_state, CrystalState)
        a_c = r_c.array
        n = min(len(a_u), len(a_c))
        err[D_c] = np.abs(a_c[:n, 1] - a_u[:n, 1]).max()
    assert err[1e-6] < 1e-7, err
    assert err[1e-9] / err[1e-6] == pytest.approx(1000.0, rel=0.01), err


CYCLE = "cc C=2 Vmin=3.0; rest t=600; cc C=-1 Vmax=4.2; cv V=4.2 Imin=0.05; rest t=600"


def test_crystal_discharge_ends_at_cutoff():
    r = run(crystal(C_rate=1.0))
    assert isinstance(r.final_state, CrystalState)            # the crystal model actually ran
    assert r.exit_reason == "cutoff_low"
    assert r.array[-1, 1] <= 2.5 and r.array[-2, 1] > 2.5


def test_capacity_falls_with_slower_diffusion():
    caps = [run(crystal(C_rate=2.0, D_c=D)).array[-1, 2] for D in (1e-10, 5e-12, 5e-13)]
    assert caps[0] > caps[1] > caps[2]


def test_tiny_diffusivity_still_reaches_cutoff():
    r = run(crystal(C_rate=1.0, D_c=1e-15))
    assert isinstance(r.final_state, CrystalState)
    assert r.exit_reason == "cutoff_low"


def test_crystal_cycle_completes():
    """The CV hold drains the crystal surfaces to theta ~ 1e-6, where a linearized log-odds update is
    huge; Newton must still converge (divergence is judged on the damped electrode update)."""
    r = run(crystal(steps=CYCLE))
    assert isinstance(r.final_state, CrystalState)
    assert r.exit_reason == "end_of_protocol"


def _theta_spread(m, st):
    """Largest spread of theta inside any one crystal."""
    th = m.cs(st) / m.cs_max
    return float(np.max(th.max(axis=1) - th.min(axis=1)))


def test_rest_relaxes_monotonically():
    """After a 2C partial discharge the crystal profiles flatten at I = 0; V moves monotonically
    toward the OCP of the mean theta, and the spread of theta inside the crystals falls during the rest and
    ends below 1e-3."""
    p = crystal(steps="cc C=2 t=900; rest t=1800")
    r = run(p)
    m0 = CrystalModel(p)
    spread0 = _theta_spread(m0, run(crystal(steps="cc C=2 t=900")).final_state)
    spread = _theta_spread(m0, r.final_state)
    assert spread < spread0 and spread < 1e-3, (spread0, spread)
    a = r.array
    v = a[a[:, 7] == 2][:, 1]
    dv = np.diff(v)
    assert np.all(dv >= -1e-9) or np.all(dv <= 1e-9)
    m = CrystalModel(p)
    st = r.final_state
    th_mean = m.cs_mean(st) / m.cs_max
    s_mean = np.log(th_mean / (1 - th_mean))
    u_ocp = m.kin.ocp(st.x[m.nodes, 0], s_mean)
    np.testing.assert_allclose(st.x[m.nodes, 1] - st.x[m.nodes, 2], u_ocp, atol=2e-3)


def test_crystal_time_first_order():
    vals = []
    for dt in (1.0, 0.5, 0.25, 0.125):
        p = crystal(C_rate=2.0, n_steps=int(36000 / dt), newton_tol=1e-13)
        vals.append(run(p, max_steps=int(8 / dt)).final_state.x[0, 0])
    d = np.diff(vals)
    orders = np.log2(np.abs(d[:-1] / d[1:]))
    assert np.all((orders > 0.9) & (orders < 1.1)), orders


def test_crystal_electrode_mesh_second_order():
    vals = []
    for sep, nj in ((7, 26), (12, 51), (22, 101), (42, 201)):
        vals.append(run(crystal(C_rate=2.0, sep_node=sep, nj=nj), max_steps=900).final_state.x[-1, 1])
    d = np.diff(vals)
    orders = np.log2(np.abs(d[:-1] / d[1:]))
    assert np.all((orders > 1.8) & (orders < 2.6)), orders


def test_cccv_hold_after_the_crystal_surfaces_drain():
    """A fast charge from mostly empty crystals drains their surfaces, then a 4.2 V hold must converge: Newton
    divergence is judged on the damped electrode update."""
    from nmc_model.uniform import kinetics
    p = crystal(cs_init=0.1 * kinetics.cs_max(crystal()), n_steps=7200,
                steps="cc C=-2 Vmax=4.2; cv V=4.2 Imin=0.5")
    r = run(p, max_steps=200)
    assert r.exit_reason == "end_of_protocol" and r.steps < 100, (r.exit_reason, r.steps)
    assert r.array[-1, 1] == pytest.approx(4.2, abs=1e-6)


def test_crystal_columns_in_the_output(tmp_path):
    """The crystal model writes the surface and mean theta of the first and last crystals: on discharge
    lithium enters at the surface, so the surface is ahead of the mean, and the columns are in the header."""
    r = run(crystal(steps="cc C=2 t=600"))
    out = tmp_path / "o.txt"
    r.write(out)
    head = out.read_text().splitlines()[0].split()
    assert head[-4:] == ["Th_surf_sep", "Th_mean_sep", "Th_surf_col", "Th_mean_col"]
    a = r.array
    ts_sep, tm_sep, ts_col, tm_col = a[-1, -4:]
    assert ts_sep > tm_sep > 0 and ts_col > tm_col > 0
    u = run(Params(mode="corrected", steps="cc C=2 t=600")).array
    assert u.shape[1] == a.shape[1] - 4                 # the uniform output is unchanged


# ------------------------------------------------------------------ the Fortran and C++ programs
CRYSTAL_GROUPS = {"model": {"particle_model": "crystal"}, "numerics": {"mode": "corrected"}}


def _numeric(path):
    return np.array([[float(x) for x in ln.split()[1:]] for ln in path.read_text().splitlines()[2:]])


@pytest.mark.parametrize("steps,extra", [("cc C=2 Vmin=3.0; rest t=600", {}), (CYCLE, {}),
                                         ("cc C=2 Vmin=3.0; rest t=600", {"dt_mode": "adaptive"})],
                         ids=["discharge", "cycle", "adaptive"])
def test_crystal_native_programs_agree(fortran_exe, cpp_exe, run_native, tmp_path, steps, extra):
    """The Fortran and C++ programs write identical files (Li_Nernst to round-off), and Python agrees to the
    printed precision, crystal columns included."""
    from .test_cross_language import assert_same
    g = {**CRYSTAL_GROUPS, "protocol": {"steps": steps}}
    g["numerics"] = {**g["numerics"], **extra}
    f = run_native(fortran_exe, g, name="f.txt")
    c = run_native(cpp_exe, g, name="c.txt")
    assert_same(f, c)
    r = run(crystal(steps=steps, **extra))
    py = tmp_path / "p.txt"
    r.write(py)
    assert f.read_text().splitlines()[:2] == py.read_text().splitlines()[:2]
    fa, pa = _numeric(f), _numeric(py)
    assert fa.shape == pa.shape
    np.testing.assert_allclose(fa[:, :3], pa[:, :3], rtol=0, atol=1e-5)          # 5 printed decimals
    np.testing.assert_allclose(fa[:, 9:], pa[:, 9:], rtol=1e-5, atol=1e-9)      # the crystal columns


def test_crystal_input_file_runs_in_python(tmp_path):
    """particle_model = 'crystal' in &model selects the crystal model in Python's input reader too."""
    from nmc_model.inputs import load
    nml = tmp_path / "in.nml"
    nml.write_text("&model\n  particle_model = 'crystal'\n/\n&numerics\n  mode = 'corrected'\n/\n"
                   "&active\n  D_c = 1.0e-11, crystal_shape = 'slab'\n/\n&cell\n  nj_crystal = 11\n/\n")
    p, model, _ = load(nml)
    assert model == p.particle_model == "crystal"
    assert (p.D_c, p.crystal_shape, p.nj_crystal) == (1e-11, "slab", 11)
