"""Independent residual of the corrected model, written face by face from docs/model.md.

It shares only the kinetics and the parameters with the package, not the assembly. The
tests compare it with the assembled right-hand side (G = -F), which catches defects in the
residual that a Jacobian check cannot see (such as deviation D-11).
"""
import numpy as np

from nmc_model.uniform import kinetics
from nmc_model.uniform.model import make_mesh
from nmc_model.uniform.params import Params


def residual(p: Params, c: np.ndarray, c_old: np.ndarray, dt: float) -> np.ndarray:
    """F(c) for one backward-Euler step, same row order and scaling as the assembly."""
    m = make_mesh(p)
    nj, s, dx = m.nj, m.s, m.dx
    F, a, I = p.F, p.spec_a, p.i_app
    t_an = 1.0 - p.t_plus
    d_cat = p.D * (1.0 + t_an / p.t_plus) / (2.0 * t_an / p.t_plus)
    d_an = d_cat * t_an / p.t_plus
    rt = p.R * p.T

    # face k joins nodes k and k+1; faces 0..s-1 lie in the separator, s..nj-2 in the cathode
    k = np.arange(nj - 1)
    sep = k < s
    e = np.where(sep, p.eps_sep, p.eps)
    tau = np.where(sep, p.tau_sep, p.tortuosity)
    e_solid = np.where(sep, 1.0 - p.eps_sep, 1.0 - p.eps)
    w = dx[:-1] / (dx[:-1] + dx[1:])
    beta = 2.0 / (dx[:-1] + dx[1:])
    cf = w * c[1:, 0] + (1.0 - w) * c[:-1, 0]
    g = beta[:, None] * (c[1:] - c[:-1])

    N_plus = -e * (d_cat / tau) * g[:, 0] - e * p.z_plus * (d_cat / rt / tau) * F * cf * g[:, 2]
    i2 = (-e * F * (p.z_plus * d_cat + p.z_minus * d_an) / tau * g[:, 0]
          - e * F ** 2 * (p.z_plus ** 2 * d_cat + p.z_minus ** 2 * d_an) / rt / tau * cf * g[:, 2])
    i1 = -e_solid * p.sigma * g[:, 1]

    i_n = kinetics.reaction_rate(p, c[:, 0], c[:, 3], c[:, 1], c[:, 2])
    dcdt = (c - c_old) / dt
    R = np.zeros_like(c)

    # Li-foil face: N+ = I/F, i1 = 0, phi2 = 0, cs frozen
    R[0, 0] = I / F - N_plus[0]
    R[0, 1] = 0.0 - i1[0]
    R[0, 2] = c[0, 2]
    R[0, 3] = -p.vf_AM * dcdt[0, 3]

    j = np.arange(1, s)                       # separator interior
    R[j, 0] = N_plus[j - 1] - N_plus[j] - p.eps_sep * dx[j] * dcdt[j, 0]
    R[j, 1] = i1[j - 1] - i1[j]
    R[j, 2] = i2[j - 1] - i2[j]
    R[j, 3] = -(1.0 - p.eps_sep) * dcdt[j, 3]

    R[s, 0] = N_plus[s - 1] - N_plus[s]        # interface (zero volume)
    R[s, 1] = i1[s - 1] - i1[s]
    R[s, 2] = i2[s - 1] - i2[s]

    j = np.arange(s + 1, nj - 1)              # cathode interior
    R[j, 0] = N_plus[j - 1] - N_plus[j] - p.eps * dx[j] * dcdt[j, 0] + a * i_n[j] * dx[j] / F
    R[j, 1] = i1[j - 1] - i1[j] - a * i_n[j] * dx[j]
    R[j, 2] = i2[j - 1] - i2[j] + a * i_n[j] * dx[j]

    n = nj - 1                                # current collector: N+ = 0, i1 = I, i2 = 0
    R[n, 0] = N_plus[n - 1]
    R[n, 1] = i1[n - 1] - I
    R[n, 2] = i2[n - 1]

    j = np.arange(s, nj)                      # uniform-particle solid balance
    R[j, 3] = -p.vf_AM * dcdt[j, 3] - a * i_n[j] / F
    return R


def faces(p: Params, c: np.ndarray):
    """Face cation flux, ionic current and anion flux (for conservation checks)."""
    m = make_mesh(p)
    nj, s, dx = m.nj, m.s, m.dx
    F = p.F
    t_an = 1.0 - p.t_plus
    d_cat = p.D * (1.0 + t_an / p.t_plus) / (2.0 * t_an / p.t_plus)
    d_an = d_cat * t_an / p.t_plus
    rt = p.R * p.T
    k = np.arange(nj - 1)
    e = np.where(k < s, p.eps_sep, p.eps)
    tau = np.where(k < s, p.tau_sep, p.tortuosity)
    w = dx[:-1] / (dx[:-1] + dx[1:])
    beta = 2.0 / (dx[:-1] + dx[1:])
    cf = w * c[1:, 0] + (1.0 - w) * c[:-1, 0]
    g = beta[:, None] * (c[1:] - c[:-1])
    N_plus = -e * d_cat / tau * g[:, 0] - e * p.z_plus * d_cat / rt / tau * F * cf * g[:, 2]
    N_minus = -e * d_an / tau * g[:, 0] - e * p.z_minus * d_an / rt / tau * F * cf * g[:, 2]
    i2 = F * (p.z_plus * N_plus + p.z_minus * N_minus)
    return N_plus, i2, N_minus


# ------------------------------------------------------------------ corrected mode (log variables)

def _bern(x):
    return 1.0 - x / 2.0 + x * x / 12.0 if abs(x) < 1e-3 else x / np.expm1(x)


def to_physical(p: Params, x: np.ndarray) -> np.ndarray:
    """A corrected-mode state (u, phi1, phi2, s) as (c, phi1, phi2, cs)."""
    out = x.copy()
    out[:, 0] = p.c_bulk * np.exp(x[:, 0])
    out[:, 3] = kinetics.cs_max(p) / (1.0 + np.exp(-x[:, 3]))
    return out


def sg_residual(p: Params, x: np.ndarray, x_old: np.ndarray, dt: float, I: float) -> np.ndarray:
    """The corrected model's residual written node by node from docs/model.md section 8.

    Scharfetter-Gummel ion fluxes, the background conductivity, the package's c-form kinetics.
    Same rows, signs and scaling as nmc_model.logcore.Electrode.
    """
    m = make_mesh(p)
    nj, s, dx = m.nj, m.s, m.dx
    F, f = p.F, p.F / (p.R * p.T)
    t_an = 1.0 - p.t_plus
    Dp = p.D * (1.0 + t_an / p.t_plus) / (2.0 * t_an / p.t_plus)
    Dm = Dp * t_an / p.t_plus
    y, yo = to_physical(p, x), to_physical(p, x_old)
    c, cs, co, cso = y[:, 0], y[:, 3], yo[:, 0], yo[:, 3]
    Np, i1, i2 = np.zeros(nj - 1), np.zeros(nj - 1), np.zeros(nj - 1)
    for k in range(nj - 1):
        sep = k < s
        eps, tau = (p.eps_sep, p.tau_sep) if sep else (p.eps, p.tortuosity)
        h = (dx[k] + dx[k + 1]) / 2.0
        g = eps / tau / h
        d = f * (x[k + 1, 2] - x[k, 2])
        Np[k] = g * Dp * (_bern(d) * c[k] - _bern(-d) * c[k + 1])
        Nm = g * Dm * (_bern(-d) * c[k] - _bern(d) * c[k + 1])
        i2[k] = F * (Np[k] - Nm) - g * p.kappa_bg * (x[k + 1, 2] - x[k, 2])
        i1[k] = -(1.0 - (p.eps_sep if sep else p.eps)) * p.sigma * (x[k + 1, 1] - x[k, 1]) / h
    i_n = kinetics.reaction_rate(p, c, cs, x[:, 1], x[:, 2])
    a = p.spec_a
    R = np.zeros_like(x)
    R[0] = [Np[0] - I / F, x[1, 1] - x[0, 1], x[0, 2], x[0, 3] - x_old[0, 3]]
    for j in range(1, nj - 1):
        eps = p.eps_sep if j < s else p.eps
        R[j, 0] = eps * dx[j] * (c[j] - co[j]) / dt + Np[j] - Np[j - 1]
        R[j, 1] = i1[j] - i1[j - 1]
        R[j, 2] = i2[j] - i2[j - 1]
        R[j, 3] = x[j, 3] - x_old[j, 3]
    R[-1, :3] = [-Np[-1], i1[-1] - I, i2[-1]]
    R[-1, 3] = 0.0
    for j in range(s, nj):                    # particles, with their reaction in the electrode rows
        R[j, 0] -= a * i_n[j] * dx[j] / F
        R[j, 1] += a * i_n[j] * dx[j]
        R[j, 2] -= a * i_n[j] * dx[j]
        R[j, 3] = p.vf_AM * (cs[j] - cso[j]) / dt + a * i_n[j] / F
    return R
