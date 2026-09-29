"""Corrected-mode core shared by both models: log variables and Scharfetter-Gummel fluxes.

Unknowns per node are (u, phi1, phi2, s):

- u = ln(c/c_bulk): the electrolyte concentration c = c_bulk e^u is positive by construction;
- s = ln(theta/(1-theta)): the particle lithiation theta = 1/(1 + e^-s) lies in (0, 1).

With these, the OCP's Nernst terms are linear in the unknowns, and the exchange current is
bounded, so a particle approaches full or empty, and the electrolyte approaches exhaustion,
smoothly and without clipping (docs/model.md section 8). Ion fluxes use exponential fitting
(Scharfetter-Gummel), which is positivity-preserving and accurate where migration dominates.
A small background conductivity kappa_bg (the solvent's own ionic conductivity) keeps phi2
defined where the salt is exhausted.

Residual convention: R(x) = 0 with blocks A, B, D of dR/dx (A: w.r.t. the previous node,
D: the next); the Newton step solves J dx = -R.
"""
from __future__ import annotations

import math

import numpy as np

U, P1, P2, S = 0, 1, 2, 3
N = 4


def bernoulli(x):
    """B(x) = x / (e^x - 1), with its series near 0 (the same formulas in Fortran and C++)."""
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-3
    xs = np.where(small | (x > 700.0), 1.0, x)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        big = xs / (np.exp(xs) - 1.0)
        huge = x * np.exp(-np.minimum(x, 745.0) * (x > 700.0))
    return np.where(small, 1.0 - x / 2.0 + x * x / 12.0, np.where(x > 700.0, huge, big))


def bernoulli_prime(x):
    """B'(x) = B(x) (1 - B(-x)) / x, with its series near 0 (no overflow for large |x|)."""
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-3
    xs = np.where(small, 1.0, x)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        big = bernoulli(xs) * (1.0 - bernoulli(-xs)) / xs
    return np.where(small, -0.5 + x / 6.0 - x ** 3 / 180.0, big)


def softplus(x):
    """ln(1 + e^x), computed as max(x, 0) + ln(1 + e^-|x|) (the same formula in Fortran and C++)."""
    return np.maximum(x, 0.0) + np.log(1.0 + np.exp(-np.abs(x)))


def sigmoid(x):
    """1 / (1 + e^-x)."""
    return 0.5 * (1.0 + np.tanh(0.5 * x))


def logit(theta):
    return math.log(theta / (1.0 - theta))


class LogKinetics:
    """Butler-Volmer insertion kinetics in the log variables.

    U(u, s) = U_ref + (RT/F)(u - s) + RK(theta), with the Redlich-Kister coefficients `ak`;
    i0 = F k c^aa (cs_max - cs)^aa cs^ac, evaluated as a logarithm so that it stays bounded.
    """

    def __init__(self, *, U_ref, ak, cs_max, k, alpha_a, alpha_c, c_bulk, R, T, F):
        self.U_ref, self.ak, self.cs_max = U_ref, tuple(ak), cs_max
        self.aa, self.ac = alpha_a, alpha_c
        self.rtf = R * T / F
        self.f = 1.0 / self.rtf
        self.F = F
        self.ln_pre = math.log(F * k * c_bulk ** alpha_a * cs_max ** (alpha_a + alpha_c))

    def _rk(self, th):
        x = 2 * th - 1
        rk = np.zeros_like(th)
        drk = np.zeros_like(th)
        for k, a in enumerate(self.ak):
            with np.errstate(divide="ignore", invalid="ignore"):
                rk = rk + a * (x ** (k + 1) - (2 * th * k * (1 - th)) / x ** (1 - k))
            term = 2.0 * (2 * k + 1) * x ** k
            if k >= 2:
                term = term - 4.0 * k * (k - 1) * th * (1 - th) * x ** (k - 2)
            drk = drk + a * term
        return rk, drk

    def ocp(self, u, s):
        rk, _ = self._rk(sigmoid(s))
        return self.U_ref + self.rtf * (u - s) + rk

    def rate(self, x):
        """i_n [A/cm2, anodic positive] and d i_n / d(u, phi1, phi2, s), stacked on the last axis."""
        u, s = x[..., U], x[..., S]
        th = sigmoid(s)
        rk, drk = self._rk(th)
        U_ = self.U_ref + self.rtf * (u - s) + rk
        dU_ds = -self.rtf + drk * th * (1 - th)
        ln_i0 = self.ln_pre + self.aa * u - self.aa * softplus(s) - self.ac * softplus(-s)
        eta = x[..., P1] - x[..., P2] - U_
        with np.errstate(over="ignore"):
            ea = np.exp(ln_i0 + self.aa * self.f * eta)
            ec = np.exp(ln_i0 - self.ac * self.f * eta)
        i = ea - ec
        di_deta = self.aa * self.f * ea + self.ac * self.f * ec
        d_u = i * self.aa - di_deta * self.rtf
        d_s = i * (-self.aa * th + self.ac * (1 - th)) - di_deta * dU_ds
        return i, np.stack([d_u, di_deta, -di_deta, d_s], axis=-1)


class Transport:
    """A binary electrolyte (D+, D- from D0 and t+) and the background conductivity."""

    def __init__(self, *, D0, t_plus, c_bulk, kappa_bg, R, T, F):
        t_an = 1.0 - t_plus
        self.Dp = D0 * (1.0 + t_an / t_plus) / (2.0 * t_an / t_plus)
        self.Dm = self.Dp * t_an / t_plus
        self.c_bulk, self.kappa_bg, self.F = c_bulk, kappa_bg, F
        self.f = F / (R * T)

    def fluxes(self, xa_, xb_, g, gs):
        """Face fluxes from state a to state b and their derivatives.

        g = eps/(tau h) per face, gs = (1 - eps) sigma / h per face. Returns F (..., 3) with
        (N+, i1, i2), and dFa, dFb (..., 3, 4) w.r.t. the states on either side.
        """
        ca, cb = self.c_bulk * np.exp(xa_[..., U]), self.c_bulk * np.exp(xb_[..., U])
        dphi2 = xb_[..., P2] - xa_[..., P2]
        d = self.f * dphi2
        Bp, Bm = bernoulli(d), bernoulli(-d)
        dBp, dBm = bernoulli_prime(d), bernoulli_prime(-d)
        g = g * np.ones_like(ca)
        gs = gs * np.ones_like(ca)
        Np = g * self.Dp * (Bp * ca - Bm * cb)
        Nm = g * self.Dm * (Bm * ca - Bp * cb)
        dNp_dd = g * self.Dp * (dBp * ca + dBm * cb)
        dNm_dd = g * self.Dm * (-dBm * ca - dBp * cb)
        kb = g * self.kappa_bg
        i2 = self.F * (Np - Nm) - kb * dphi2
        i1 = -gs * (xb_[..., P1] - xa_[..., P1])
        Fv = np.stack([Np, i1, i2], axis=-1)
        dFa = np.zeros(ca.shape + (3, N)); dFb = np.zeros(ca.shape + (3, N))
        dFa[..., 0, U] = g * self.Dp * Bp * ca
        dFb[..., 0, U] = -g * self.Dp * Bm * cb
        dFa[..., 0, P2] = -self.f * dNp_dd
        dFb[..., 0, P2] = self.f * dNp_dd
        dFa[..., 1, P1] = gs
        dFb[..., 1, P1] = -gs
        dFa[..., 2, U] = self.F * (dFa[..., 0, U] - g * self.Dm * Bm * ca)
        dFb[..., 2, U] = self.F * (dFb[..., 0, U] + g * self.Dm * Bp * cb)
        dFa[..., 2, P2] = -self.F * self.f * (dNp_dd - dNm_dd) + kb
        dFb[..., 2, P2] = self.F * self.f * (dNp_dd - dNm_dd) - kb
        return Fv, dFa, dFb


class Electrode:
    """The electrode scale: Li foil | separator | cathode | collector, in the log variables.

    `mesh` is the finite-volume mesh (nodes 0 .. nj-1, interface node s). The S column is either
    held fixed (agglomerate model) or carries the particles' log-odds at nodes s .. nj-1 with a
    local reaction (uniform model: `particles` = (kinetics, a, vf)).
    """

    def __init__(self, mesh, *, eps, eps_sep, tau, tau_sep, sigma, transport, particles=None):
        self.m = mesh
        s, nj = mesh.s, mesh.nj
        k = np.arange(nj - 1)
        sep = k < s
        h = (mesh.dx[:-1] + mesh.dx[1:]) / 2.0
        self.g = np.where(sep, eps_sep / tau_sep, eps / tau) / h
        self.gs = np.where(sep, 1.0 - eps_sep, 1.0 - eps) * sigma / h
        self.eps_node = np.where(np.arange(nj) < s, eps_sep, eps)
        self.tr = transport
        self.particles = particles

    def residual_and_blocks(self, x, xold, dt, I):
        m, tr = self.m, self.tr
        nj, s = m.nj, m.s
        F = tr.F
        Fv, dFa, dFb = tr.fluxes(x[:-1], x[1:], self.g, self.gs)
        R = np.zeros_like(x)
        A = np.zeros((nj, N, N)); B = np.zeros_like(A); D = np.zeros_like(A)
        c, cold = tr.c_bulk * np.exp(x[:, U]), tr.c_bulk * np.exp(xold[:, U])
        # foil face: N+ = I/F, zero electronic current, phi2 = 0 (the gauge)
        R[0, U] = Fv[0, 0] - I / F
        B[0, U] = dFa[0, 0]; D[0, U] = dFb[0, 0]
        R[0, P1] = x[1, P1] - x[0, P1]
        B[0, P1, P1] = -1.0; D[0, P1, P1] = 1.0
        R[0, P2] = x[0, P2]
        B[0, P2, P2] = 1.0
        # separator, interface and cathode: flux differences plus storage
        js = np.arange(1, nj - 1)
        for f, row in enumerate((U, P1, P2)):
            R[js, row] = Fv[js, f] - Fv[js - 1, f]
            B[js, row] += dFa[js, f] - dFb[js - 1, f]
            D[js, row] += dFb[js, f]
            A[js, row] += -dFa[js - 1, f]
        R[js, U] += self.eps_node[js] * m.dx[js] * (c[js] - cold[js]) / dt
        B[js, U, U] += self.eps_node[js] * m.dx[js] * c[js] / dt
        # collector: no salt flux, no ionic current, electronic current I
        R[-1, U] = -Fv[-1, 0]; A[-1, U] = -dFa[-1, 0]; B[-1, U] = -dFb[-1, 0]
        R[-1, P1] = Fv[-1, 1] - I; A[-1, P1] = dFa[-1, 1]; B[-1, P1] = dFb[-1, 1]
        R[-1, P2] = Fv[-1, 2]; A[-1, P2] = dFa[-1, 2]; B[-1, P2] = dFb[-1, 2]
        # the S column
        R[:, S] = x[:, S] - xold[:, S]
        B[:, S, :] = 0.0; B[:, S, S] = 1.0
        if self.particles is not None:
            kin, a, vf = self.particles
            jr = np.arange(s, nj)
            i, di = kin.rate(x[jr])
            dxr = m.dx[jr]
            R[jr, U] -= a * i * dxr / F
            B[jr, U, :] -= a * di * dxr[:, None] / F
            R[jr, P1] += a * i * dxr
            B[jr, P1, :] += a * di * dxr[:, None]
            R[jr, P2] -= a * i * dxr
            B[jr, P2, :] -= a * di * dxr[:, None]
            th, thold = sigmoid(x[jr, S]), sigmoid(xold[jr, S])
            R[jr, S] = vf * kin.cs_max * (th - thold) / dt + a * i / F
            B[jr, S, :] = a * di / F
            B[jr, S, S] += vf * kin.cs_max * th * (1 - th) / dt
        return R, A, B, D


def equilibrate(A, B, D, *rhs):
    """Scale every equation by the largest entry of its row in B (the solution is unchanged)."""
    sc = np.abs(B).max(axis=-1)
    sc[sc == 0] = 1.0
    out = [A / sc[..., None], B / sc[..., None], D / sc[..., None]]
    return out + [r / sc for r in rhs]


def bounded(dxs, caps=((U, 1.0), (P1, 0.1), (P2, 0.1), (S, 2.0))):
    """Newton step length <= 1 limiting |du| <= 1, |dphi| <= 0.1 V and |ds| <= 2 per iteration."""
    lam = 1.0
    for arr in dxs:
        for col, cap in caps:
            mx = float(np.max(np.abs(arr[..., col])))
            if mx > cap:
                lam = min(lam, cap / mx)
    return lam


DIVERGED = 1.0e3   # an update this large in the O(1) log variables and potentials means no solution: stop early


def physical_update(x, dx, frozen_s=False):
    """The Newton update measured in the physical variables: max of |dc|/c_bulk = e^u |du|, |dphi| [V]
    and |dtheta| = theta(1 - theta) |ds|.

    Near an exhausted electrolyte or a full or empty particle the log variables are ill-conditioned,
    and their round-off (in u or s) is physically irrelevant; convergence is judged on c, phi and
    theta instead. With frozen_s the S column is not a particle (the agglomerate model's electrode).
    """
    th = sigmoid(x[..., S])
    parts = [np.exp(x[..., U]) * np.abs(dx[..., U]), np.abs(dx[..., P1]), np.abs(dx[..., P2])]
    if not frozen_s:
        parts.append(th * (1 - th) * np.abs(dx[..., S]))
    return max(float(np.max(p)) for p in parts)


def converged(upd, prev, tol):
    """The scaled update is below tol, or it has stagnated at the round-off floor."""
    return upd <= tol or (upd <= 1.0e3 * tol and upd >= 0.5 * prev)
