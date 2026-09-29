"""Open-circuit potential and Butler-Volmer kinetics of the uniform-particle NMC111 model.

The OCP is a Redlich-Kister expansion plus a Nernst term in the electrolyte concentration
(docs/model.md section 2; deviation D-16 concerns the Nernst term). All functions are
vectorized. In faithful mode U_ref, the coefficients and x_max = 0.55 are float32-rounded,
the Redlich-Kister sum is accumulated in single precision, integer powers follow libgcc's
__powidf2, and the exchange current is rounded to float32 on every evaluation (D-4).
"""
from __future__ import annotations

import numpy as np

from ..fnum import powi, powr
from .params import Params, f32

_U_REF = 3.8685682447595453
_AK = (-0.2018059457910574, 0.1123408808723528, -0.0483699097647364, 0.0231624989428732,
       -0.0377897311905149, -0.3307806975105846, 0.2392976745148739, 0.7787126945566982,
       -0.2599451275866008, -0.5898456896544948, 0.0520147453263591)
_X_MAX = 0.55


def _faithful(p: Params) -> bool:
    return p.mode == "faithful"


def _x_max(p: Params) -> float:
    return f32(_X_MAX) if _faithful(p) else _X_MAX


def cs_max(p: Params) -> float:
    """Maximum lithium concentration in the active material [mol/cm3]."""
    return _x_max(p) * p.mol_vol


def theta(p: Params, cs):
    """Lithiation fraction x/x_max, with x = cs/mol_vol."""
    return (cs / p.mol_vol) / _x_max(p)


def _rk(p: Params, th):
    """Redlich-Kister sum (accumulated in single precision in faithful mode)."""
    faithful = _faithful(p)
    ak = [f32(a) for a in _AK] if faithful else list(_AK)
    x = 2 * th - 1
    v = np.zeros_like(th, dtype=np.float32) if faithful else np.zeros_like(th)
    for k, a in enumerate(ak):
        with np.errstate(divide="ignore", invalid="ignore"):
            term = a * (powi(x, k + 1) - (2 * th * k * (1 - th)) / powi(x, 1 - k))
        w = v.astype(np.float64) + term
        v = w.astype(np.float32) if faithful else w
    return v.astype(np.float64)


def ocp(p: Params, c, cs):
    """U(theta, c) = U_ref + (RT/F) ln[(c/c_bulk)(1-theta)/theta] + Redlich-Kister sum [V]."""
    th = theta(p, cs)
    u_ref = f32(_U_REF) if _faithful(p) else _U_REF
    with np.errstate(divide="ignore", invalid="ignore"):
        nernst = p.R * p.T / p.F * np.log(c / p.c_bulk * (1.0 - th) / th)
    return u_ref + nernst + _rk(p, th)


def ocp_slopes(p: Params, c, cs):
    """(dU/dc, dU/dcs) for corrected mode."""
    th = theta(p, cs)
    rtf = p.R * p.T / p.F
    x = 2 * th - 1
    drk = np.zeros_like(th)
    for k, a in enumerate(_AK):
        term = 2.0 * (2 * k + 1) * x ** k
        if k >= 2:
            term = term - 4.0 * k * (k - 1) * th * (1 - th) * x ** (k - 2)
        drk = drk + a * term
    with np.errstate(divide="ignore", invalid="ignore"):
        du_dth = rtf * (-1.0 / (1.0 - th) - 1.0 / th) + drk
        return rtf / c, du_dth / (p.mol_vol * _x_max(p))


THETA_REG = 1.0e-6   # corrected mode: below this lithiation (or vacancy) fraction, x^alpha is regularized


def _power_reg(x, alpha, delta):
    """x**alpha for x >= delta; below it the C1 quadratic delta**alpha*((2-alpha)u + (alpha-1)u**2), u = x/delta.

    It keeps g(0) = 0 with a finite slope, so Newton's method stays well posed when a particle
    empties or fills completely (deviation D-13). Returns (g, dg/dx).
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        u = x / delta
        lo = x < delta
        g = np.where(lo, delta ** alpha * ((2.0 - alpha) * u + (alpha - 1.0) * u * u), x ** alpha)
        dg = np.where(lo, delta ** (alpha - 1.0) * ((2.0 - alpha) + 2.0 * (alpha - 1.0) * u), alpha * x ** (alpha - 1.0))
    return g, dg


def exchange_current(p: Params, c, cs):
    """i0 = F k c^aa (cs_max - cs)^aa cs^ac [A/cm2]."""
    if p.mode == "faithful":
        with np.errstate(invalid="ignore"):
            i0 = p.F * p.k_rxn * powr(c, p.alpha_a) * powr(cs_max(p) - cs, p.alpha_a) * powr(cs, p.alpha_c)
        return np.asarray(i0, dtype=np.float32).astype(np.float64)
    return exchange_current_and_slope(p, c, cs)[0]


def exchange_current_and_slope(p: Params, c, cs):
    """Corrected-mode i0 and d(i0)/dcs, with the solid-concentration powers regularized (D-13)."""
    delta = THETA_REG * cs_max(p)
    gv, dgv = _power_reg(cs_max(p) - cs, p.alpha_a, delta)
    gs, dgs = _power_reg(cs, p.alpha_c, delta)
    with np.errstate(invalid="ignore"):
        pre = p.F * p.k_rxn * (c ** p.alpha_a)
    return pre * gv * gs, pre * (gs * -dgv + gv * dgs)


def reaction_rate(p: Params, c, cs, phi1, phi2):
    """Butler-Volmer current density per unit interfacial area [A/cm2]; anodic positive."""
    eta = phi1 - phi2 - ocp(p, c, cs)
    i0 = exchange_current(p, c, cs)
    rt = p.R * p.T
    with np.errstate(invalid="ignore", over="ignore"):
        return i0 * (np.exp(p.alpha_a * p.F * eta / rt) - np.exp(-(p.alpha_c * p.F * eta / rt)))


def reaction_derivatives(p: Params, c, cs, phi1, phi2):
    """Rate and its finite-difference derivatives (d/dc, d/dcs, d/dphi1, d/dphi2).

    Central differences with an absolute step, except forward differences when
    c or cs is at or below the step (as in the original, deviation D-6).
    """
    h = p.fd_step
    f = lambda c_, cs_, p1_, p2_: reaction_rate(p, c_, cs_, p1_, p2_)
    i = f(c, cs, phi1, phi2)
    with np.errstate(invalid="ignore"):
        d_c = np.where(c <= h,
                       (f(c + h, cs, phi1, phi2) - i) / h,
                       (f(c + h, cs, phi1, phi2) - f(c - h, cs, phi1, phi2)) / (2.0 * h))
        d_cs = np.where(cs <= h,
                        (f(c, cs + h, phi1, phi2) - i) / h,
                        (f(c, cs + h, phi1, phi2) - f(c, cs - h, phi1, phi2)) / (2.0 * h))
        d_p1 = (f(c, cs, phi1 + h, phi2) - f(c, cs, phi1 - h, phi2)) / (2.0 * h)
        d_p2 = (f(c, cs, phi1, phi2 + h) - f(c, cs, phi1, phi2 - h)) / (2.0 * h)
    return i, d_c, d_cs, d_p1, d_p2


def reaction_derivatives_analytic(p: Params, c, cs, phi1, phi2):
    """Rate and exact derivatives (d/dc, d/dcs, d/dphi1, d/dphi2); used in corrected mode (fixes D-6)."""
    rt = p.R * p.T
    A_, B_ = p.alpha_a * p.F / rt, p.alpha_c * p.F / rt
    eta = phi1 - phi2 - ocp(p, c, cs)
    i0, di0_dcs = exchange_current_and_slope(p, c, cs)
    du_dc, du_dcs = ocp_slopes(p, c, cs)
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        ea, ec = np.exp(A_ * eta), np.exp(-B_ * eta)
        i = i0 * (ea - ec)
        di_deta = i0 * (A_ * ea + B_ * ec)
        d_c = p.alpha_a * i / c - di_deta * du_dc
        d_cs = di0_dcs * (ea - ec) - di_deta * du_dcs
    return i, d_c, d_cs, di_deta, -di_deta


def li_exchange_current(p: Params, c):
    """Exchange current density of the lithium counter electrode [A/cm2]."""
    return p.F * p.k_Li * (c ** 0.5) * (p.c_Li_ref ** 0.5)


def li_foil(p: Params, c, I):
    """The lithium foil's Nernst potential U_Li and overpotential eta_Li [V] (corrected mode).

    With the foil metal as the 0 V reference, phi2 at the foil is -(U_Li + eta_Li), where
    U_Li = (RT/F) ln(c/c_Li_ref) and eta_Li = (RT/(alpha F)) asinh(I/(2 i0)) (symmetric
    Butler-Volmer, alpha = 0.5; I > 0 on discharge, when the foil is oxidized).
    """
    alpha = 0.5
    rtf = p.R * p.T / p.F
    i0 = li_exchange_current(p, c)
    return rtf * np.log(c / p.c_Li_ref), rtf / alpha * np.arcsinh(I / (2.0 * i0))
