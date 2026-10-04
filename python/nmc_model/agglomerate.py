"""Agglomerate model: porous electrode with porous spherical agglomerates (docs/model.md section 4).

Faithful mode reproduces the original program `NMC111_agg.f95`: the same equations, the same
single-precision constants, the same order of floating-point operations where rounding depends
on it, one linearized solve per step, and the electrode and agglomerate scales solved one after
the other.

Arrays: the electrode state `c` has shape (nj, 4) with columns (c, phi1, phi2, cs); the
agglomerate state `ca` has shape (nj, nja, 4) with the same columns along each agglomerate
radius (node 0 = center, nja-1 = surface). Agglomerates are solved at nodes s..nj-1, where
s = sep_node - 1 is the separator/cathode interface.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np

import bandsolver

from .fnum import f32, powi, powr

C, P1, P2, CS = 0, 1, 2, 3
N = 4


# ---------------------------------------------------------------------------------- parameters

@dataclass(frozen=True)
class AggParams:
    """Parameters of the agglomerate model (docs/parameters.md).

    Defaults are generic, documented values. The original program's fitted values (agglomerate
    diffusivity, rate constant, electrode tortuosity, loading) are not distributed; faithful mode
    takes them as arguments (see `faithful`).
    """
    # run-specific values (substituted by the original run generator)
    nj: int = 75
    nja: int = 33                      # agglomerate nodes
    C_rate: float = 1.0
    eps: float = 0.4                   # cathode porosity
    thickness_um: float = 100.0        # cathode thickness [um]
    sigma: float = 0.1                 # electronic conductivity [S/cm]
    time_mod: float = 20.0             # number of steps = 3600 * C_rate * time_mod
    # fixed values of the template
    sep_node: int = 22
    L_sep: float = 25.0e-4
    eps_sep: float = 0.39
    tau_sep: float = 4.0
    D0: float = 2.89e-6                # salt diffusivity [cm2/s]
    t_plus: float = 0.375
    tortuosity: Optional[float] = None # faithful mode only: electrode tortuosity (D_e = eps D0 / tortuosity)
    c_bulk: float = 1.0e-3
    c0_init: float = 1.0e-3
    mass_loading: float = 0.020        # [g/cm2]
    percent_active: float = 0.95
    M: float = 96.46
    rho: float = 4.7
    mol_vol: float = 0.0476881609      # a literal in the original, not rho/M
    Q_th: float = 0.155
    R_agg: float = 1.0e-4 * 5.0        # agglomerate radius [cm]
    R_xtal: float = 200.0e-7           # crystal radius inside agglomerates [cm]
    eps_agg: float = 0.2
    D_agg: Optional[float] = None      # faithful mode only: salt diffusivity in agglomerates [cm2/s]
    rxn_k: float = 2.5e-6              # rate constant (i0 ~ 0.1 mA/cm2 at theta = 0.5, 1 M)
    alpha_a: float = 0.5
    alpha_c: float = 0.5
    phi1_init: float = 4.3
    cs_init: float = 1.0e-5
    t_max: float = 72000.0
    R: float = 8.314
    T: float = 298.0
    F: float = 96485.0
    PI: float = 3.141592654
    ramp_start: float = 50.0           # the current starts at I/ramp_start ...
    ramp_factor: float = 1.5           # ... and grows by this factor per step
    fd_step: float = 1.0e-6
    # --- corrected mode only (docs/model.md section 6) ---
    bruggeman: float = -0.5            # electrode tortuosity = eps**bruggeman (the fitted `tortuosity` is faithful only)
    tau_agg: Optional[float] = None    # agglomerate-pore tortuosity; None: eps_agg**-0.5 (Bruggeman)
    sigma_agg: Optional[float] = None  # agglomerate solid conductivity; None: sigma
    k_Li: float = 1.0e-6               # lithium foil rate constant
    c_Li_ref: float = 1.0e-3
    dt_s: float = 1.0                  # time step [s]
    kappa_bg: float = 1.0e-8           # background (solvent) ionic conductivity [S/cm], keeps phi2 defined where c = 0
    newton_tol: float = 1.0e-10
    newton_max_iter: int = 25
    V_min: float = 3.0                 # [V]
    V_max: float = 4.4                 # [V]
    steps: str = ""                    # protocol (docs/protocol.md); empty: discharge at C_rate
    cycles: int = 1
    write_interval: float = 18.0       # [s]
    mode: str = "faithful"

    @classmethod
    def faithful(cls, *, tortuosity: float, mass_loading: float, D_agg_mult: Optional[float] = None,
                 k_exp: Optional[float] = None, D_agg: Optional[float] = None, rxn_k: Optional[float] = None,
                 **run) -> "AggParams":
        """Parameters exactly as the original stored them (single-precision literals, D-4).

        The original's fitted values are required and are not part of this package: the electrode
        tortuosity, the mass loading [g/cm2], the agglomerate diffusivity (D_agg, or D_agg_mult with
        D_agg = 1e-10 * float32(D_agg_mult) as the original computed it) and the rate constant
        (rxn_k, or k_exp with k = float32(10**-float32(k_exp))).
        """
        if (D_agg is None) == (D_agg_mult is None) or (rxn_k is None) == (k_exp is None):
            raise TypeError("give D_agg or D_agg_mult, and rxn_k or k_exp")
        if D_agg is None:
            D_agg = 1.0e-10 * f32(D_agg_mult)
        if rxn_k is None:
            rxn_k = f32(10.0 ** (-f32(k_exp)))
        p = cls(tortuosity=tortuosity, mass_loading=mass_loading, **run)
        th = p.thickness_um
        return replace(
            p,
            C_rate=f32(p.C_rate), eps=f32(p.eps), sigma=f32(p.sigma), time_mod=f32(p.time_mod),
            eps_sep=f32(p.eps_sep), c_bulk=f32(p.c_bulk), c0_init=f32(p.c0_init),
            mass_loading=f32(p.mass_loading), percent_active=f32(p.percent_active),
            M=f32(p.M), rho=f32(p.rho), mol_vol=f32(p.mol_vol), Q_th=f32(p.Q_th),
            eps_agg=f32(p.eps_agg), D_agg=D_agg, rxn_k=rxn_k, phi1_init=f32(p.phi1_init),
            R=f32(p.R), PI=f32(p.PI), mode="faithful",
            thickness_um=th,
        )

    # derived values, in the original's order of evaluation
    @property
    def L_cath(self) -> float:
        """THICKNESS/10000.0 is a single-precision division in the original."""
        return f32(self.thickness_um / 10000.0) if self.mode == "faithful" else self.thickness_um * 1.0e-4

    @property
    def v_AM(self) -> float:
        return self.percent_active * self.mass_loading / (self.rho * self.L_cath)

    @property
    def spec_a(self) -> float:
        """Electrode-scale area: agglomerates treated as solid spheres (D-15)."""
        return 3.0 * self.v_AM / self.R_agg

    @property
    def spec_a_agg(self) -> float:
        return 3.0 * (1.0 - self.eps_agg) / self.R_xtal

    @property
    def n_steps(self) -> int:
        return int(3.6e3 * self.C_rate * self.time_mod)

    @property
    def i_specific(self) -> float:
        return self.Q_th * self.C_rate

    @property
    def dt(self) -> float:
        """Corrected-mode time step (faithful mode computes its own from t_max and n_steps)."""
        return self.dt_s

    @property
    def mass_area(self) -> float:
        """Active-material loading [g/cm2]."""
        return self.L_cath * self.v_AM * self.rho

    @property
    def i_1C(self) -> float:
        return self.Q_th * self.mass_area

    @property
    def i_final(self) -> float:
        """Applied current density after the ramp [A/cm2]."""
        return self.i_specific * self.L_cath * self.v_AM * self.rho


# ---------------------------------------------------------------------------------- kinetics

# OCP fit of the agglomerate template: U_ref and 12 Redlich-Kister coefficients
_U_REF = 3.8637058886774844
_AK = (-0.255139064974728, 0.0691287746986728, -0.1178158454270744, -0.0444434841626702,
       0.243569591966704, 0.0775338354167729, -1.0934643144519782, -0.8893166395840808,
       1.7690915896916977, 1.8213923583001588, -1.2074949744867922, -1.3952076583801158)


class Kinetics:
    def __init__(self, p: AggParams):
        self.p = p
        faithful = p.mode == "faithful"
        self.U_ref = f32(_U_REF) if faithful else _U_REF        # default-real literals (D-4)
        self.AK = [f32(a) for a in _AK] if faithful else list(_AK)
        self.x_max = p.M * p.Q_th * 3600 / p.F
        self.cs_max = p.mol_vol * p.M * p.Q_th * 3600 / p.F
        self.rt = p.R * p.T

    def ocp(self, c0, cs):
        p = self.p
        theta = (cs / p.mol_vol) / self.x_max
        x = 2 * theta - 1
        vint = np.zeros_like(theta, dtype=np.float32) if p.mode == "faithful" else np.zeros_like(theta)
        for k, ak in enumerate(self.AK):
            with np.errstate(divide="ignore", invalid="ignore"):
                term = ak * (powi(x, k + 1) - (2 * theta * k * (1 - theta)) / powi(x, 1 - k))
            v = vint.astype(np.float64) + term
            vint = v.astype(np.float32) if p.mode == "faithful" else v   # Vint is REAL (D-4)
        with np.errstate(divide="ignore", invalid="ignore"):
            nernst = self.rt / p.F * np.log(c0 / p.c_bulk * (1.0 - theta) / theta)
        return self.U_ref + nernst + vint.astype(np.float64)

    def exchange_current(self, c0, cs):
        p = self.p
        with np.errstate(invalid="ignore"):
            return p.F * p.rxn_k * powr(c0, p.alpha_a) * powr(self.cs_max - cs, p.alpha_a) * powr(cs, p.alpha_c)

    def rate(self, c0, cs, phi1, phi2):
        p = self.p
        eta = phi1 - phi2 - self.ocp(c0, cs)
        i0 = self.exchange_current(c0, cs)
        if p.mode == "faithful":
            i0 = f32(i0)                                   # ex_curr is REAL (D-4)
        with np.errstate(over="ignore", invalid="ignore"):
            return i0 * (np.exp(p.alpha_a * p.F * eta / self.rt) - np.exp(-(p.alpha_c * p.F * eta / self.rt)))

    def rate_and_fd(self, c0, cs, phi1, phi2):
        """Rate and finite-difference derivatives (d/dc, d/dphi1, d/dphi2, d/dcs) (D-6)."""
        h = self.p.fd_step
        f = self.rate
        i = f(c0, cs, phi1, phi2)
        with np.errstate(invalid="ignore"):
            d_c = np.where(c0 <= h, (f(c0 + h, cs, phi1, phi2) - i) / h,
                           (f(c0 + h, cs, phi1, phi2) - f(c0 - h, cs, phi1, phi2)) / (2.0 * h))
            d_cs = np.where(cs <= h, (f(c0, cs + h, phi1, phi2) - i) / h,
                            (f(c0, cs + h, phi1, phi2) - f(c0, cs - h, phi1, phi2)) / (2.0 * h))
            d_p1 = (f(c0, cs, phi1 + h, phi2) - f(c0, cs, phi1 - h, phi2)) / (2.0 * h)
            d_p2 = (f(c0, cs, phi1, phi2 + h) - f(c0, cs, phi1, phi2 - h)) / (2.0 * h)
        return i, np.stack([d_c, d_p1, d_p2, d_cs], axis=-1)


# ---------------------------------------------------------------------------------- grids

def _faces(dx):
    """Face interpolation weights and gradient factors for a line of control volumes."""
    aW = np.zeros_like(dx); bW = np.zeros_like(dx); aE = np.zeros_like(dx); bE = np.zeros_like(dx)
    with np.errstate(invalid="ignore", divide="ignore"):
        aW[..., 1:] = dx[..., :-1] / (dx[..., :-1] + dx[..., 1:])
        bW[..., 1:] = 2.0 / (dx[..., :-1] + dx[..., 1:])
        aE[..., :-1] = dx[..., :-1] / (dx[..., 1:] + dx[..., :-1])
        bE[..., :-1] = 2.0 / (dx[..., :-1] + dx[..., 1:])
    return aW, aE, bW, bE


def _face_values(c, aW, aE, bW, bE):
    """West/east face values and gradients of state c (..., n, N) along axis -2."""
    cW = np.zeros_like(c); cE = np.zeros_like(c); gW = np.zeros_like(c); gE = np.zeros_like(c)
    cW[..., 1:, :] = aW[..., 1:, None] * c[..., 1:, :] + (1.0 - aW[..., 1:, None]) * c[..., :-1, :]
    gW[..., 1:, :] = bW[..., 1:, None] * (c[..., 1:, :] - c[..., :-1, :])
    cE[..., :-1, :] = aE[..., :-1, None] * c[..., 1:, :] + (1.0 - aE[..., :-1, None]) * c[..., :-1, :]
    gE[..., :-1, :] = bE[..., :-1, None] * (c[..., 1:, :] - c[..., :-1, :])
    return cW, cE, gW, gE


def _blocks(rj, dW, dE, fW, fE, aW, aE, bW, bE):
    """Control-volume coefficients -> BAND blocks along axis -3 (A dc[j-1] + B dc[j] + D dc[j+1] = G)."""
    aW, aE, bW, bE = (v[..., None, None] for v in (aW, aE, bW, bE))
    A = (1.0 - aW) * fW - bW * dW
    B = rj + bW * dW + aW * fW - (1.0 - aE) * fE + bE * dE
    D = -(aE * fE) - bE * dE
    B[..., 0, :, :] = rj[..., 0, :, :] - (1.0 - aE[..., 0, :, :]) * fE[..., 0, :, :] + bE[..., 0, :, :] * dE[..., 0, :, :]
    B[..., -1, :, :] = rj[..., -1, :, :] + bW[..., -1, :, :] * dW[..., -1, :, :] + aW[..., -1, :, :] * fW[..., -1, :, :]
    A[..., 0, :, :] = 0.0
    D[..., -1, :, :] = 0.0
    return A, B, D


class Model:
    """Assembly of both scales for one time step."""

    def __init__(self, p: AggParams):
        self.p = p
        self.kin = Kinetics(p)
        nj, s = p.nj, p.sep_node - 1
        self.s = s
        h_sep = p.L_sep / float(p.sep_node - 2)
        h_cat = p.L_cath / float(nj - p.sep_node - 1)
        dx = np.zeros(nj)
        dx[1:s] = h_sep
        dx[s + 1:nj - 1] = h_cat
        self.dx = dx
        self.faces = _faces(dx)
        # agglomerate grid
        na = p.nja
        h_c = p.R_agg / float(na - 2)
        xa = np.zeros(na)
        j = np.arange(na)
        xa[1:na - 1] = h_c * (j[1:na - 1]).astype(float) - h_c / 2.0
        xa[na - 1] = p.R_agg
        dxa = np.zeros(na)
        dxa[1:na - 1] = h_c
        self.xa, self.dxa = xa, dxa
        self.faces_a = _faces(dxa)
        faithful = p.mode == "faithful"
        pi4 = float(np.float32(4.0) * np.float32(p.PI)) if faithful else 4.0 * p.PI
        pi43 = float(np.float32(pi4) / np.float32(3.0)) if faithful else 4.0 * p.PI / 3.0
        r_W = xa - dxa / 2.0
        r_E = xa + dxa / 2.0
        self.A_W = pi4 * powr(r_W, 2.0)
        self.A_E = pi4 * powr(r_E, 2.0)
        self.dV = pi43 * (powr(r_E, 3.0) - powr(r_W, 3.0))
        # coupling factor V_agg / A_agg with the original's single-precision 4.0/3.0*PI and 4.0*PI
        pi43b = float(np.float32(np.float32(4.0) / np.float32(3.0)) * np.float32(p.PI)) if faithful else 4.0 / 3.0 * p.PI
        self.vol_agg = pi43b * powi(p.R_agg, 3)
        self.area_agg = pi4 * powi(p.R_agg, 2)
        # transport coefficients of the electrode (evaluated as in fillmat)
        t_an = 1.0 - p.t_plus
        diff_e = p.eps * p.D0 / p.tortuosity
        self.D_cat = diff_e * (1.0 + (t_an / p.t_plus)) / (2.0 * t_an / p.t_plus)
        self.D_an = self.D_cat * t_an / p.t_plus
        rt = p.R * p.T
        self.u0 = p.D0 / rt
        self.u_cat = self.D_cat / rt
        self.u_an = self.D_an / rt
        self.u_agg = p.D_agg / rt

    # -------------------------------------------------------------------- electrode scale
    def assemble_electrode(self, c, ca, dt, I):
        p, s, dx = self.p, self.s, self.dx
        nj = p.nj
        F, a, eps, eps_s, sig = p.F, p.spec_a, p.eps, p.eps_sep, p.sigma
        z, zn = 1.0, -1.0
        F2 = F * F
        aW, aE, bW, bE = self.faces
        cW, cE, gW, gE = _face_values(c, aW, aE, bW, bE)
        # kinetics use the electrode potentials and the agglomerate-surface crystal concentration
        cs_k = ca[:, -1, CS]
        i, dI = self.kin.rate_and_fd(c[:, C], cs_k, c[:, P1], c[:, P2])

        dW = np.zeros((nj, N, N)); dE = np.zeros((nj, N, N)); fW = np.zeros((nj, N, N)); fE = np.zeros((nj, N, N))
        rj = np.zeros((nj, N, N)); g = np.zeros((nj, N))
        Ds, us = p.D0 / p.tau_sep, self.u0 / p.tau_sep     # separator: the same D0 for both ions (D-14)
        k_sep = z ** 2 * us + zn ** 2 * us
        k_cat = z ** 2 * self.u_cat + zn ** 2 * self.u_an

        # node 0: Li-foil face
        dE[0, C, C] = -(eps_s * Ds)
        fE[0, C, C] = -(eps_s * z * us * F * gE[0, P2])
        dE[0, C, P2] = -(eps_s * z * us * F * cE[0, C])
        g[0, C] = -I / F + (dE[0, C, C] * gE[0, C] + fE[0, C, C] * cE[0, C])
        rj[0, CS, CS] = 0.0 - 1.0 * (1.0 - eps_s) / dt
        dE[0, P1, P1] = -(1.0 - eps_s) * sig
        g[0, P1] = 0.0 - dE[0, P1, P1] * gE[0, P1]                     # D-1
        rj[0, P2, P2] = 1.0
        g[0, P2] = 0.0 - c[0, P2]

        # separator interior
        js = slice(1, s)
        dW[js, C, C] = -(eps_s * Ds); dE[js, C, C] = -(eps_s * Ds)
        fW[js, C, C] = -(eps_s * z * us * F * gW[js, P2]); fE[js, C, C] = -(eps_s * z * us * F * gE[js, P2])
        dW[js, C, P2] = -(eps_s * z * us * F * cW[js, C]); dE[js, C, P2] = -(eps_s * z * us * F * cE[js, C])
        rj[js, C, C] = -(eps_s / dt * dx[js])
        g[js, C] = 0.0 - (fW[js, C, C] * cW[js, C] + dW[js, C, C] * gW[js, C]) + (fE[js, C, C] * cE[js, C] + dE[js, C, C] * gE[js, C])
        rj[js, CS, CS] = -(1.0 * (1.0 - eps_s) / dt)
        dW[js, P1, P1] = -(1.0 - eps_s) * sig; dE[js, P1, P1] = -(1.0 - eps_s) * sig
        g[js, P1] = 0.0 - (fW[js, P1, P1] * cW[js, P1] + dW[js, P1, P1] * gW[js, P1]) + (fE[js, P1, P1] * cE[js, P1] + dE[js, P1, P1] * gE[js, P1])
        dW[js, P2, C] = -(eps_s * F * (z * Ds + zn * Ds)); dE[js, P2, C] = -(eps_s * F * (z * Ds + zn * Ds))
        ksep_b = z ** 2 * self.u0 + zn ** 2 * self.u0
        fW[js, P2, C] = -((eps_s / p.tau_sep) * F2 * ksep_b * gW[js, P2]); fE[js, P2, C] = -((eps_s / p.tau_sep) * F2 * ksep_b * gE[js, P2])
        dW[js, P2, P2] = -(eps_s * F2 * k_sep * cW[js, C]); dE[js, P2, P2] = -(eps_s * F2 * k_sep * cE[js, C])
        g[js, P2] = 0.0 - (fW[js, P2, P2] * cW[js, P2] + dW[js, P2, P2] * gW[js, P2]) + (fE[js, P2, P2] * cE[js, P2] + dE[js, P2, P2] * gE[js, P2])

        # separator/cathode interface
        ji = slice(s, s + 1)
        dW[ji, C, C] = -(eps_s * Ds); dE[ji, C, C] = -self.D_cat
        fW[ji, C, C] = -(eps_s * z * us * F * gW[ji, P2]); fE[ji, C, C] = -(z * self.u_cat * F * gE[ji, P2])
        dW[ji, C, P2] = -(eps_s * z * us * F * cW[ji, C]); dE[ji, C, P2] = -(z * self.u_cat * F * cE[ji, C])
        g[ji, C] = 0.0 - (fW[ji, C, C] * cW[ji, C] + dW[ji, C, C] * gW[ji, C]) + (fE[ji, C, C] * cE[ji, C] + dE[ji, C, C] * gE[ji, C])
        dW[ji, P1, P1] = -(1.0 - eps) * sig; dE[ji, P1, P1] = -(1.0 - eps) * sig
        g[ji, P1] = 0.0 - (fW[ji, P1, P1] * cW[ji, P1] + dW[ji, P1, P1] * gW[ji, P1]) + (fE[ji, P1, P1] * cE[ji, P1] + dE[ji, P1, P1] * gE[ji, P1])
        dW[ji, P2, C] = -(eps_s * F * (z * Ds + zn * Ds)); dE[ji, P2, C] = -(F * (z * self.D_cat + zn * self.D_an))
        fW[ji, P2, C] = -(eps_s * F2 * k_sep * gW[ji, P2]); fE[ji, P2, C] = -(F2 * k_cat * gE[ji, P2])
        dW[ji, P2, P2] = -(eps_s * F2 * k_sep * cW[ji, C]); dE[ji, P2, P2] = -(F2 * k_cat * cE[ji, C])
        g[ji, P2] = 0.0 - (fW[ji, P2, P2] * cW[ji, P2] + dW[ji, P2, P2] * gW[ji, P2]) + (fE[ji, P2, P2] * cE[ji, P2] + dE[ji, P2, P2] * gE[ji, P2])

        # solid balance at the interface, cathode interior and collector
        jr = slice(s, nj)
        rj[jr, CS, :] = -(a * dI[jr] / F)
        rj[jr, CS, CS] = -(a * dI[jr, CS] / F) - 1.0 * p.v_AM / dt
        g[jr, CS] = +(a * i[jr] / F)

        # cathode interior
        jc = slice(s + 1, nj - 1)
        dxc = dx[jc]
        dW[jc, C, C] = -self.D_cat; dE[jc, C, C] = -self.D_cat
        fW[jc, C, C] = -(z * self.u_cat * F * gW[jc, P2]); fE[jc, C, C] = -(z * self.u_cat * F * gE[jc, P2])
        dW[jc, C, P2] = -(z * self.u_cat * F * cW[jc, C]); dE[jc, C, P2] = -(z * self.u_cat * F * cE[jc, C])
        rj[jc, C, :] = (a * dI[jc] / F) * dxc[:, None]
        rj[jc, C, C] = (a * dI[jc, C] / F) * dxc - (eps / dt) * dxc
        g[jc, C] = -((a * i[jc] / F) * dxc) - (fW[jc, C, C] * cW[jc, C] + dW[jc, C, C] * gW[jc, C]) + (fE[jc, C, C] * cE[jc, C] + dE[jc, C, C] * gE[jc, C])
        dW[jc, P1, P1] = -(1.0 - eps) * sig; dE[jc, P1, P1] = -(1.0 - eps) * sig
        rj[jc, P1, :] = -((a * dI[jc]) * dxc[:, None])
        g[jc, P1] = (a * i[jc]) * dxc - (fW[jc, P1, P1] * cW[jc, P1] + dW[jc, P1, P1] * gW[jc, P1]) + (fE[jc, P1, P1] * cE[jc, P1] + dE[jc, P1, P1] * gE[jc, P1])
        dW[jc, P2, C] = -(F * (z * self.D_cat + zn * self.D_an)); dE[jc, P2, C] = -(F * (z * self.D_cat + zn * self.D_an))
        fW[jc, P2, C] = -(F2 * k_cat * gW[jc, P2]); fE[jc, P2, C] = -(F2 * k_cat * gE[jc, P2])
        dW[jc, P2, P2] = -(F2 * k_cat * cW[jc, C]); dE[jc, P2, P2] = -(F2 * k_cat * cE[jc, C])
        rj[jc, P2, :] = (a * dI[jc]) * dxc[:, None]
        g[jc, P2] = -((a * i[jc]) * dxc) - (fW[jc, P2, P2] * cW[jc, P2] + dW[jc, P2, P2] * gW[jc, P2]) + (fE[jc, P2, P2] * cE[jc, P2] + dE[jc, P2, P2] * gE[jc, P2])   # D-11

        # current collector
        jn = slice(nj - 1, nj)
        dW[jn, C, C] = -self.D_cat
        fW[jn, C, C] = -(z * self.u_cat * F * gW[jn, P2])
        dW[jn, C, P2] = -(z * self.u_cat * F * cW[jn, C])
        g[jn, C] = 0.0 - dW[jn, C, C] * gW[jn, C] - fW[jn, C, C] * cW[jn, C]
        dW[jn, P1, P1] = -(1.0 - eps) * sig
        g[jn, P1] = I - dW[jn, P1, P1] * gW[jn, P1]
        dW[jn, P2, C] = -(F * (z * self.D_cat + zn * self.D_an))
        fW[jn, P2, C] = -(F2 * k_cat * gW[jn, P2])
        dW[jn, P2, P2] = -(F2 * k_cat * cW[jn, C])
        g[jn, P2] = 0.0 - dW[jn, P2, C] * gW[jn, C] - fW[jn, P2, C] * cW[jn, C]

        A, B, D = _blocks(rj, dW, dE, fW, fE, aW, aE, bW, bE)
        return A, B, D, g, i

    # -------------------------------------------------------------------- agglomerate scale
    def assemble_agglomerates(self, ca, phi2_electrode, dcs_dt):
        """Blocks for the agglomerates at every node in ca (shape (nl, nja, 4))."""
        p = self.p
        F, a, ea, sig = p.F, p.spec_a_agg, p.eps_agg, p.sigma
        z, zn = 1.0, -1.0
        F2 = F * F
        dt = self._dt
        nl, na = ca.shape[:2]
        aW, aE, bW, bE = (np.broadcast_to(v, (nl, na)) for v in self.faces_a)
        cW, cE, gW, gE = _face_values(ca, aW, aE, bW, bE)
        i, dI = self.kin.rate_and_fd(ca[..., C], ca[..., CS], ca[..., P1], ca[..., P2])
        # coupling current at the agglomerate surface
        c_spec_agg = dcs_dt * F / p.rho
        i_agg = c_spec_agg * self.vol_agg * (1 - ea) * p.rho / self.area_agg

        dW = np.zeros((nl, na, N, N)); dE = np.zeros_like(dW); fW = np.zeros_like(dW); fE = np.zeros_like(dW)
        rj = np.zeros_like(dW); g = np.zeros((nl, na, N))
        AW, AE, dV = self.A_W, self.A_E, self.dV
        k_agg = z ** 2 * self.u_agg + zn ** 2 * self.u_agg

        # center (r = 0): zero-gradient rows with the doubling sign (D-1)
        dE[:, 0, C, C] = 1.0
        g[:, 0, C] = 0.0 - dE[:, 0, C, C] * gE[:, 0, C]
        dE[:, 0, P1, P1] = -(1.0 - ea) * sig
        g[:, 0, P1] = 0.0 - dE[:, 0, P1, P1] * gE[:, 0, P1]
        dE[:, 0, P2, P2] = 1.0
        g[:, 0, P2] = 0.0 - dE[:, 0, P2, P2] * gE[:, 0, P2]
        rj[:, 0, CS, :] = -(a * dI[:, 0] / F)
        rj[:, 0, CS, CS] = -(a * dI[:, 0, CS] / F) - 1.0 * (1.0 - ea) / dt
        g[:, 0, CS] = +(a * i[:, 0] / F)

        # interior (spherical control volumes)
        ji = slice(1, na - 1)
        dW[:, ji, C, C] = AW[ji] * (-(ea * p.D_agg)); dE[:, ji, C, C] = AE[ji] * (-(ea * p.D_agg))
        rj[:, ji, C, :] = (a * dI[:, ji] / F) * dV[ji, None]
        rj[:, ji, C, C] = (a * dI[:, ji, C] / F) * dV[ji] - (ea / dt) * dV[ji]
        g[:, ji, C] = -((a * i[:, ji] / F) * dV[ji]) - (fW[:, ji, C, C] * cW[:, ji, C] + dW[:, ji, C, C] * gW[:, ji, C]) + (fE[:, ji, C, C] * cE[:, ji, C] + dE[:, ji, C, C] * gE[:, ji, C])
        rj[:, ji, CS, :] = -(a * dI[:, ji] / F)
        rj[:, ji, CS, CS] = -(a * dI[:, ji, CS] / F) - (1.0 - ea) / dt
        g[:, ji, CS] = +(a * i[:, ji] / F)
        dW[:, ji, P1, P1] = AW[ji] * (-(1.0 - ea) * sig); dE[:, ji, P1, P1] = AE[ji] * (-(1.0 - ea) * sig)
        rj[:, ji, P1, :] = -((a * dI[:, ji]) * dV[ji, None])
        g[:, ji, P1] = (a * i[:, ji]) * dV[ji] - (fW[:, ji, P1, P1] * cW[:, ji, P1] + dW[:, ji, P1, P1] * gW[:, ji, P1]) + (fE[:, ji, P1, P1] * cE[:, ji, P1] + dE[:, ji, P1, P1] * gE[:, ji, P1])
        dW[:, ji, P2, C] = AW[ji] * (-(ea * F * (z * p.D_agg + zn * p.D_agg))); dE[:, ji, P2, C] = AE[ji] * (-(ea * F * (z * p.D_agg + zn * p.D_agg)))
        fW[:, ji, P2, C] = AW[ji] * (-(ea * F2 * k_agg * gW[:, ji, P2])); fE[:, ji, P2, C] = AE[ji] * (-(ea * F2 * k_agg * gE[:, ji, P2]))
        dW[:, ji, P2, P2] = AW[ji] * (-(ea * F2 * k_agg * cW[:, ji, C])); dE[:, ji, P2, P2] = AE[ji] * (-(ea * F2 * k_agg * cE[:, ji, C]))
        fW[:, ji, P2, P2] = AW[ji] * 0.0; fE[:, ji, P2, P2] = AE[ji] * 0.0
        rj[:, ji, P2, :] = (a * dI[:, ji]) * dV[ji, None]
        g[:, ji, P2] = -((a * i[:, ji]) * dV[ji]) - (fW[:, ji, P2, P2] * cW[:, ji, P2] + dW[:, ji, P2, P2] * gW[:, ji, P2]) + (fE[:, ji, P2, P2] * cE[:, ji, P2] + dE[:, ji, P2, P2] * gE[:, ji, P2])

        # surface (r = R): c = c_bulk (D-17), imposed solid current, phi2 of the electrode
        js = na - 1
        rj[:, js, C, C] = 1.0
        g[:, js, C] = p.c_bulk - ca[:, js, C]
        rj[:, js, CS, :] = -(a * dI[:, js] / F)
        rj[:, js, CS, CS] = -(a * dI[:, js, CS] / F) - (1.0 - ea) / dt
        g[:, js, CS] = +(a * i[:, js] / F)
        dW[:, js, P1, P1] = -(1.0 - ea) * sig
        g[:, js, P1] = i_agg - dW[:, js, P1, P1] * gW[:, js, P1]
        rj[:, js, P2, P2] = 1.0
        g[:, js, P2] = phi2_electrode - ca[:, js, P2]

        A, B, D = _blocks(rj, dW, dE, fW, fE, aW, aE, bW, bE)
        return A, B, D, g


# ---------------------------------------------------------------------------------- output

HEADER = ("State", "Time", "Voltage", "mAhg", "Equivalence", "Solid_Conc", "current_density", "iloc",
          "Solution_Pot", "c0", "cs_edge", "cs", "eta_contact", "U", "eta_rxn", "i0", "current", "ramp_current")
UNITS = ("CDR", "hours", "Volts", "mAh/g", "LixNMC", "LixNMC", "A/cm2", "A/cm2", "Volts", "mol/cm3",
         "mol/cm3", "mol/cm3", "Volts", "Volts", "Volts", "A/cm2", "A/cm2", ",")


def _fixed(v):
    return f"{'NaN':>12}" if math.isnan(v) else f"{v:12.5f}"


def _sci(v):
    return f"{'NaN':>15}" if math.isnan(v) else f"{v:15.5E}"


def format_header():
    def row(cols):
        return (f"{cols[0][:5]:>5} " + " ".join(f"{c[:12]:>12}" for c in cols[1:3]) + " "
                + " ".join(f"{c[:15]:>15}" for c in cols[3:]))
    return row(HEADER) + "\n" + row(UNITS) + "\n"


def format_row(r):
    state, t, v, *rest = r
    return f"{state:>5} {_fixed(t)} {_fixed(v)} " + " ".join(_sci(x) for x in rest) + "\n"


@dataclass
class Result:
    rows: list = field(default_factory=list)
    exit_reason: str = ""
    steps: int = 0
    state: Optional[np.ndarray] = None
    agglomerates: Optional[np.ndarray] = None

    @property
    def array(self):
        return np.array([r[1:] for r in self.rows], dtype=float)

    def write(self, path):
        with open(path, "w") as fh:
            fh.write(format_header())
            for r in self.rows:
                fh.write(format_row(r))


# ---------------------------------------------------------------------------------- time loop

def _solve_faithful(A, B, D, G, backend):
    """The archival MATINV: legacy pivot, and a block is singular only when no nonzero pivot is left.

    The original then printed DETERM=0 and continued with an undefined result; here the
    update is NaN, which ends the run as the original's NaN did.
    """
    try:
        return bandsolver.solve(A, B, D, G, pivot="legacy", singular="exact", backend=backend)
    except (bandsolver.NonFiniteError, bandsolver.SingularBlockError):
        return np.full(G.shape, np.nan)


def run(p: AggParams, *, max_steps: Optional[int] = None, backend: str = "fortran", on_step=None) -> Result:
    """The original program's constant-current discharge (faithful mode)."""
    m = Model(p)
    kin = m.kin
    nj, s, na = p.nj, m.s, p.nja
    c = np.empty((nj, N))
    c[:, C] = p.c_bulk; c[:, P1] = p.phi1_init; c[:, P2] = 0.0; c[:, CS] = p.cs_init
    ca = np.empty((nj, na, N))
    ca[..., C] = p.c0_init; ca[..., CS] = p.cs_init; ca[..., P1] = p.phi1_init; ca[..., P2] = 0.0
    n_steps = p.n_steps
    t_max = p.t_max
    # tmax is an implicitly REAL parameter, so these divisions are single precision (D-4)
    f32t = np.float32(t_max)
    dt_nominal = float(f32t / np.float32(n_steps))
    write_every = float(f32t / np.float32(n_steps))
    dt = dt_nominal
    t = 0.0
    mAhg = 0.0
    last_write = 0                          # uninitialized in the original (D-19); 0 as on the HPC cluster
    state = "D"
    ramp = 1.0
    current = 0.0
    dc = np.zeros_like(c)
    to_electrons = 1.0 / p.rho * p.M
    res = Result()
    thr_depleted = f32(0.0001)
    thr_x = f32(0.55)
    lit36 = f32(3.6)
    n = n_steps if max_steps is None else max_steps

    def output():
        phi1, phi2, c0, cs = c[-1, P1], c[-1, P2], c[-1, C], ca[-1, -1, CS]
        u = float(kin.ocp(np.array([c0]), np.array([cs]))[0])
        iloc = float(kin.rate(np.array([c0]), np.array([cs]), np.array([phi1]), np.array([phi2]))[0])
        i0 = float(kin.exchange_current(np.array([c0]), np.array([cs]))[0])
        eta_contact = current * 0.0
        return (state, t / float(3600), phi1 - eta_contact, mAhg, mAhg * p.M * lit36 / p.F, cs * to_electrons,
                p.i_final, iloc, phi2, c0, cs, c[-1, CS], eta_contact, u, phi1 - phi2 - u, i0, current, ramp)

    for it in range(1, n + 1):
        if it == 1:
            res.rows.append(output())
        elif (t - last_write) >= write_every:
            res.rows.append(output())
            last_write = int(t - dt)
        elif it >= n_steps:
            res.rows.append(output())
        if c[-1, P1] >= 99.0 and state == "C":
            res.rows.append(output()); res.exit_reason = "end_of_charge"; break
        elif ca[s, -1, CS] * to_electrons >= thr_x:
            res.rows.append(output()); res.exit_reason = "x_limit"; break
        elif math.isnan(dc[0, C]):
            res.rows.append(output()); res.exit_reason = "nan"; break
        elif t >= 99.0 * 3600.0:
            res.rows.append(output()); res.exit_reason = "max_time"; break

        if state == "R":
            dt = dt * 1.0001
        elif c[-1, C] <= thr_depleted:
            dt = float(f32t / np.float32(np.float32(n_steps) * np.float32(10.0)))
        elif c[-1, C] <= f32(0.00001):
            dt = float(f32t / np.float32(np.float32(n_steps) * np.float32(100.0)))
        else:
            dt = dt_nominal
        t = t + dt

        # current ramp
        if ramp == 1.0:
            current = p.i_final / 50.0
            ramp = 2.0
        elif ramp == 2.0 and abs(current * 1.5) < abs(p.i_final):
            current = current * 1.5
        else:
            current = p.i_final
            ramp = 0.0

        # electrode scale (Coulomb counting at the collector node uses the final current)
        A, B, D, G, _ = m.assemble_electrode(c, ca, dt, current)
        mAhg = mAhg + 1000.0 * p.i_specific * dt / 3600.0
        dc = _solve_faithful(A, B, D, G, backend)
        c = c + dc

        # agglomerate scale, with the electrode solid uptake of this step (D-18)
        m._dt = dt
        Aa, Ba, Da, Ga = m.assemble_agglomerates(ca[s:], c[s:, P2], dc[s:, CS] / dt)
        for k in range(nj - s):
            dca = _solve_faithful(Aa[k], Ba[k], Da[k], Ga[k], backend)
            ca[s + k] = ca[s + k] + dca
        res.steps = it
        if on_step is not None:
            on_step(it, t, c, ca)
    else:
        res.exit_reason = res.exit_reason or "max_steps"
    res.state, res.agglomerates = c, ca
    return res
