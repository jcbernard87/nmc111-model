"""Finite-volume discretization and block assembly (docs/model.md sections 1, 5-7).

Unknowns per node (0-based column): 0 = c, 1 = phi1, 2 = phi2, 3 = cs.
Nodes (0-based): 0 = Li-foil face, 1..s-1 = separator, s = separator/cathode
interface, s+1..nj-2 = cathode, nj-1 = current collector, with s = sep_node - 1.

Every row is written in control-volume form. For each node we collect face-flux
coefficients (d = coefficient of the face gradient, f = coefficient of the face
value), a local Jacobian `rj` and the right-hand side `g`, then convert them into
the BAND blocks A, B, D, G with A dc[j-1] + B dc[j] + D dc[j+1] = G.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import kinetics
from .params import Params

N = 4
C, P1, P2, CS = 0, 1, 2, 3


@dataclass(frozen=True)
class Mesh:
    nj: int
    s: int              # 0-based interface node
    dx: np.ndarray      # control-volume widths (0 at the boundary nodes)
    x: np.ndarray       # node positions [cm]
    aW: np.ndarray      # west-face interpolation weight
    aE: np.ndarray
    bW: np.ndarray      # west-face gradient factor
    bE: np.ndarray


def make_mesh(p: Params) -> Mesh:
    nj, s = p.nj, p.sep_node - 1
    h_sep = p.L_sep / float(p.sep_node - 2)
    h_cat = p.L_cath / float(nj - p.sep_node - 1)
    dx = np.zeros(nj)
    dx[1:s] = h_sep
    dx[s + 1:nj - 1] = h_cat
    x = np.zeros(nj)
    j = np.arange(nj)
    x[1:s] = h_sep * (j[1:s]) - h_sep / 2.0
    x[s] = p.L_sep
    x[s + 1:nj - 1] = p.L_sep + h_cat * (j[s + 1:nj - 1] - s) - h_cat / 2.0
    x[nj - 1] = p.L_cath + p.L_sep
    aW = np.zeros(nj); bW = np.zeros(nj); aE = np.zeros(nj); bE = np.zeros(nj)
    with np.errstate(invalid="ignore", divide="ignore"):
        aW[1:] = dx[:-1] / (dx[:-1] + dx[1:])
        bW[1:] = 2.0 / (dx[:-1] + dx[1:])
        aE[:-1] = dx[:-1] / (dx[1:] + dx[:-1])
        bE[:-1] = 2.0 / (dx[:-1] + dx[1:])
    return Mesh(nj, s, dx, x, aW, aE, bW, bE)


@dataclass(frozen=True)
class Transport:
    """Effective ion transport coefficients for one region (already divided by tortuosity)."""
    d_cat: float
    d_an: float
    u_cat: float
    u_an: float


def transport(p: Params, tau: float) -> Transport:
    t_an = 1.0 - p.t_plus
    d_cat = p.D * (1.0 + (t_an / p.t_plus)) / (2.0 * t_an / p.t_plus)
    d_an = d_cat * t_an / p.t_plus
    rt = p.R * p.T
    return Transport(d_cat / tau, d_an / tau, (d_cat / rt) / tau, (d_an / rt) / tau)


class Assembler:
    """Builds the linearized block system for one time step."""

    def __init__(self, p: Params):
        self.p = p
        self.mesh = make_mesh(p)
        self.tr_sep = transport(p, p.tau_sep)
        self.tr_cat = transport(p, p.tortuosity)
        faithful = p.mode == "faithful"
        # D-2: separator faces adjacent to the boundary nodes use the cathode porosity in faithful mode
        self.eps_sep_face = p.eps if faithful else p.eps_sep
        # D-1: sign of the Li-foil face solid-potential row
        self.phi1_row_sign = -1.0 if faithful else 1.0
        # D-11: include the diffusion current in the ionic-current residual
        self.full_current = not faithful

    # -- helpers -----------------------------------------------------------------
    def _cation(self, eps, tr, cface, gphi2):
        """Cation-flux coefficients: d11, f11, d13 (f13 = 0)."""
        p = self.p
        d11 = -(eps * tr.d_cat)
        f11 = -(eps * p.z_plus * tr.u_cat * p.F * gphi2)
        d13 = -(eps * p.z_plus * tr.u_cat * p.F * cface)
        return d11, f11, d13

    def _current(self, eps, tr, cface, gphi2):
        """Ionic-current coefficients: d31, f31, d33 (f33 = 0)."""
        p = self.p
        k = p.z_plus ** 2 * tr.u_cat + p.z_minus ** 2 * tr.u_an
        d31 = -(eps * p.F * (p.z_plus * tr.d_cat + p.z_minus * tr.d_an))
        f31 = -(eps * p.F ** 2 * k * gphi2)
        d33 = -(eps * p.F ** 2 * k * cface)
        return d31, f31, d33

    def time_terms(self, dt: float) -> np.ndarray:
        """Storage coefficients T (nj, 4): the time-derivative part of each row is T * (c - c_old)."""
        p, m = self.p, self.mesh
        s, nj = m.s, m.nj
        T = np.zeros((nj, N))
        T[0, CS] = -(p.vf_AM / dt)
        T[1:s, C] = -(p.eps_sep / dt * m.dx[1:s])
        T[1:s, CS] = -((1.0 - p.eps_sep) / dt)
        T[s + 1:nj - 1, C] = -((p.eps / dt) * m.dx[s + 1:nj - 1])
        T[s:, CS] = -(p.vf_AM / dt)
        return T

    # -- assembly ----------------------------------------------------------------
    def assemble(self, c: np.ndarray, dt: float, I: Optional[float] = None):
        """Return (A, B, D, G, i_rxn) for state c with shape (nj, 4) and applied current I [A/cm2]."""
        p, m = self.p, self.mesh
        nj, s = m.nj, m.s
        F, a = p.F, p.spec_a
        if I is None:
            I = p.i_app

        # face values and gradients
        cW = np.zeros((nj, N)); cE = np.zeros((nj, N)); gW = np.zeros((nj, N)); gE = np.zeros((nj, N))
        cW[1:] = m.aW[1:, None] * c[1:] + (1.0 - m.aW[1:, None]) * c[:-1]
        gW[1:] = m.bW[1:, None] * (c[1:] - c[:-1])
        cE[:-1] = m.aE[:-1, None] * c[1:] + (1.0 - m.aE[:-1, None]) * c[:-1]
        gE[:-1] = m.bE[:-1, None] * (c[1:] - c[:-1])

        rates = kinetics.reaction_derivatives if p.mode == "faithful" else kinetics.reaction_derivatives_analytic
        i, di_c, di_cs, di_p1, di_p2 = rates(p, c[:, C], c[:, CS], c[:, P1], c[:, P2])
        dI = np.stack([di_c, di_p1, di_p2, di_cs], axis=1)   # derivative w.r.t. each unknown column

        dW = np.zeros((nj, N, N)); dE = np.zeros((nj, N, N))
        fW = np.zeros((nj, N, N)); fE = np.zeros((nj, N, N))
        rj = np.zeros((nj, N, N)); g = np.zeros((nj, N))

        def west(rows, eps, tr, sl):
            """Fill west-face cation (row 0) and current (row 2) coefficients for node slice sl."""
            d11, f11, d13 = self._cation(eps, tr, cW[sl, C], gW[sl, P2])
            d31, f31, d33 = self._current(eps, tr, cW[sl, C], gW[sl, P2])
            dW[sl, C, C] = d11; fW[sl, C, C] = f11; dW[sl, C, P2] = d13
            dW[sl, P2, C] = d31; fW[sl, P2, C] = f31; dW[sl, P2, P2] = d33

        def east(rows, eps, tr, sl):
            d11, f11, d13 = self._cation(eps, tr, cE[sl, C], gE[sl, P2])
            d31, f31, d33 = self._current(eps, tr, cE[sl, C], gE[sl, P2])
            dE[sl, C, C] = d11; fE[sl, C, C] = f11; dE[sl, C, P2] = d13
            dE[sl, P2, C] = d31; fE[sl, P2, C] = f31; dE[sl, P2, P2] = d33

        def flux_W(sl, row, col_c):
            """Face flux value on the west face: f*c_face + d*grad (cation row uses c, solid row uses phi1)."""
            return fW[sl, row, col_c] * cW[sl, col_c] + dW[sl, row, col_c] * gW[sl, col_c]

        def flux_E(sl, row, col_c):
            return fE[sl, row, col_c] * cE[sl, col_c] + dE[sl, row, col_c] * gE[sl, col_c]

        def current_W(sl):
            v = fW[sl, P2, P2] * cW[sl, P2] + dW[sl, P2, P2] * gW[sl, P2]
            if self.full_current:
                v = v + dW[sl, P2, C] * gW[sl, C]
            return v

        def current_E(sl):
            v = fE[sl, P2, P2] * cE[sl, P2] + dE[sl, P2, P2] * gE[sl, P2]
            if self.full_current:
                v = v + dE[sl, P2, C] * gE[sl, C]
            return v

        eps, eps_sep, sig = p.eps, p.eps_sep, p.sigma
        tS, tC = self.tr_sep, self.tr_cat

        # ---- node 0: Li-foil face -------------------------------------------------
        j0 = slice(0, 1)
        d11, f11, d13 = self._cation(self.eps_sep_face, tS, cE[j0, C], gE[j0, P2])
        dE[j0, C, C] = d11; fE[j0, C, C] = f11; dE[j0, C, P2] = d13
        g[j0, C] = -I / F + (dE[j0, C, C] * gE[j0, C] + fE[j0, C, C] * cE[j0, C])
        rj[j0, CS, CS] = 0.0 - 1.0 * p.vf_AM / dt
        dE[j0, P1, P1] = -(1.0 - self.eps_sep_face) * sig
        g[j0, P1] = self.phi1_row_sign * dE[j0, P1, P1] * gE[j0, P1]
        rj[j0, P2, P2] = 1.0
        g[j0, P2] = 0.0 - c[j0, P2]

        # ---- separator interior ----------------------------------------------------
        js = slice(1, s)
        west(None, eps_sep, tS, js); east(None, eps_sep, tS, js)
        rj[js, C, C] = -(eps_sep / dt * m.dx[js])
        g[js, C] = 0.0 - flux_W(js, C, C) + flux_E(js, C, C)
        rj[js, CS, CS] = -(1.0 * (1.0 - eps_sep) / dt)
        dW[js, P1, P1] = -(1.0 - eps_sep) * sig; dE[js, P1, P1] = -(1.0 - eps_sep) * sig
        g[js, P1] = 0.0 - flux_W(js, P1, P1) + flux_E(js, P1, P1)
        g[js, P2] = 0.0 - current_W(js) + current_E(js)

        # ---- separator/cathode interface ------------------------------------------
        ji = slice(s, s + 1)
        west(None, self.eps_sep_face, tS, ji); east(None, eps, tC, ji)
        g[ji, C] = 0.0 - flux_W(ji, C, C) + flux_E(ji, C, C)
        dW[ji, P1, P1] = -(1.0 - self.eps_sep_face) * sig; dE[ji, P1, P1] = -(1.0 - eps) * sig
        g[ji, P1] = 0.0 - flux_W(ji, P1, P1) + flux_E(ji, P1, P1)
        g[ji, P2] = 0.0 - current_W(ji) + current_E(ji)

        # ---- solid concentration: interface, cathode interior and collector ---------
        jr = slice(s, nj)
        rj[jr, CS, :] = -(a * dI[jr] / F)
        rj[jr, CS, CS] = -(a * dI[jr, CS] / F) - 1.0 * p.vf_AM / dt
        g[jr, CS] = +(a * i[jr] / F)

        # ---- cathode interior -------------------------------------------------------
        jc = slice(s + 1, nj - 1)
        dxc = m.dx[jc]
        west(None, eps, tC, jc); east(None, eps, tC, jc)
        rj[jc, C, :] = (a * dI[jc] / F) * dxc[:, None]
        rj[jc, C, C] = (a * dI[jc, C] / F) * dxc - (eps / dt) * dxc
        g[jc, C] = -((a * i[jc] / F) * dxc) - flux_W(jc, C, C) + flux_E(jc, C, C)
        dW[jc, P1, P1] = -(1.0 - eps) * sig; dE[jc, P1, P1] = -(1.0 - eps) * sig
        rj[jc, P1, :] = -((a * dI[jc]) * dxc[:, None])
        g[jc, P1] = (a * i[jc]) * dxc - flux_W(jc, P1, P1) + flux_E(jc, P1, P1)
        rj[jc, P2, :] = (a * dI[jc]) * dxc[:, None]
        g[jc, P2] = -((a * i[jc]) * dxc) - current_W(jc) + current_E(jc)

        # ---- current collector --------------------------------------------------------
        jn = slice(nj - 1, nj)
        west(None, eps, tC, jn)
        g[jn, C] = 0.0 - dW[jn, C, C] * gW[jn, C] - fW[jn, C, C] * cW[jn, C]
        dW[jn, P1, P1] = -(1.0 - eps) * sig
        g[jn, P1] = I - dW[jn, P1, P1] * gW[jn, P1]
        g[jn, P2] = 0.0 - dW[jn, P2, C] * gW[jn, C] - fW[jn, P2, C] * cW[jn, C]
        # the collector row has no east face; clear the west-only helpers' unused rows
        # (west() also set cation/current coefficients, which is what this row needs)

        # ---- convert to BAND blocks -------------------------------------------------------
        aW, aE, bW, bE = (v[:, None, None] for v in (m.aW, m.aE, m.bW, m.bE))
        A = (1.0 - aW) * fW - bW * dW
        B = rj + bW * dW + aW * fW - (1.0 - aE) * fE + bE * dE
        D = -(aE * fE) - bE * dE
        # boundary nodes: only one face (their other-face coefficients are zero, but the
        # interpolation factors there are undefined, so drop those terms explicitly)
        B[0] = rj[0] - (1.0 - m.aE[0]) * fE[0] + m.bE[0] * dE[0]
        B[-1] = rj[-1] + m.bW[-1] * dW[-1] + m.aW[-1] * fW[-1]
        A[0] = 0.0
        D[-1] = 0.0
        return A, B, D, g, i
