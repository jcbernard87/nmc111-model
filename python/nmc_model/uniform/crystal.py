"""Crystal scale: solid diffusion in the crystals, coupled to the electrode (docs/model.md, crystal model; #4).

Corrected mode with particle_model = 'crystal'. Each interior cathode node carries one crystal
(slab, cylinder or sphere, k = 0, 1, 2), discretized by vertex-centred finite volumes in r. The unknown at every
crystal node is the log-odds s = ln(theta/(1-theta)), and the flux is Fick's law in theta with a constant D_c.
Each Newton iteration condenses the crystals into the electrode's diagonal blocks (the agglomerate model's
pattern). Ported from lfp-model's crystal model; the kinetics are NMC111's (Redlich-Kister OCP, no tails).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import bandsolver

from ..driver import SolverFailure, limit_reason
from ..logcore import (DIVERGED, N, P1, P2, S, U, Electrode, LogKinetics, Transport, bounded, converged,
                       equilibrate, logit, physical_update, sigmoid)
from . import kinetics
from .model import make_mesh
from .params import CRYSTAL_SHAPES, Params


class CrystalMesh:
    """Vertex-centred finite volumes in r: nodes at r_j = j h (center and surface included), each owning
    the control volume [r_j - h/2, r_j + h/2] within [0, R]. Every node stores lithium, so the surface
    can take up or give off lithium as it fills or empties, as in the continuum."""

    def __init__(self, R: float, nc: int, k: int):
        h = R / float(nc - 1)
        self.r = h * np.arange(nc)
        self.r[-1] = R
        rW = np.maximum(self.r - h / 2.0, 0.0)
        rE = np.minimum(self.r + h / 2.0, R)
        self.k, self.R, self.h = k, R, h
        self.V = (rE ** (k + 1) - rW ** (k + 1)) / (k + 1)   # node volumes (per unit geometric constant)
        self.A = (h * (np.arange(nc - 1) + 0.5)) ** k          # face areas between nodes f and f+1
        self.A_R = R ** k                                      # surface area
        self.Vtot = R ** (k + 1) / (k + 1)


def diffusion_rows(cm: CrystalMesh, D: float, cs_max: float, xc, xcold, dt):
    """Residual and tridiagonal dR/ds (sub An, diagonal B, super Dn) of the crystal rows, without the
    reaction. Every array is (nl, nc). Row j is storage plus the net outward diffusive flux of its
    control volume; the caller adds the flux out through the surface, A_R i/F, to the last row."""
    th, thold = sigmoid(xc), sigmoid(xcold)
    dth = th * sigmoid(-xc)
    g = D * cs_max * cm.A / cm.h                           # (nc-1,)
    Jf = -g * (th[:, 1:] - th[:, :-1])                     # outward flux through each internal face
    R = cs_max * cm.V * (th - thold) / dt
    R[:, :-1] += Jf
    R[:, 1:] -= Jf
    B = cs_max * cm.V * dth / dt
    B[:, :-1] += g * dth[:, :-1]
    B[:, 1:] += g * dth[:, 1:]
    Dn = np.zeros_like(xc); An = np.zeros_like(xc)
    Dn[:, :-1] = -g * dth[:, 1:]
    An[:, 1:] = -g * dth[:, :-1]
    return R, B, Dn, An


@dataclass
class CrystalState:
    x: np.ndarray     # electrode (nj, 4): u, phi1, phi2, s (s frozen)
    xc: np.ndarray    # crystals (nl, nc): log-odds at each crystal node

    def copy(self):
        return CrystalState(self.x.copy(), self.xc.copy())


class CrystalModel:
    """Corrected mode with a crystal at every interior cathode node (docs/model.md, crystal model)."""

    def __init__(self, p: Params):
        if p.particle_model != "crystal":
            raise ValueError("CrystalModel needs particle_model='crystal'")
        self.p = p
        self.mesh = make_mesh(p)
        m = self.mesh
        self.nodes = np.arange(m.s + 1, m.nj - 1)          # cathode volumes carry crystals
        self.nl = len(self.nodes)
        k = CRYSTAL_SHAPES[p.crystal_shape]
        self.cm = CrystalMesh(p.R_p, p.nj_crystal, k)
        self.a = (k + 1) * p.vf_AM / p.R_p                  # crystal surface per electrode volume
        self.cs_max = kinetics.cs_max(p)
        self.kin = LogKinetics(U_ref=kinetics._U_REF, ak=kinetics._AK, cs_max=self.cs_max, k=p.k_rxn,
                               alpha_a=p.alpha_a, alpha_c=p.alpha_c, c_bulk=p.c_bulk, R=p.R, T=p.T, F=p.F)
        self.tr = Transport(D0=p.D, t_plus=p.t_plus, c_bulk=p.c_bulk, kappa_bg=p.kappa_bg, R=p.R, T=p.T, F=p.F)
        self.el = Electrode(m, eps=p.eps, eps_sep=p.eps_sep, tau=p.tortuosity, tau_sep=p.tau_sep,
                            sigma=p.sigma, transport=self.tr)   # particles=None: the S column is frozen

    # ------------------------------------------------------------------ state helpers
    def initial_state(self) -> CrystalState:
        p = self.p
        x = np.zeros((p.nj, N))
        x[:, P1], x[:, P2] = p.phi1_init, p.phi2_init
        x[:, S] = logit(p.cs_init / self.cs_max)
        return CrystalState(x, np.full((self.nl, p.nj_crystal), logit(p.cs_init / self.cs_max)))

    def electrode(self, st: CrystalState) -> np.ndarray:
        return st.x

    def conc(self, st: CrystalState):
        """Electrolyte concentration [mol/cm3] (nj,)."""
        return self.p.c_bulk * np.exp(st.x[:, U])

    def cs(self, st: CrystalState):
        """Crystal concentrations [mol/cm3] (nl, nc)."""
        return self.cs_max * sigmoid(st.xc)

    def cs_mean(self, st: CrystalState):
        """Volume-mean crystal concentration at each crystal node's electrode node (nl,)."""
        return self.cs(st) @ (self.cm.V / self.cm.Vtot)

    def cs_surface(self, st: CrystalState):
        return self.cs(st)[:, -1]

    def surface_rate(self, st: CrystalState):
        """i_n (nl,) and d i_n / d(u, phi1, phi2, s_surface) (nl, 4) at each crystal surface."""
        y = st.x[self.nodes].copy()
        y[:, S] = st.xc[:, -1]
        return self.kin.rate(y)

    # ------------------------------------------------------------------ residuals
    def _pieces(self, st, old, dt, I):
        """Electrode residual/blocks with the reaction's (u, phi1, phi2) terms, crystal rows with the
        surface reaction, and the couplings: g (nl, 3), the electrode rows' dependence on the surface s,
        and Jce (nl, 3), the surface row's dependence on the electrode's (u, phi1, phi2)."""
        p, F = self.p, self.p.F
        Re, Ae, Be, De = self.el.residual_and_blocks(st.x, old.x, dt, I)
        i, di = self.surface_rate(st)
        w = self.a * self.mesh.dx[self.nodes]
        sgn = np.array([-1.0 / F, 1.0, -1.0])                  # rows U, P1, P2
        Re[self.nodes, :3] += (w * i)[:, None] * sgn
        Be[self.nodes, :3, :3] += (w[:, None, None] * sgn[None, :, None]) * di[:, None, :3]
        Rc, Bc, Dc, Ac = diffusion_rows(self.cm, p.D_c, self.cs_max, st.xc, old.xc, dt)
        AR = self.cm.A_R
        Rc[:, -1] += AR * i / F                                 # flux out through the surface
        Bc[:, -1] += AR * di[:, S] / F
        g = (w * di[:, S])[:, None] * sgn
        Jce = AR * di[:, :3] / F
        return Re, Ae, Be, De, Rc, Ac, Bc, Dc, g, Jce

    def residuals(self, st, old, dt, I):
        """Full (uncondensed) residuals, for tests: (electrode (nj, 4), crystals (nl, nc))."""
        pc = self._pieces(st, old, dt, I)
        return pc[0], pc[4]

    def jacobian_dense(self, st, old, dt, I):
        """The full uncondensed Jacobian (electrode x, then crystals xc, both row-major), for tests."""
        Re, Ae, Be, De, Rc, Ac, Bc, Dc, g, Jce = self._pieces(st, old, dt, I)
        nj, nl, nc = self.p.nj, self.nl, self.p.nj_crystal
        ne = nj * N
        J = np.zeros((ne + nl * nc, ne + nl * nc))
        for j in range(nj):
            r = slice(j * N, (j + 1) * N)
            J[r, r] = Be[j]
            if j > 0:
                J[r, (j - 1) * N:j * N] = Ae[j]
            if j < nj - 1:
                J[r, (j + 1) * N:(j + 2) * N] = De[j]
        for l in range(nl):
            for q in range(nc):
                row = ne + l * nc + q
                J[row, row] = Bc[l, q]
                if q > 0:
                    J[row, row - 1] = Ac[l, q]
                if q < nc - 1:
                    J[row, row + 1] = Dc[l, q]
            jn = self.nodes[l]
            J[jn * N:jn * N + 3, ne + l * nc + nc - 1] = g[l]
            J[ne + l * nc + nc - 1, jn * N:jn * N + 3] = Jce[l]
        return J

    # ------------------------------------------------------------------ Newton
    def newton_step(self, old: CrystalState, dt: float, I: float, *, backend: str = "fortran",
                    start=None, history=None) -> CrystalState:
        """One backward-Euler step from `old`, solved by the condensed Newton iteration."""
        p = self.p
        nl, nc = self.nl, p.nj_crystal
        st = (start or old).copy()
        prev = math.inf
        for _ in range(p.newton_max_iter):
            Re, Ae, Be, De, Rc, Ac, Bc, Dc, g, Jce = self._pieces(st, old, dt, I)
            # crystals: factor once (all stacked, block size 1); solve for the update and the three
            # responses of the surface to the electrode's u, phi1, phi2
            rhs = np.zeros((4, nl, nc))
            rhs[0] = -Rc
            rhs[1:, :, -1] = -Jce.T
            A1, B1, D1 = (M.reshape(-1, 1, 1) for M in (Ac, Bc, Dc))
            A1, B1, D1, *r = equilibrate(A1, B1, D1, *(q.reshape(-1, 1) for q in rhs))
            try:
                fac = bandsolver.factor(A1, B1, D1, backend=backend)
                sol = [fac.solve(q).reshape(nl, nc) for q in r]
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(f"crystal solve: {e}") from e
            dc0, Z = sol[0], np.stack(sol[1:], axis=-1)        # Z (nl, nc, 3)
            # electrode, with the crystals condensed into its diagonal blocks
            Be[self.nodes, :3, :3] += g[:, :, None] * Z[:, -1, None, :]
            G = -Re
            G[self.nodes, :3] -= g * dc0[:, -1, None]
            Ae, Be, De, G = equilibrate(Ae, Be, De, G)
            try:
                dxe = bandsolver.solve(Ae, Be, De, G, backend=backend)
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(f"electrode solve: {e}") from e
            dxc = dc0 + np.einsum("lnk,lk->ln", Z, dxe[self.nodes, :3])
            dxc4 = np.zeros((nl, nc, N))
            dxc4[..., S] = dxc
            lam = bounded((dxe, dxc4))
            # divergence is judged on the damped electrode update: near theta = 0 or 1 a linearized
            # log-odds update of the crystals is legitimately huge, and through the condensation it also
            # inflates the undamped electrode update; `bounded` scales both down
            raw = lam * float(np.max(np.abs(dxe)))
            st = CrystalState(st.x + lam * dxe, st.xc + lam * dxc)
            th = sigmoid(st.xc)
            upd = max(physical_update(st.x, dxe, frozen_s=True),
                      float(np.max(th * sigmoid(-st.xc) * np.abs(dxc))))
            if history is not None:
                history.append(upd)
            if not math.isfinite(raw) or raw > DIVERGED:
                break
            if lam == 1.0 and converged(upd, prev, p.newton_tol):
                return st
            prev = upd
        raise SolverFailure("Newton did not converge")

    def limit_reason(self, st: CrystalState):
        th = sigmoid(st.xc[:, -1])                              # the surfaces, where the reaction is
        return limit_reason(float(self.conc(st).min()), self.p.c_bulk, float(th.min()), float(th.max()))
