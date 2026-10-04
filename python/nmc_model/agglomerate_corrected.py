"""Corrected mode of the agglomerate model (docs/model.md sections 6 and 8).

A double-porosity model with one reaction description:

- **Electrode scale** (macro-pores between agglomerates): the electrolyte concentration, phi1 and
  phi2, with one binary electrolyte everywhere. It has no reaction term of its own: the
  agglomerates at each cathode volume draw salt, ionic current and electronic current through
  their surfaces (D-15).
- **Agglomerate scale** (porous spheres of uniform-concentration crystals): the pore electrolyte,
  phi1, phi2 and the crystal lithiation along the radius, with Butler-Volmer kinetics on the
  crystals, and the electrode's values at the surface (D-17).
- **Fully coupled Newton** (D-18): each Newton iteration condenses the agglomerates onto the
  electrode's diagonal blocks. The agglomerates' linear response to their surface values is
  computed with one factorization of all agglomerate systems stacked into one block-tridiagonal
  system, so the electrode system stays block-tridiagonal (BAND).

The unknowns are the log variables of nmc_model.logcore: u = ln(c/c_bulk) and the crystals'
log-odds s = ln(theta/(1-theta)), with Scharfetter-Gummel fluxes, so electrolyte exhaustion and
full or empty crystals are reached smoothly (docs/model.md section 8). The electrode's fourth
unknown has no equation of its own here and is held fixed. The time driver (protocol, cutoffs,
CV) is nmc_model.driver.

A flux `q_in` is into an agglomerate through its surface.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

import bandsolver

from .agglomerate import AggParams, _AK, _U_REF
from .driver import SolverFailure, limit_reason, run_protocol
from .logcore import (DIVERGED, N, P1, P2, S, U, Electrode, LogKinetics, Transport, bounded, converged, equilibrate,
                      logit, physical_update, sigmoid)
from .uniform import kinetics as ukin
from .uniform.model import make_mesh
from .uniform.params import Params as _UParams
from .uniform.simulate import HEADER, UNITS, Result, format_row

# column names kept for readers of the state arrays
C, CS = U, S


# ---------------------------------------------------------------------------------- state

@dataclass
class AggState:
    c: np.ndarray        # electrode (nj, 4): u, phi1, phi2, (fixed)
    ca: np.ndarray       # agglomerates at the cathode volumes (nl, nja, 4): u, phi1, phi2, s

    def copy(self):
        return AggState(self.c.copy(), self.ca.copy())


# ---------------------------------------------------------------------------------- model

def electrode_params(p: AggParams) -> _UParams:
    """The electrode's mesh and lithium-foil parameters in the uniform model's form."""
    return _UParams(
        L_cath_um=p.thickness_um, L_sep=p.L_sep, nj=p.nj, sep_node=p.sep_node, eps=p.eps, eps_AM=0.0,
        eps_sep=p.eps_sep, tau_sep=p.tau_sep, bruggeman=p.bruggeman, D=p.D0, t_plus=p.t_plus,
        c_bulk=p.c_bulk, sigma=p.sigma, M=p.M, rho=p.rho, Q_th=p.Q_th, k_Li=p.k_Li, c_Li_ref=p.c_Li_ref,
        R=p.R, T=p.T, F=p.F, phi1_init=p.phi1_init, cs_init=p.cs_init, mode="corrected",
    )


class CorrectedModel:
    """Residuals, Jacobians and the condensed Newton step of the corrected agglomerate model."""

    def __init__(self, p: AggParams):
        if p.mode == "faithful":
            raise ValueError("CorrectedModel needs mode='corrected'")
        self.p = p
        self.up = electrode_params(p)
        self.mesh = make_mesh(self.up)
        m = self.mesh
        self.nodes = np.arange(m.s + 1, m.nj - 1)           # cathode volumes carry agglomerates
        self.nl = len(self.nodes)
        x_max = p.M * p.Q_th * 3600 / p.F
        self.cs_max = p.mol_vol * x_max
        self.kin = LogKinetics(U_ref=_U_REF, ak=_AK, cs_max=self.cs_max, k=p.rxn_k, alpha_a=p.alpha_a,
                               alpha_c=p.alpha_c, c_bulk=p.c_bulk, R=p.R, T=p.T, F=p.F)
        self.tr = Transport(D0=p.D0, t_plus=p.t_plus, c_bulk=p.c_bulk, kappa_bg=p.kappa_bg, R=p.R, T=p.T, F=p.F)
        self.el = Electrode(m, eps=p.eps, eps_sep=p.eps_sep, tau=p.eps ** p.bruggeman, tau_sep=p.tau_sep,
                            sigma=p.sigma, transport=self.tr)
        # radial grid (zero-volume center and surface nodes)
        na = p.nja
        h = p.R_agg / float(na - 2)
        j = np.arange(na)
        xa = np.zeros(na)
        xa[1:na - 1] = h * j[1:na - 1] - h / 2.0
        xa[na - 1] = p.R_agg
        dxa = np.zeros(na)
        dxa[1:na - 1] = h
        self.xa, self.dxa = xa, dxa
        r_face = xa[:-1] + dxa[:-1] / 2.0                   # face k is on the east side of node k
        r_face[0] = 0.0
        self.a_area = 4.0 * math.pi * r_face ** 2
        a_h = (dxa[:-1] + dxa[1:]) / 2.0
        r_W, r_E = xa - dxa / 2.0, xa + dxa / 2.0
        r_W[0] = r_E[0] = 0.0
        r_W[-1] = r_E[-1] = p.R_agg
        self.dV = 4.0 * math.pi / 3.0 * (r_E ** 3 - r_W ** 3)
        tau = p.tau_agg if (p.tau_agg is not None and p.tau_agg > 0) else p.eps_agg ** -0.5
        sig = p.sigma_agg if (p.sigma_agg is not None and p.sigma_agg > 0) else p.sigma
        self.a_g = p.eps_agg / tau / a_h                     # per agglomerate face
        self.a_gs = (1.0 - p.eps_agg) * sig / a_h
        self.a_x = 3.0 * (1.0 - p.eps_agg) / p.R_xtal          # crystal area per agglomerate volume
        self.v_agg = p.v_AM / (1.0 - p.eps_agg)                 # agglomerate volume fraction
        self.s_agg = 3.0 * self.v_agg / p.R_agg                 # agglomerate surface per electrode volume
        if p.eps + self.v_agg > 1.0 + 1e-12:
            raise ValueError(f"porosity {p.eps} + agglomerate volume fraction {self.v_agg:.3f} exceed 1")

    # ------------------------------------------------------------------ state helpers
    def initial_state(self) -> AggState:
        p = self.p
        c = np.zeros((p.nj, N))
        c[:, P1] = p.phi1_init
        ca = np.zeros((self.nl, p.nja, N))
        ca[..., U] = math.log(p.c0_init / p.c_bulk)
        ca[..., P1] = p.phi1_init
        ca[..., S] = logit(p.cs_init / self.cs_max)
        return AggState(c, ca)

    def conc(self, st: AggState):
        """Electrolyte concentrations [mol/cm3]: (electrode (nj,), agglomerate pores (nl, nja))."""
        return self.p.c_bulk * np.exp(st.c[:, U]), self.p.c_bulk * np.exp(st.ca[..., U])

    def cs(self, st: AggState) -> np.ndarray:
        """Crystal concentrations [mol/cm3] (nl, nja)."""
        return self.cs_max * sigmoid(st.ca[..., S])

    def mean_cs(self, st: AggState) -> np.ndarray:
        """Volume-average crystal concentration of each agglomerate (nl,)."""
        return self.cs(st) @ (self.dV / self.dV.sum())

    # ------------------------------------------------------------------ agglomerate scale
    def agglomerates(self, xa, xaold, dt, xe):
        """Residual and blocks of every agglomerate (nl, na, ...); surface values xe (nl, 4)."""
        p = self.p
        nl, na = xa.shape[:2]
        ea, F = p.eps_agg, p.F
        Fv, dFa, dFb = self.tr.fluxes(xa[:, :-1], xa[:, 1:], self.a_g, self.a_gs)
        i, di = self.kin.rate(xa)
        Ar = self.a_area
        R = np.zeros_like(xa)
        A = np.zeros((nl, na, N, N)); B = np.zeros_like(A); D = np.zeros_like(A)
        # center: zero gradient
        for col in (U, P1, P2):
            R[:, 0, col] = xa[:, 1, col] - xa[:, 0, col]
            B[:, 0, col, col] = -1.0
            D[:, 0, col, col] = 1.0
        # interior volumes
        ji = np.arange(1, na - 1)
        V = self.dV[ji]
        for f, row in enumerate((U, P1, P2)):
            R[:, ji, row] = Ar[ji] * Fv[:, ji, f] - Ar[ji - 1] * Fv[:, ji - 1, f]
            B[:, ji, row] += Ar[ji, None] * dFa[:, ji, f] - Ar[ji - 1, None] * dFb[:, ji - 1, f]
            D[:, ji, row] += Ar[ji, None] * dFb[:, ji, f]
            A[:, ji, row] += -Ar[ji - 1, None] * dFa[:, ji - 1, f]
        c = p.c_bulk * np.exp(xa[..., U])
        cold = p.c_bulk * np.exp(xaold[..., U])
        R[:, ji, U] += ea * V * (c[:, ji] - cold[:, ji]) / dt - self.a_x * i[:, ji] * V / F
        B[:, ji, U, U] += ea * V * c[:, ji] / dt
        B[:, ji, U, :] -= self.a_x * di[:, ji] * V[None, :, None] / F
        R[:, ji, P1] += self.a_x * i[:, ji] * V
        B[:, ji, P1, :] += self.a_x * di[:, ji] * V[None, :, None]
        R[:, ji, P2] -= self.a_x * i[:, ji] * V
        B[:, ji, P2, :] -= self.a_x * di[:, ji] * V[None, :, None]
        # surface: the electrode's u, phi1, phi2
        for col in (U, P1, P2):
            R[:, -1, col] = xa[:, -1, col] - xe[:, col]
            B[:, -1, col, col] = 1.0
        # crystals, every node: (1 - eps_agg) dcs/dt = -a i/F
        th, thold = sigmoid(xa[..., S]), sigmoid(xaold[..., S])
        R[..., S] = (1.0 - ea) * self.cs_max * (th - thold) / dt + self.a_x * i / F
        B[..., S, :] = self.a_x * di / F
        B[..., S, S] += (1.0 - ea) * self.cs_max * th * sigmoid(-xa[..., S]) / dt
        return R, A, B, D

    def surface_flux(self, xa):
        """q_in (nl, 3) = (N+, i1, i2) into each agglomerate through r = R, and dq/dx at its last two
        nodes (nl, 3, 2, 4)."""
        Fv, dFa, dFb = self.tr.fluxes(xa[:, -2], xa[:, -1], self.a_g[-1], self.a_gs[-1])
        return -Fv, -np.stack([dFa, dFb], axis=2)

    # ------------------------------------------------------------------ Newton
    def residuals(self, st: AggState, old: AggState, dt, I):
        """Full (uncondensed) residuals, for tests: (electrode (nj, 4), agglomerates (nl, na, 4))."""
        Re = self.el.residual_and_blocks(st.c, old.c, dt, I)[0]
        q, _ = self.surface_flux(st.ca)
        Re[self.nodes, :3] += (self.s_agg * self.mesh.dx[self.nodes])[:, None] * q
        Ra = self.agglomerates(st.ca, old.ca, dt, st.c[self.nodes])[0]
        return Re, Ra

    def newton_step(self, old: AggState, dt: float, I: float, *, backend: str = "cpp",
                    start: Optional[AggState] = None, history: Optional[list] = None) -> AggState:
        """One backward-Euler step from `old`, solved by the condensed Newton iteration.

        `start` is the first iterate (default: `old`); `history` collects each iteration's update.
        """
        p = self.p
        nl, na = self.nl, p.nja
        st = (start or old).copy()
        w = self.s_agg * self.mesh.dx[self.nodes]
        prev = math.inf
        for _ in range(p.newton_max_iter):
            # agglomerates: factor once (all stacked); solve for the update and the three surface responses
            Ra, Aa, Ba, Da = self.agglomerates(st.ca, old.ca, dt, st.c[self.nodes])
            rhs = np.zeros((4, nl, na, N))
            rhs[0] = -Ra
            for k, col in enumerate((U, P1, P2), start=1):
                rhs[k, :, na - 1, col] = 1.0
            Aa, Ba, Da, *r = equilibrate(Aa, Ba, Da, *rhs)
            try:
                fac = bandsolver.factor(Aa.reshape(-1, N, N), Ba.reshape(-1, N, N), Da.reshape(-1, N, N),
                                        backend=backend)
                sol = [fac.solve(rk.reshape(-1, N)).reshape(nl, na, N) for rk in r]
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(f"agglomerate solve: {e}") from e
            dxa0 = sol[0]
            Z = np.stack(sol[1:], axis=-1)                     # (nl, na, 4, 3): response to the electrode's u, phi1, phi2
            # electrode, with the condensed agglomerate sources: R_e + w q, linearized in the agglomerates
            Re, Ae, Be, De = self.el.residual_and_blocks(st.c, old.c, dt, I)
            q, Jq = self.surface_flux(st.ca)
            Re[self.nodes, :3] += w[:, None] * q
            Jf = Jq.reshape(nl, 3, 2 * N)
            Be[self.nodes, :3, :3] += w[:, None, None] * np.einsum("lrk,lkc->lrc", Jf,
                                                                    Z[:, na - 2:].reshape(nl, 2 * N, 3))
            G = -Re
            G[self.nodes, :3] -= w[:, None] * np.einsum("lrk,lk->lr", Jf, dxa0[:, na - 2:].reshape(nl, 2 * N))
            Ae, Be, De, G = equilibrate(Ae, Be, De, G)
            try:
                dce = bandsolver.solve(Ae, Be, De, G, backend=backend)
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(f"electrode solve: {e}") from e
            dca = dxa0 + np.einsum("lnkc,lc->lnk", Z, dce[self.nodes, :3])
            lam = bounded((dce, dca))
            # divergence is judged on the damped electrode update: near theta = 0 or 1 a linearized log-odds
            # update of the agglomerates is legitimately huge, and through the condensation it also inflates
            # the undamped electrode update; `bounded` scales both down
            raw = lam * float(np.max(np.abs(dce)))
            st = AggState(st.c + lam * dce, st.ca + lam * dca)
            upd = max(physical_update(st.c, dce, frozen_s=True), physical_update(st.ca, dca))
            if history is not None:
                history.append(upd)
            if not math.isfinite(raw) or raw > DIVERGED:
                break
            if lam == 1.0 and converged(upd, prev, p.newton_tol):
                return st
            prev = upd
        raise SolverFailure("Newton did not converge")

    # ------------------------------------------------------------------ outputs
    def foil(self, st: AggState, I: float):
        u_li, eta_li = ukin.li_foil(self.up, float(self.p.c_bulk * math.exp(st.c[0, U])), I)
        return float(u_li), float(eta_li)

    def voltage(self, st: AggState, I: float) -> float:
        """Cell voltage against the lithium foil (0 V): phi1(collector) - U_Li - eta_Li."""
        u_li, eta_li = self.foil(st, I)
        return float(st.c[-1, P1]) - u_li - eta_li


# ---------------------------------------------------------------------------------- driver glue

HEADER_EXTRA = ("Current", "Step", "Li_Nernst", "x_front", "c_collector")
UNITS_EXTRA = ("mA/cm2", "#", "mV", "LixNMC", "mol/cm3")


class AggStepper:
    """The corrected agglomerate model as a stepper for nmc_model.driver."""

    def __init__(self, p: AggParams, *, backend: str = "cpp"):
        self.p = p
        self.model = CorrectedModel(p)
        self.backend = backend

    def initial_state(self):
        return self.model.initial_state()

    def newton_step(self, st, h, I):
        return self.model.newton_step(st, h, I, backend=self.backend)

    def voltage(self, st, I):
        return self.model.voltage(st, I)

    def row(self, t, st, mAhg, I, k):
        m, p = self.model, self.p
        state = "D" if I > 0 else ("C" if I < 0 else "R")
        ce, _ = m.conc(st)
        u_li, eta_li = m.foil(st, I)
        i0_li = float(ukin.li_exchange_current(m.up, ce[0]))
        x_front = float(m.cs(st)[0, -1] / p.mol_vol)
        return (state, t / 3600.0, m.voltage(st, I), mAhg * p.M * 3.6 / p.F, -eta_li * 1.0e3, i0_li * 1.0e3,
                float(ce[0]), I * 1.0e3, k, u_li * 1.0e3, x_front, float(ce[-1]))

    @staticmethod
    def finite(st):
        return bool(np.all(np.isfinite(st.c)) and np.all(np.isfinite(st.ca)))

    def limit_reason(self, st):
        m = self.model
        ce, ca = m.conc(st)
        th = m.cs(st) / m.cs_max
        return limit_reason(min(float(ce.min()), float(ca.min())), self.p.c_bulk, float(th.min()), float(th.max()))


class AggResult(Result):
    def write(self, path) -> None:
        def row(cols):
            return (f"{cols[0][:5]:>5} " + " ".join(f"{c[:12]:>12}" for c in cols[1:3]) + " "
                    + " ".join(f"{c[:15]:>15}" for c in cols[3:]))
        with open(path, "w") as fh:
            fh.write(row(HEADER + HEADER_EXTRA) + "\n" + row(UNITS + UNITS_EXTRA) + "\n")
            for r in self.rows:
                fh.write(format_row(r[0], r[1:]))


def run_corrected(p: AggParams, *, max_steps: Optional[int] = None, backend: str = "cpp") -> AggResult:
    return run_protocol(AggStepper(p, backend=backend), max_steps=max_steps, result=AggResult())
