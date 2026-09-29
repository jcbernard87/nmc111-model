"""Corrected mode of the uniform-particle model in the log variables (docs/model.md section 8).

The unknowns per node are (u, phi1, phi2, s): u = ln(c/c_bulk) and the particles' log-odds
s = ln(theta/(1-theta)) at the cathode nodes, with Scharfetter-Gummel ion fluxes
(nmc_model.logcore). A particle then fills or empties, and the electrolyte runs out, smoothly
and without clipping. Faithful mode keeps the original formulation (model.py).
"""
from __future__ import annotations

import math

import numpy as np

import bandsolver

from ..driver import SolverFailure
from ..logcore import DIVERGED, physical_update, N, P1, P2, S, U, Electrode, LogKinetics, Transport, bounded, converged, equilibrate, logit
from ..logcore import sigmoid
from . import kinetics
from .model import make_mesh
from .params import Params


class LogUniformModel:
    def __init__(self, p: Params):
        self.p = p
        self.mesh = make_mesh(p)
        self.cs_max = kinetics.cs_max(p)
        self.kin = LogKinetics(U_ref=kinetics._U_REF, ak=kinetics._AK, cs_max=self.cs_max, k=p.k_rxn,
                               alpha_a=p.alpha_a, alpha_c=p.alpha_c, c_bulk=p.c_bulk, R=p.R, T=p.T, F=p.F)
        self.tr = Transport(D0=p.D, t_plus=p.t_plus, c_bulk=p.c_bulk, kappa_bg=p.kappa_bg, R=p.R, T=p.T, F=p.F)
        self.el = Electrode(self.mesh, eps=p.eps, eps_sep=p.eps_sep, tau=p.tortuosity, tau_sep=p.tau_sep,
                            sigma=p.sigma, transport=self.tr, particles=(self.kin, p.spec_a, p.vf_AM))

    def initial_state(self) -> np.ndarray:
        p = self.p
        x = np.zeros((p.nj, N))
        x[:, U] = 0.0
        x[:, P1] = p.phi1_init
        x[:, P2] = p.phi2_init
        x[:, S] = logit(p.cs_init / self.cs_max)
        return x

    def conc(self, x):
        """Electrolyte concentration [mol/cm3] (nj,)."""
        return self.p.c_bulk * np.exp(x[:, U])

    def cs(self, x):
        """Particle concentration [mol/cm3] (nj,); meaningful at the cathode nodes s .. nj-1."""
        return self.cs_max * sigmoid(x[:, S])

    def newton_step(self, old: np.ndarray, dt: float, I: float, *, backend: str = "fortran",
                    start=None, history=None) -> np.ndarray:
        p = self.p
        x = (old if start is None else start).copy()
        prev = math.inf
        for _ in range(p.newton_max_iter):
            R, A, B, D = self.el.residual_and_blocks(x, old, dt, I)
            A, B, D, G = equilibrate(A, B, D, -R)
            try:
                dx = bandsolver.solve(A, B, D, G, backend=backend)
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(str(e)) from e
            lam = bounded((dx,))
            raw = float(np.max(np.abs(dx)))
            x = x + lam * dx
            upd = physical_update(x, dx)
            if history is not None:
                history.append(upd)
            if not math.isfinite(raw) or raw > DIVERGED:
                break
            if lam == 1.0 and converged(upd, prev, p.newton_tol):
                return x
            prev = upd
        raise SolverFailure("Newton did not converge")
