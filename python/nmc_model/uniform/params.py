"""Parameters of the uniform-particle NMC111 model.

Defaults are the values of the original research code `NMC111_Electrode_Crystal.f95`
(see docs/parameters.md).
`Params.faithful()` returns the values that program actually stored: several
constants were written as single-precision Fortran literals and are rounded to
float32 before being used in double-precision arithmetic (deviation D-4).
"""
from __future__ import annotations

from dataclasses import dataclass, fields, replace
from typing import Optional

import numpy as np

MODES = ("faithful", "corrected")
DT_MODES = ("fixed", "adaptive")
PARTICLE_MODELS = ("uniform", "crystal")   # the uniform-model family; the agglomerate model is AggParams
CRYSTAL_SHAPES = {"slab": 0, "cylinder": 1, "sphere": 2}


def check_dt_mode(p):
    """dt_mode is 'fixed' or 'adaptive', and adaptive only in corrected mode (both particle models); likewise
    restart files."""
    if (p.restart_read or p.restart_write) and p.mode != "corrected":
        raise ValueError("restart files need mode='corrected'")
    if p.dt_mode not in DT_MODES:
        raise ValueError(f"dt_mode must be one of {DT_MODES}, got {p.dt_mode!r}")
    if p.dt_mode == "adaptive" and p.mode != "corrected":
        raise ValueError("dt_mode='adaptive' needs mode='corrected'")


def f32(x: float) -> float:
    """Round to the nearest single-precision value, returned as a Python float."""
    return float(np.float32(x))


@dataclass(frozen=True)
class Params:
    # --- geometry and mesh ---
    L_cath_um: float = 24.0        # cathode thickness [um]; L_cath = L_cath_um * 1e-4 cm, as in the original
    L_sep: float = 25.0e-4         # separator thickness [cm]
    nj: int = 43                   # total nodes
    sep_node: int = 22             # separator/cathode interface node (1-based, as in the original)
    eps: float = 0.5               # cathode porosity
    eps_AM: float = 0.4            # active-material volume fraction (with eps = 0.5 the volumes add up to 0.9)
    eps_sep: float = 0.39          # separator porosity
    tau_sep: float = 4.8           # separator tortuosity
    bruggeman: float = -0.5        # cathode tortuosity = eps**bruggeman
    # --- electrolyte ---
    D: float = 2.0e-6              # salt diffusivity [cm2/s]
    t_plus: float = 0.25           # cation transference number
    c_bulk: float = 1.0e-3         # initial electrolyte concentration [mol/cm3]
    z_plus: float = 1.0
    z_minus: float = -1.0
    # --- active material ---
    sigma: float = 0.1             # electronic conductivity [S/cm] (generic; the original's value is private)
    M: float = 96.46               # molar mass [g/mol]
    rho: float = 4.6               # density [g/cm3]
    Q_th: float = 0.150            # theoretical capacity [Ah/g]
    R_p: float = 200.0e-7          # particle radius [cm] (crystal model: crystal radius, or a slab's half-thickness)
    particle_model: str = "uniform"  # 'uniform' particles, or 'crystal' (solid diffusion; corrected mode; #4)
    crystal_shape: str = "sphere"  # crystal model: 'sphere', 'cylinder' or 'slab'
    D_c: float = 5.162e-12         # crystal model: solid diffusivity [cm2/s] (the original's diff_c, C:L97)
    nj_crystal: int = 21           # crystal model: nodes across a crystal (center and surface included)
    k_rxn: float = 2.5e-6          # rate constant (generic, i0 ~ 0.1 mA/cm2 at theta = 0.5, 1 M)
    alpha_a: float = 0.5
    alpha_c: float = 0.5
    # --- lithium counter electrode (output only) ---
    k_Li: float = 1.0e-6
    c_Li_ref: float = 1.0e-3
    # --- constants ---
    R: float = 8.314
    T: float = 298.0
    F: float = 96485.0
    # --- operation ---
    C_rate: float = 1.0            # [1/h], positive = discharge
    phi1_init: float = 4.1         # [V]
    phi2_init: float = 0.0         # [V]
    cs_init: float = 1.0e-5        # [mol/cm3]
    t_max: float = 36000.0         # [s]
    n_steps: int = 36000
    V_min: float = 2.5             # [V] lower cutoff (corrected mode; the original has none, D-3)
    V_max: float = 4.2             # [V] upper cutoff (corrected mode)
    # --- numerics ---
    fd_step: float = 1.0e-6        # absolute step of the finite-difference reaction derivatives
    newton_tol: float = 1.0e-10    # corrected mode: scaled update tolerance per time step (D-7)
    newton_max_iter: int = 25
    dt_mode: str = "fixed"         # corrected mode: 'fixed' time steps, or 'adaptive' (docs/protocol.md)
    rtol: float = 1.0e-4           # adaptive: local error tolerance, relative ...
    atol: float = 1.0e-4           # ... and absolute, on the state (ln c, potentials [V], log-odds of theta)
    dt_min: float = 1.0e-6         # [s] adaptive: smallest step
    dt_max: float = 0.0            # [s] adaptive: largest step (0: write_interval)
    kappa_bg: float = 1.0e-8       # corrected mode: background (solvent) ionic conductivity [S/cm]
    # --- protocol and output (corrected mode; docs/protocol.md) ---
    steps: str = ""                # empty: classic discharge at C_rate
    cycles: int = 1
    write_interval: float = 18.0   # [s]
    restart_read: str = ""         # corrected mode: start from this restart file (#13); '' = the initial state
    restart_write: str = ""        # corrected mode: write the final state to this restart file; '' = none
    mode: str = "faithful"

    def __post_init__(self):
        if self.mode not in MODES:
            raise ValueError(f"mode must be 'faithful' or 'corrected', got {self.mode!r}")
        check_dt_mode(self)
        if self.particle_model not in PARTICLE_MODELS:
            raise ValueError(f"particle_model must be one of {PARTICLE_MODELS}, got {self.particle_model!r}")
        if self.crystal_shape not in CRYSTAL_SHAPES:
            raise ValueError(f"crystal_shape must be one of {tuple(CRYSTAL_SHAPES)}, got {self.crystal_shape!r}")
        if self.particle_model == "crystal":
            if self.mode != "corrected":
                raise ValueError("particle_model='crystal' needs mode='corrected'")
            if self.nj_crystal < 4:
                raise ValueError(f"nj_crystal must be at least 4, got {self.nj_crystal}")

    @classmethod
    def faithful(cls, *, sigma: float, k_exp: Optional[float] = None, k_rxn: Optional[float] = None,
                 **overrides) -> "Params":
        """Parameters exactly as stored by the original program (float32-rounded literals).

        The original's fitted values are required and are not part of this package: sigma [S/cm]
        and the rate constant, either as k_exp (k = 1e-8 * 10**k_exp, evaluated in single
        precision as the original did) or as the resulting value k_rxn.
        A `C_rate` override is also rounded to float32, because the original run generator
        wrote it into the source as a single-precision literal (for example `C_rate = 0.1`).
        """
        if (k_exp is None) == (k_rxn is None):
            raise TypeError("give exactly one of k_exp and k_rxn")
        if k_rxn is None:
            k_rxn = 1.0e-8 * f32(10.0 ** f32(k_exp))
        if "C_rate" in overrides:
            overrides["C_rate"] = f32(overrides["C_rate"])
        p = cls(
            R=f32(8.314), c_bulk=f32(0.001), Q_th=f32(0.150), M=f32(96.46), rho=f32(4.6),
            phi1_init=f32(4.1), eps_sep=f32(0.39), eps_AM=f32(0.4), tau_sep=f32(4.8),
            k_rxn=k_rxn, sigma=sigma, mode="faithful",
        )
        return replace(p, **overrides)

    def with_(self, **changes) -> "Params":
        return replace(self, **changes)

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    # --- derived quantities (evaluated in the same order as the original) ---
    @property
    def mol_vol(self) -> float:
        return self.rho / self.M

    @property
    def L_cath(self) -> float:
        """Cathode thickness [cm], computed as the original did (24 * 1.0d-4 is not 24.0e-4)."""
        return self.L_cath_um * 1.0e-4

    @property
    def vf_AM(self) -> float:
        """Active-material volume fraction of the electrode."""
        return self.eps_AM

    @property
    def spec_a(self) -> float:
        """Specific interfacial area a = 3 vf_AM / R_p [1/cm]."""
        return 3 * self.vf_AM / self.R_p

    @property
    def tortuosity(self) -> float:
        return self.eps ** self.bruggeman

    @property
    def mass_area(self) -> float:
        """Active-material loading [g/cm2]."""
        return self.L_cath * self.vf_AM * self.rho

    @property
    def i_1C(self) -> float:
        """1C current density [A/cm2]."""
        return self.Q_th * self.mass_area

    @property
    def i_specific(self) -> float:
        """Applied specific current [A/g]."""
        return self.Q_th * self.C_rate

    @property
    def i_app(self) -> float:
        """Applied current density [A/cm2]."""
        return self.i_specific * self.L_cath * self.vf_AM * self.rho

    @property
    def dt(self) -> float:
        return self.t_max / float(self.n_steps)
