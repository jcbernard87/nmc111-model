# Changelog

## 0.1.0 (2026-09-29)

First public version.

- Two particle models of an NMC111 half cell: uniform particles, and porous agglomerates of crystals.
- Fortran (`fortran/nmc.f90`), C++ (`cpp/nmc.cpp`) and Python (`nmc_model`) implementations that read one namelist input file and write the same output.
- Faithful mode reproduces the original research programs' output byte for byte, given the original's fitted parameter values, which are not distributed.
- Corrected mode fixes the defects listed in `docs/deviations.md`. For the agglomerate model:
  - a consistent double-porosity coupling with one reaction description
  - the local electrolyte at the agglomerate surfaces
  - one electrolyte (t₊ = 0.375) in the separator, the cathode and the agglomerate pores
  - a fully coupled, condensed Newton step for both scales
- Both models, corrected mode:
  - the cell voltage against the lithium foil
  - Newton iteration with exact kinetic derivatives
  - voltage cutoffs
  - the diffusion current in the ionic-current residual
  - a Butler–Volmer lithium counter electrode
- Physical limits handled smoothly in corrected mode: log-concentration and log-odds unknowns, Scharfetter–Gummel ion fluxes and a small background conductivity, so exhausted electrolyte (for example in agglomerate cores) and full or empty particles are reached without clipping or solver failures; runs that cannot continue stop with a physical exit reason (`electrolyte_depleted`, `particles_full`, `particles_empty`).
- Cycling protocols in corrected mode: constant current, constant voltage and rest steps, repeated `cycles` times (`docs/protocol.md`).
- Generic, documented public defaults (`docs/parameters.md`).
- Test suite (conservation of lithium and salt, independent residual, Jacobian, time and mesh convergence, protocols, cross-language agreement), CI, and three tutorial notebooks.
