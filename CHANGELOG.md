# Changelog

## Unreleased

- **Restart files** (corrected mode, both models, all three languages; #13): `restart_write` writes the final state, time, capacity and current; `restart_read` starts a run from it. A continued run repeats the run it continues exactly (byte for byte in Fortran and C++); the implementations read each other's files.
- **CI** (#17): every supported Python version (3.9-3.13), coverage in the run summary.

## 0.3.0 (2026-10-14)

- **Fix:** a '/' inside a quoted value (an output path with directories) ended the namelist group in the C++ and Python readers, so the C++ program stopped with "cannot write" and Python lost the value. Both now treat quoted text as text, as the Fortran program does.
- **Crystal model** (`particle_model = 'crystal'`, corrected mode, all three languages; #4): the uniform model with solid-state diffusion in each cathode volume's crystal (sphere, cylinder or slab), ported from lfp-model: Fick's law in the log-odds, vertex-centred finite volumes in r, the crystals condensed into the electrode's Newton step. Default `D_c` = 5.162 × 10⁻¹² cm²/s, the original program's value. Four more output columns (surface and mean θ of the crystals next to the separator and the collector); `input/crystal.nml`. Fortran and C++ byte-identical, Python to the printed precision; validated as in lfp-model (exact constant-flux solution, second order in r and x, first order in t, conservation, the full Jacobian, the uniform limit).
- **Adaptive time stepping** (corrected mode, both models, all three languages; #8): `dt_mode = 'adaptive'` chooses each time step from backward Euler's local error estimate (linear predictor), against `rtol` and `atol`, within `dt_min` and `dt_max`; steps land on the output times. Agglomerate runs stay byte-identical across Python, Fortran and C++. The default stays `dt_mode = 'fixed'`, so existing inputs give the same output.

## 0.2.1 (2026-10-05)

- **Fix (D-20):** the OCP's Redlich-Kister sum evaluated 0/0 at θ = 1/2 exactly, so U was NaN there and a corrected run started at θ = 1/2 stopped at once (all three languages). Corrected mode evaluates the k = 0 term as 2θ − 1 alone; output elsewhere is bit-identical. Faithful mode keeps the original's arithmetic.
- **Input checks:** mesh sizes are validated (`sep_node` ≥ 3, `nj` − `sep_node` ≥ 3, `nja` ≥ 4); the Fortran program rejects unknown namelist groups (it skipped them silently) and names the offending entry; the C++ program rejects a name in the wrong group (it accepted any group); Python rejects a misspelled `mode` for the agglomerate model too; all three reject a group that appears twice (Fortran read the first, C++ merged them, Python kept the last). Messages agree across the three languages.
- **Tests and CI:** bad input rejected by all three implementations (#5), the OCP electrolyte term (#6), monotonic relaxation of an agglomerate rest (#7); ruff lint and a required Windows job (#9, #10); tests read text files as UTF-8.

## 0.2.0 (2026-10-04)

- **Driver fixes** (corrected mode, both models, all three languages): a cc discharge ends only at its `Vmin` and a charge only at its `Vmax`, as documented (both bounds were applied, so a step starting beyond the other bound stopped at once as a cutoff); when a time step cannot be solved, the exit row is the last converged sub-step, not the state at the start of the time step (the capacity of the partial step was lost). Faithful mode is unchanged. Found while porting the Zn/MnO₂ model.
- **Constant voltage** (corrected mode, both models, all three languages): a collapsed current bracket is accepted only within 10⁻⁶ V of the set voltage (#1); a CV step proceeds in sub-steps when no current holds the voltage for a whole time step (#3): a 4.2 V hold right after a 2C discharge stopped with `particles_full` and now runs to its current limit. The agglomerate Newton judges divergence on the damped electrode update (#2).

## 0.1.0 (2026-10-01)

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
