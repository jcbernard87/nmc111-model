# Validation

Measured on macOS arm64 (gfortran 14.2, Apple clang 16, Python 3.13, numpy 2.5, bandsolver 0.1.2) on 2026-09-29. The public checks (`pytest tests`) recompute everything, and no results are stored in the repository.

## 1. Reproduction of the original programs (private)

The original programs and their archived outputs are kept privately, together with the fitted values the agglomerate model needs (see parameters.md). The agglomerate model's archived runs were produced on an x86-64 Linux cluster, and their outputs stop at the cluster's 20-minute job limit. The uniform-particle model's archived run was produced on the author's Mac.

- **Reference builds.** The original agglomerate snapshots were rebuilt with `gfortran -O0 -finit-local-zero`: zero-initialization reproduces the cluster's behaviour for an uninitialized output timer (deviation D-19). Over the archived prefixes, the rebuilt runs agree with the cluster outputs to within 0.01 mV in cell voltage (the last printed digit) for all ten runs.
- **Exact reproduction requirements.**
  - Reproducing the original programs' arithmetic bit for bit requires bandsolver to be built **without fused multiply-add contraction** (`-ffp-contract=off`). The original programs, built at `-O0`, never fuse. In the agglomerate model a last-bit difference at the first step grows about tenfold per step.
  - It also needs bandsolver's archival singular-block rule (`pivot="legacy", singular="exact"`, bandsolver ≥ 0.1.2). At step 290 of some runs a block is nearly singular; the original divides by its tiny pivot, where bandsolver's default rule would stop.

| Model | Implementation | Result |
|---|---|---|
| Uniform-particle | Python, Fortran, C++ | byte-identical to the archived output (211 lines) |
| Agglomerate | Fortran, C++ | byte-identical to the rebuilt original over the full length of all 10 archived snapshots |
| Agglomerate | Python | byte-identical to the rebuilt original for the first 1500 steps of all 10 snapshots (full runs take hours in Python); the internal state is bitwise identical over the first 14 steps |

## 2. Corrected mode (public, `tests/`)

Both models' corrected mode uses the log variables and exponential-fitting fluxes of docs/model.md section 8.

| Check | Result |
|---|---|
| Agglomerate model: lithium in the crystals vs the charge passed, 1C discharge to 3.0 V | equal to 10⁻¹² (test); measured 2 × 10⁻¹⁶ |
| Agglomerate model: salt in the macro-pores plus the agglomerate pores | constant to 10⁻¹³ (test) |
| Full (uncondensed) residuals of both scales after one condensed Newton step | ≤ 10⁻⁹ of their size at a perturbed state, every equation |
| Newton convergence from a perturbed iterate, including a state with drained agglomerate cores | quadratic (for example 2.1 × 10⁻³ → 1.3 × 10⁻⁶ → 6.4 × 10⁻¹³) |
| Charge balance of each agglomerate, i₁,in + i₂,in | zero to the Newton tolerance (≤ 10⁻⁷ of i₁,in) |
| Drained agglomerate cores (τ_agg = 22 400 at 1C, 2 240 at 3C) | the discharge runs smoothly to the 3.0 V cutoff, with the pore salt down to 2 × 10⁻¹² c_bulk; lithium conserved to round-off |
| Charge to nearly empty crystals (cycle with a 4.4 V CV step) | completes; θ < 10⁻³ at the end, with no clipping |
| Uniform model: assembled residual vs an independent node-by-node implementation (tests/uniform_reference_residual.py) | agree to 10⁻¹⁰ relative |
| Uniform model: `bandsolver.check_jacobian` at a mid-discharge state | max error < 10⁻⁴ (finite-difference noise) |
| Uniform model: salt and lithium over 1800 s at 1C; separator current | conserved to 10⁻¹²; i₂ = I to 10⁻⁹, no anion flux |
| Uniform model: discharge to the cutoff | reaches 2.5 V at 0.54979 electron equivalents (1C, 2C), the same capacity as the earlier c-variable formulation |
| Rest at equilibrium (no current, Φ₁ = U everywhere), both models | no change to 10⁻¹² |
| Time-step convergence: agglomerate (voltage after 600 s at 2C, Δt = 20 … 2.5 s); uniform (Li-face concentration at 8 s) | orders 0.92, 0.96; 0.96, 0.98 (backward Euler: 1) |
| Mesh convergence: agglomerate radial (4 … 32 volumes) and electrode (8 … 64 cathode volumes); uniform (26 … 201 nodes) | orders 1.99, 2.00; 2.00, 2.00; 2.47, 2.28 (approaching 2) |
| Fortran vs Python vs C++, agglomerate model: 1C, 2C, a full cycle, drained cores | identical output files (apart from round-off, about 10⁻¹⁰ mV, in the Li_Nernst column while c at the foil relaxes to c_ref) |
| Fortran vs C++, uniform model, 1C and 2C | identical files; Python agrees to rtol 10⁻⁵ |

A 1C discharge of the agglomerate model with 30 electrode nodes, 10 agglomerate nodes and 10 s steps takes about 0.4 s in Fortran or C++ and 1.5 s in Python.

## 3. Not covered yet

- Comparison with experimental data. This is done privately and is never published with the repository.
