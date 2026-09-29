# NMC111 model: equations and discretization

Sections 1–5 describe the two models **as implemented in the original research code**, which is the reference for the faithful mode. Sections 6 and 7 describe the corrected mode. Suspected defects are marked **[D-n]** and tracked in [deviations.md](deviations.md). Line references point into the private reference copies of the original programs, which are not distributed:

- `A:Lnnn`: the agglomerate model, conductivity-study template (`NMC111_agg.f95`)
- `C:Lnnn`: the uniform-particle ("electrode–crystal") model (`NMC111_Electrode_Crystal.f95`)

Both are one-dimensional half cells: lithium foil | separator | porous NMC111 cathode | current collector. They share the finite-volume layout and conventions of [lfp-model](https://github.com/jcbernard87/lfp-model), whose `docs/model.md` explains the grid, face interpolation and block assembly in detail. This document states what the NMC models do and where they differ.

| | Uniform-particle model (`particle_model = 'uniform'`) | Agglomerate model (`particle_model = 'agglomerate'`) |
|---|---|---|
| Original program | `NMC111_Electrode_Crystal.f95` | `NMC111_agg.f95` |
| Particle scale | none solved: each node's particles have one uniform concentration. The original's crystal scale is present but switched off (`crystal_scale = 0`, C:L104), as in LFP | porous spherical agglomerates (radius 5 µm) of uniform-concentration crystals (radius 200 nm), solved at every cathode node |
| Unknowns | electrode: c, Φ₁, Φ₂, c_s | electrode: c, Φ₁, Φ₂, c_s; each agglomerate: c, Φ₁, Φ₂, c_s along its radius |
| Nodes | 43 (separator 20 volumes, cathode 20) | 68–99 in the archived runs (template placeholder `NUMJ`); agglomerate 20–48 (`NUMC`) |

## 1. Shared electrode scale

**Grid.** The nodes are node 1 at the lithium foil, the separator volumes, the zero-volume interface node `SEP_NODE = 22`, the cathode volumes, and the zero-volume collector node `NJ`. Face values and gradients use distance weighting, exactly as in lfp-model.

**Rows per node** (the same row form as lfp-model: flux differences plus a local term; one linearized solve per time step **[D-7]**):

| Row | Separator interior | Cathode interior |
|---|---|---|
| 1 | ε_sep ∂c/∂t = −∂N₊/∂x | ε ∂c/∂t = −∂N₊/∂x + a i_n/F |
| 2 | ∂/∂x[(1−ε_sep)σ ∂Φ₁/∂x] = 0 (a solid potential in the separator has no physical meaning) | ∂i₁/∂x = −a i_n |
| 3 | ∂i₂/∂x = 0 | ∂i₂/∂x = +a i_n |
| 4 | ∂c_s/∂t = 0 | ε_AM ∂c_s/∂t = −a i_n/F |

N₊ is the cation flux and i₂ the ionic current of a binary 1:1 electrolyte (dilute-solution theory); i₁ = −(1−ε)σ ∂Φ₁/∂x.

**Boundary and interface rows** follow lfp-model:
- **Li foil:** N₊ = I/F, zero electronic-current row (Φ₁), Φ₂ = 0, c_s frozen.
- **Interface:** flux and current continuity; the solid balance holds at the zero-volume node.
- **Collector:** N₊ = 0, i₁ = I, i₂ = 0, solid balance.

Recurring defects:
- The Li-foil Φ₁ row doubles the old gradient instead of zeroing it **[D-1]**.
- The ionic-current residual in the separator, interface and cathode rows leaves out the diffusion current; only the Jacobian has it **[D-11]**.

**Differences between the two models at the electrode scale:**

| | Uniform-particle (C) | Agglomerate (A) |
|---|---|---|
| Transport coefficients | lfp-model form: flux = −ε(D_i/τ)…, cathode τ = ε^(−1/2), separator τ_sep = 4.8 | cathode: D_e = εD₀/τ with a fitted τ (private; see parameters.md), and D₊, D₋ from t₊ = 0.375 (the flux has no separate ε factor); separator: ε_sep·D₀/τ_sep with **the same D₀ for both ions** (t₊ = 0.5 in the separator) **[D-14]** |
| Porosity on the separator faces next to the boundary nodes | cathode ε on separator faces **[D-2]** | ε_sep throughout (D-2 does not occur) |
| Kinetic solid concentration c_s | the electrode unknown c_s | the **agglomerate outer-node** crystal concentration c_s(l, r = R) (A:L955); the electrode c_s is still solved but used only for coupling and output |
| Specific area a | 3ε_AM/R_crystal | 3 v_AM/R_agglomerate: agglomerates are treated as solid spheres (A:L155) **[D-15]** |

## 2. Kinetics and open-circuit potential

- **Butler–Volmer** as in lfp-model: `i_n = i₀[exp(α_a Fη/RT) − exp(−α_c Fη/RT)]`, α_a = α_c = 0.5, `i₀ = F k c^½ (c_s,max − c_s)^½ c_s^½`. The exchange current is rounded to single precision (implicit REAL `ex_1`/`ex_curr`) **[D-4]**. Derivatives are taken by finite differences with absolute steps of 10⁻⁶ **[D-6]**.
- **OCP** in both models: a Redlich–Kister expansion **plus a Nernst term in the electrolyte concentration**:

  U(θ, c) = U_ref + (RT/F)·ln[(c/c_bulk)·(1−θ)/θ] + Σₖ Aₖ[(2θ−1)^(k+1) − 2θk(1−θ)/(2θ−1)^(1−k)]

  The Nernst term in c comes on top of the c^½ dependence in i₀ **[D-16]**. The sum is accumulated in an implicitly single-precision variable (`Vint`) **[D-4]**.

  | | Uniform-particle (C:L145–L179) | Agglomerate (A:L208–L245) |
  |---|---|---|
  | Terms | 11, U_ref = 3.8686 | 12, U_ref = 3.8637 |
  | θ | x/0.55 | x/x_max, with x_max = M·Q_th·3600/F |
  | c_s,max | 0.55·ρ/M | ρ/M · M·Q_th·3600/F |

## 3. Uniform-particle ("electrode–crystal") model

This is the lfp-model program with NMC parameters (see [parameters.md](parameters.md)). All of lfp-model's `docs/model.md` §§1–9 apply, with these changes: the OCP above; θ and c_s,max above; 43 nodes (so a 20-volume cathode); τ_sep = 4.8; ε_AM = 0.4 with ε = 0.5 (so, unlike LFP, the volume fractions add up to 0.9); a fitted σ and k (private; see parameters.md); Φ₁,init = 4.1 V; 1C run of 36 000 steps of 1 s. The output file `Time_Voltage.txt` has the same seven columns, including the lithium-foil overpotential written as (RT/F)·ln(I/i₀) **[D-12]**. The run ends on a NaN when a particle overfills **[D-3]**.

## 4. Agglomerate model

**Agglomerate scale** (A:L1707–L2114). At every cathode node l (A:L790: nodes `SEP_NODE … NJ`, which includes the two zero-volume nodes), a porous sphere of radius R = 5 µm is solved on a radial finite-volume grid with zero-volume center and surface nodes. Spherical face areas 4πr² and volumes (4π/3)(r_E³ − r_W³) use the **single-precision** `PI` **[D-4]**.

| Row | Interior (volume ΔV) |
|---|---|
| 1 | ε_agg ∂c/∂t = ∇·(ε_agg D_agg ∇c) + a_agg i_n/F. Salt diffusion only: both ions share D_agg, and there is **no migration term** |
| 2 | ∇·((1−ε_agg)σ∇Φ₁) = a_agg i_n |
| 3 | ∇·(ε_agg F²(u₊+u₋) c ∇Φ₂) = −a_agg i_n. Migration only (the diffusion current vanishes because D₊ = D₋) |
| 4 | (1−ε_agg) ∂c_s/∂t = −a_agg i_n/F. Crystals are uniform-concentration, with a_agg = 3(1−ε_agg)/r_crystal and r_crystal = 200 nm |

- **Center (r = 0):** zero-gradient conditions on c, Φ₁ and Φ₂. All three rows use the sign that doubles the old gradient **[D-1]**. The solid balance is also solved at the center node.
- **Surface (r = R):**
  - c = **c_bulk**. It is not the local electrode concentration, so the agglomerates never see electrolyte depletion **[D-17]**; the coupling to the electrode value is commented out (A:L1718).
  - Solid current (1−ε_agg)σ ∂Φ₁/∂r = i_agg (A:L1938–L1939).
  - Φ₂ = Φ₂ of the electrode node.
  - The solid balance holds at the surface node.
- **Coupling current** (A:L1728–L1730): i_agg = (∂c_s,electrode/∂t)·(F/ρ)·(V_agg(1−ε_agg)ρ/A_agg), with `∂c_s/∂t = δc_s(l)/Δt` from the electrode solve of the same step. The electrode scale fixes each agglomerate's total uptake with its own kinetics; the agglomerate redistributes it internally with its own potentials. There are **two separate reaction descriptions** for one particle **[D-15]**.
- **Splitting.** Each step solves the electrode first, using the agglomerate state from the previous step, then every agglomerate. The two scales are not solved together **[D-18]**.

**Time stepping and control** (A:L604–L722, L558–L593):
- **Step size:** N_steps = ⌊3600·C_rate·Time_mod⌋ (e.g. about 347 000), Δt = t_max/N_steps with t_max = 72 000 s (about 0.2 s). Δt is reduced ×10 when the collector-node electrolyte concentration is at most 10⁻⁴ mol/cm³. The branch for 10⁻⁵ can never be reached.
- **Order within a step:** the time is advanced **before** the solve. The current then ramps up from I/50 by factors of 1.5 over the first steps (`current_ramp`), and the electrode and agglomerate solves follow.
- **Coulomb counting:** mAh/g is advanced inside the collector-node fill using the *final* current, so the ramp is not counted.
- **Output and exits:**
  - An output row is written almost every step. The timer `last_write_time` is implicitly an integer, and the write interval equals Δt **[D-4]**.
  - The run exits on:
    - a charge sentinel
    - x ≥ 0.55 at the surface crystal of the interface-node agglomerate (A:L662; a "cimax" exit, not a voltage cutoff) **[D-3]**
    - NaN
    - 99 h

**Output** (`Time_Voltage.txt`, A:L467–L489): 18 columns.
- State and time [h].
- Voltage = Φ₁(NJ) − I·R_contact. There is **no** lithium-foil overpotential; R_contact = 0.
- mAh/g, equivalents, surface crystal x at the collector node, final current density, local reaction rate, Φ₂, c and c_s at the collector.
- Electrode c_s, the contact drop, U, η, i₀, the ramped current, and the ramp flag.
- Two more files are written: `Time_Voltage_Position.txt` (every node, every row) and an empty `Crystal_Conc_Position.txt`.

## 6. Corrected agglomerate model (`mode = 'corrected'`)

The corrected agglomerate model is a double-porosity model with **one** reaction description. It fixes D-14, D-15, D-17 and D-18, and uses the corrections of the uniform model (D-1, D-3, D-4, D-6, D-7, D-11, D-12). Like the uniform model's corrected mode, it is solved in the variables of section 8, so the electrolyte can run out and the crystals can fill or empty smoothly.

**Electrode scale** (the macro-pores between agglomerates): unknowns c, Φ₁, Φ₂. Transport is the uniform model's corrected form: one binary electrolyte with t₊ = 0.375 everywhere (D-14), the full dilute-solution current (D-11), and τ = ε^−½ in the cathode. The electrode has no reaction term of its own. Each cathode volume holds agglomerates, which exchange salt, ionic current and electronic current with it through their surfaces:

| Row | Cathode volume Δx |
|---|---|
| c | ε ∂c/∂t Δx = −(N₊,E − N₊,W) − s_agg N₊,in Δx |
| Φ₁ | i₁,E − i₁,W = −s_agg i₁,in Δx |
| Φ₂ | i₂,E − i₂,W = −s_agg i₂,in Δx |

s_agg = 3 v_agg/R_agg is the agglomerate surface per electrode volume, with v_agg = v_AM/(1 − ε_agg). The `*,in` values are the fluxes into an agglomerate through r = R, taken from the agglomerate's own discretization, so the exchange conserves salt, lithium and charge exactly. The electrode's fourth unknown (c_s) is held fixed; the crystals carry the solid state.

**Agglomerate scale** (radius R, spherical finite volumes with zero-volume center and surface nodes): unknowns c, Φ₁, Φ₂, c_s.

| Row | Interior (volume ΔV, face areas 4πr²) |
|---|---|
| c | ε_agg ∂c/∂t ΔV = −Σ A N₊ + a_x i_n/F ΔV |
| Φ₁ | Σ A i₁ = −a_x i_n ΔV, with i₁ = −(1−ε_agg)σ_agg ∂Φ₁/∂r |
| Φ₂ | Σ A i₂ = +a_x i_n ΔV |
| c_s | (1−ε_agg) ∂c_s/∂t = −a_x i_n/F, at every node |

- The pores carry the electrode's electrolyte: N₊ and i₂ are the full dilute-solution expressions with D₊, D₋ scaled by ε_agg/τ_agg, τ_agg = ε_agg^−½.
- a_x = 3(1−ε_agg)/r_crystal.
- Butler–Volmer kinetics use the crystal concentration, the local pore electrolyte and the local potentials. U includes the electrolyte Nernst term (D-16), and the derivatives are analytic (D-6).
- **Center:** zero gradient of c, Φ₁ and Φ₂ (D-1 fixed).
- **Surface:** c, Φ₁ and Φ₂ equal the electrode's values at that node (D-17). The surface electrolyte therefore sees depletion, and the agglomerate solid is in electronic contact with the electrode solid.

**Solving: a condensed, fully coupled Newton step (D-18).** Each backward-Euler step is solved to convergence (scaled update ≤ 10⁻¹⁰). In every Newton iteration:

1. The agglomerates of all cathode volumes are assembled and stacked into one block-tridiagonal system, which is factored once. It is solved for the Newton update at fixed surface values, δx_a⁰, and for its response Z to the three surface values (c, Φ₁, Φ₂).
2. The agglomerate sources are linearized in the agglomerate unknowns next to the surface, and the response δx_a = δx_a⁰ + Z δx_e is substituted. This adds a 3 × 3 term to each cathode node's diagonal block and a term to its right-hand side, so the electrode system stays block-tridiagonal and is solved with BAND.
3. The agglomerate updates follow by back-substitution, and one step length, limited as in section 8, is applied to both scales.

The result is the Newton step of the full two-scale system: it converges quadratically, and after convergence the uncondensed residuals of both scales vanish (see validation.md).

**Time stepping and output** are shared with the uniform model: cc / cv / rest protocols with voltage cutoffs located to within 0.1 mV ([protocol.md](protocol.md)). The output adds `x_front` (the surface-crystal lithiation of the agglomerate next to the separator) and `c_collector` (the electrolyte concentration at the current collector).

## 7. Voltage reference in corrected mode (both models)

The cell voltage is measured against the lithium foil (0 V): V = Φ₁(collector) − U_Li − η_Li, with U_Li = (RT/F) ln(c(0)/c_Li,ref) and η_Li = (RT/(αF)) asinh(I/(2 i₀,Li)), α = 0.5. The solver fixes the gauge with Φ₂ = 0 at the foil face. The equations depend only on potential differences, so this is an exact shift of the foil-referenced solution. Imposing the foil potential as a boundary condition instead makes the first block singular at rest. The output column `Li_Nernst` is U_Li in mV. This fixes the missing foil Nernst term (D-16) and follows lfp-model v0.2.0.

## 8. Physical limits: log variables and exponential fitting (corrected mode, both models)

A discharge can exhaust the electrolyte locally (for example in the cores of agglomerates with slow internal transport), and particles or crystals can approach full or empty. Physically these limits are reached smoothly:
- where the salt runs out, the local reaction stops (i₀ ∝ c^½, and U contains ln c) while diffusion keeps feeding the region;
- a particle's OCP diverges as it fills, so it approaches full asymptotically;
- in both cases the cell voltage falls toward the cutoff.

Corrected mode is formulated so that these limits are reached the same way numerically, without clipping.

**Unknowns.** At every node the unknowns are u = ln(c/c_bulk), Φ₁, Φ₂ and the log-odds s = ln(θ/(1−θ)) of the particles or crystals, with θ = c_s/c_s,max. Then c = c_bulk·eᵘ > 0 and 0 < θ < 1 by construction. The storage terms stay in the conserved quantities, ε(c − c_old)/Δt and c_s,max(θ − θ_old)/Δt, so salt and lithium are conserved to round-off.

**Kinetics.** Both singular terms of the OCP become linear in the unknowns, and the exchange current is evaluated as a logarithm:

  U = U_ref + (RT/F)(u − s) + RK(θ),  ln i₀ = ln(F k c_bulk^α_a c_s,max^(α_a+α_c)) + α_a u − α_a ln(1 + eˢ) − α_c ln(1 + e⁻ˢ)

All derivatives stay bounded. As θ → 1 the cathodic current falls off like e⁻ˢ, so a full particle stops reacting smoothly. (This replaces lfp-model's regularization D-13, which is no longer needed.)

**Ion fluxes: exponential fitting (Scharfetter–Gummel).** For ion i with charge zᵢ across a face between nodes a and b a distance h apart, with Δ = zᵢF(Φ₂,b − Φ₂,a)/(RT):

  Nᵢ = (ε/τ)(Dᵢ/h)·[B(Δ)·c_a − B(−Δ)·c_b],  B(x) = x/(eˣ − 1)

This is exact for a constant flux in a constant field between the nodes, keeps concentrations positive, and stays accurate where migration dominates (the depleted regime). For |Δ| ≪ 1 it reduces to the centered scheme. The cation row uses N₊, and the current row uses i₂ = F(N₊ − N₋), so charge and salt are consistent by construction.

**Background conductivity.** Where the salt is exhausted no current can flow, and Φ₂ is physically undefined. The solvent's own ionic conductivity κ_bg (`kappa_bg`, default 10⁻⁸ S/cm, about 10⁻⁶ of the 1 M electrolyte's) is added as an ohmic current (ε/τ)κ_bg∇Φ₂ that carries no salt. It keeps Φ₂ defined there and changes nothing measurable elsewhere.

**Newton.**
- Each iteration's step is limited to |Δu| ≤ 1, |ΔΦ| ≤ 0.1 V and |Δs| ≤ 2; the exponentials make larger steps overshoot.
- Convergence is judged on the physical variables: |Δc|/c_bulk, |ΔΦ| and |Δθ| ≤ 10⁻¹⁰. Near a limit the log variables are ill-conditioned, and their round-off is physically irrelevant.
- A raw update above 10³ means the step has no solution (for example a current the remaining capacity cannot carry for the whole step). The driver then shortens the step. Near a full particle the voltage falls about 50 mV per e-fold of the remaining capacity, so the cutoff is approached with geometrically shrinking sub-steps (down to 10⁻¹⁰ s; docs/protocol.md).

## 9. Source lineage

| Copy | Role |
|---|---|
| Conductivity-study template `NMC111_agg.f95` | the reference for the agglomerate model. The run generator substitutes `NUMJ`, `NUMC`, `TIMEMOD`, `CRATE`, `POROSITY`, `THICKNESS` and `CONDUCTIVITY`; the archived snapshots also change `c0min` (10⁻³² → 10⁻²², unused) and comment out one output file |
| Quals / Binder Paper `NMC111_agg.f95` (byte-identical copies) | the same model as a fitting template with placeholders `MASS_LOADING`, `LENGTH`, `EFF_TORTUOSITY`, `ACTIVE_FRACTION`, `DIFF_AGLOMERATE`, `REACTIONK`; NJ = 75, agglomerate 33 nodes |
| `NMC111_agg_old.f95` | earlier version with the 11-term OCP fit (the one the uniform-particle model uses; this is the OCP/OSP difference in audit §8.20) and D_agg scaled 10⁻⁹ |
| `NMC111_agg_input.f95` (two copies) | command-line variant (12 arguments) of the earlier version |
| `NMC111_Electrode_Crystal.f95`, "questions" copy (Jan 2021) | the reference for the uniform-particle model; reproduces the archived output byte for byte |
| Blackbox Testing copy | command-line variant (thickness, C-rate, label) used for black-box optimization; the same model |
| Old Work copy | debugging copy: it has the corrected Li-foil Φ₁ sign but also an interface current row that adds the two face currents instead of subtracting them, and it diverges within 36 s. Not ported |
