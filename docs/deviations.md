# Deviations from the original research code

Faithful mode reproduces the original programs, including every item below. Corrected mode applies the fixes. Where a defect is the same as one found in [lfp-model](https://github.com/jcbernard87/lfp-model/blob/main/docs/deviations.md), it keeps the same number. U = uniform-particle model, A = agglomerate model.

**Corrected mode (2026-09-29).** The fixes below are in all three implementations (docs/model.md §§6–7):

| ID | Fix in corrected mode |
|---|---|
| D-1 | zero-gradient rows impose (∂/∂x)_new = 0 (Li foil Φ₁; agglomerate center c, Φ₁, Φ₂) |
| D-2 | every separator face uses ε_sep |
| D-3 | voltage cutoffs `V_min`/`V_max`, located to within 0.1 mV; bound-preserving Newton steps |
| D-4 | intended decimal values and double precision throughout; output every `write_interval` seconds |
| D-6 | analytic reaction derivatives, including dU/dθ and dU/dc |
| D-7 | Newton's method to convergence at every step (scaled update ≤ 10⁻¹⁰, or stagnated at the round-off floor below 10⁻⁷ near a full or empty particle), with sub-step halving on failure |
| D-10 | the out-of-bounds read is not ported |
| D-11 | the ionic-current residual includes the diffusion current (electrode and agglomerate) |
| D-12 | the Li counter electrode is a symmetric Butler–Volmer interface, η = (RT/(αF))·asinh(I/(2i₀)) |
| D-13 | not needed: lfp-model regularizes c_s^α near an empty or full particle; here the log-odds unknowns and the logarithmic exchange current keep every derivative bounded (docs/model.md section 8) |
| D-14 | one electrolyte (t₊ = 0.375) in the separator, the cathode and the agglomerate pores |
| D-15 | one reaction description: the kinetics act only inside the agglomerates, and the electrode's sources are the agglomerate surface fluxes |
| D-16 | the OCP keeps its Nernst term; the voltage is measured against the lithium foil and includes the foil's Nernst term |
| D-17 | the agglomerate surface takes the local electrode c (and Φ₁, Φ₂) |
| D-18 | both scales solved together by a condensed Newton step |
| D-19 | the output timer is a double initialized to 0 |
| D-20 | the Redlich-Kister sum's k = 0 term is evaluated as 2θ − 1 alone |

Status values: **candidate** (suspected from the source), **confirmed** (demonstrated by a test or run), **fixed**, **kept** (reviewed and left as is, with the reason).

| ID | Model | Status | Summary |
|---|---|---|---|
| D-1 | U, A | fixed | Zero-gradient rows (Li-foil Φ₁; agglomerate-center c, Φ₁, Φ₂) double the old gradient instead of zeroing it |
| D-2 | U | fixed | Cathode porosity on the separator faces next to the boundary nodes |
| D-3 | U, A | fixed | No voltage cutoff: U ends on a NaN; A ends on a surface-lithiation limit (x ≥ 0.55 at one agglomerate) |
| D-4 | U, A | fixed | Single-precision literals and implicitly single-precision or integer variables (`ex_1`/`ex_curr`, `Vint`, `PI`, `last_write_time`) |
| D-6 | U, A | fixed | Finite-difference reaction derivatives with absolute steps of 10⁻⁶ |
| D-7 | U, A | fixed | One linearized solve per time step |
| D-10 | U | fixed (not ported) | Out-of-bounds read at index (SEP_NODE−NJ)/2 in the output routine |
| D-11 | U, A | fixed | Diffusion current missing from the ionic-current residual |
| D-12 | U | fixed | Lithium-foil overpotential (RT/F)·ln(I/i₀): singular at zero current, half the Butler–Volmer slope |
| D-14 | A | fixed | Separator uses the same diffusivity for both ions (t₊ = 0.5) while the cathode uses t₊ = 0.375 |
| D-15 | A | fixed (author's decision) | Two reaction descriptions for one particle: electrode kinetics on solid-sphere agglomerates fix the uptake, and the agglomerate redistributes it with its own kinetics |
| D-16 | U, A | kept (reviewed); foil term fixed | The OCP contains a Nernst term in the electrolyte concentration. This is consistent with the electrostatic Φ₂ of the Nernst–Planck equations, not double counting. The cell voltage, however, lacks the matching term of the lithium foil (fixed in corrected mode by referencing 0 V to the foil) |
| D-17 | A | fixed | Agglomerate surface electrolyte fixed at c_bulk instead of the local electrode concentration |
| D-18 | A | fixed | The electrode and agglomerate scales are solved one after the other, not together |
| D-19 | A | fixed | The output timer `last_write_time` is never initialized, so which rows are written depends on leftover memory |
| D-20 | U, A | fixed | The OCP's Redlich-Kister sum evaluates 0/0 in its k = 0 term, 2θk(1−θ)/(2θ−1)^(1−k), at θ = 1/2 exactly, so U is NaN there: a corrected run started at θ = 1/2 stopped at once with `solver_fail` (found 2026-10-04). The term is zero for every θ; corrected mode evaluates the k = 0 term as 2θ − 1 alone, which is bit-identical elsewhere |

## Evidence (T4, instrumented private copies)

**Agglomerate model.** One conductivity-study snapshot (σ = 6.4 × 10⁻⁵ S/cm, 1.29C, 200 µm), stopped after 3000 steps (t = 130 s):

- **D-1 (latent).** The center-node jumps of c, Φ₁ and Φ₂ in every agglomerate are exactly 0: the initial state is uniform and the doubling rows multiply an exact zero. As in lfp-model, it breaks an iterated (Newton) solver.
- **D-17.** The electrode electrolyte concentration is 2.6 × 10⁻³ mol/cm³ at the separator interface and **1.1 × 10⁻¹³** at the current collector: the back of the electrode is fully depleted after 130 s. The agglomerates there still see c_bulk = 10⁻³ at their surface and keep reacting. The step reduction for a depleted collector node (Δt/10) is active. Salt is not conserved globally, because the fixed-concentration surface condition is a source of salt.
- **D-11.** At a mid-cathode face, i₁ + i₂,migration = **1.000026 I**, which is the quantity the model conserves. The full dilute-solution current, i₁ + i₂,migration + i₂,diffusion, is **0.856 I**. The diffusion current (−0.14 I) is missing from the residual.
- **D-15/D-18.** At every sampled cathode node (5 nodes from the interface to the collector), the electrode-scale c_s equals the agglomerate volume-average c_s to 7 significant digits. The coupling therefore conserves lithium between the scales. The remaining issue is modelling consistency: the electrode kinetics use the agglomerate-surface crystal concentration with the *electrode* potentials and area, while the agglomerate uses its own potentials and crystal area.
- **D-3.** The run's only exits are a surface-lithiation limit at the interface-node agglomerate, a NaN, or 99 h. There is no voltage cutoff.

**Uniform-particle model.** Rebuilt from the "questions" copy, it reproduces the archived output byte for byte and ends with `EXIT BECAUSE delC ISNAN` (D-3). Its equations are the lfp-model program's, line for line, apart from the parameters and OCP. D-1, D-2, D-4, D-6, D-7, D-10, D-11 and D-12 are therefore present in the same form; their effects were measured in lfp-model and will be re-measured with the NMC parameters in T9/T10.


**D-19 (agglomerate output timer never initialized).** `last_write_time` (implicitly an integer) is read at A:L646 before it is ever assigned (A:L648 is its only assignment). On the HPC cluster it evidently started at 0, since the archived runs wrote a row every step. Rebuilt on macOS arm64 with the same compiler flags, 6 of the 10 snapshots started with a large leftover value. They wrote no rows until simulated time exceeded it, leaving output files with a single data row before the end. The simulation itself is unaffected; only the choice of rows written changes. Reference builds therefore use `-finit-local-zero`. Faithful mode starts the timer at 0, and corrected mode uses a proper output schedule.
