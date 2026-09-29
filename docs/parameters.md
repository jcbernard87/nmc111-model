# NMC111 model parameters

These are the values in the original research code. Line references are as in [model.md](model.md). "Stored" gives the single-precision rounding that the original actually used, where the value was written as a default-kind literal **[D-4]**.

**Fitted values are not distributed.** Some agglomerate-model values were fitted to experiments on a specific cell. This repository shares models, not results or cell data, so those values are not included, and the table marks them **private**. Faithful mode of the agglomerate model therefore needs them as inputs (`D_agg`, `k_rxn`, `tortuosity_e`, `mass_loading`); the programs stop if they are missing. The public defaults below are generic values.

## Uniform-particle model (`NMC111_Electrode_Crystal.f95`)

Same parameter set and meaning as lfp-model ([its parameters.md](https://github.com/jcbernard87/lfp-model/blob/main/docs/parameters.md)), except:

| Name | Value | Unit | Source | Notes |
|---|---|---|---|---|
| `NJ` | 43 | | C:L35 | 20 separator and 20 cathode volumes |
| `molar_mass_AM` | 96.46 | g/mol | C:L71 | NMC111 |
| `density_AM` | 4.6 | g/cm³ | C:L72 | |
| `Qmth` | 0.150 | Ah/g | C:L73 | |
| `sigma` | private (fitted) | S/cm | C:L122 | faithful mode needs it as input; public default 0.1 S/cm |
| `rxn_k` | private (fitted), 10^x evaluated in single precision | | C:L98 | faithful mode needs it as input (`k_rxn`); public default 2.5 × 10⁻⁶ |
| `volfrac_AM` | 0.4 | | C:L118 | with ε = 0.5 |
| `sep_tortuosity` | 4.8 | | C:L110 | |
| `Phi_1_init` | 4.1 | V | C:L79 | |
| `diff_c` | 5.162 × 10⁻¹² | cm²/s | C:L97 | used only by the inactive crystal scale |
| OCP | 11-term Redlich–Kister, U_ref = 3.8685682447595453, x_max = 0.55 | | C:L145–L179 | see model.md §2 |
| `cimax` | 0.55·ρ/M | mol/cm³ | C:L197 | |
| `C_rate` | 1.0 | 1/h | C:L123 | |

**Public defaults of the uniform model** (corrected mode): σ = 0.1 S/cm, mid-range for carbon-containing NMC electrodes; k = 2.5 × 10⁻⁶, which gives i₀ ≈ 0.1 mA/cm² at θ = 0.5 and 1 M. Other values are the original's (above).

## Agglomerate model (`NMC111_agg.f95`, conductivity-study template)

| Name | Symbol | Value (stored) | Unit | Source | Notes |
|---|---|---|---|---|---|
| `Rigc`, `Temp`, `Fconst` | R, T, F | 8.314 (8.314000129699707), 298, 96485 | | A:L36 | |
| `PI` | π | single precision (3.1415927410125732) | | A:L37 | used in the agglomerate geometry **[D-4]** |
| `electrode_area` | | 1.13 | cm² | A:L42 | energy output only |
| `eps` | ε | per run (`POROSITY`) | | A:L120 | cathode porosity |
| `eps_sep` | ε_sep | 0.39 | | A:L45 | |
| `diff_0` | D₀ | 2.89 × 10⁻⁶ | cm²/s | A:L50 | salt diffusivity |
| `tortuosity` | τ | private (fitted) | | A:L123 | cathode; D_e = εD₀/τ; faithful mode only (input `tortuosity_e`) |
| `tau_sep` | τ_sep | 4.0 | | A:L61 | |
| `transference_num_cat` | t₊ | 0.375 | | A:L59 | cathode only; the separator effectively uses 0.5 **[D-14]** |
| `cbulk`, `c0_init` | c⁰ | 0.001 (0.0010000000474974513) | mol/cm³ | A:L53, L77 | |
| `sigma` | σ | per run (`CONDUCTIVITY`, a single-precision literal) | S/cm | A:L125 | swept 10⁻⁴–8 S/cm in the study |
| `Massloading` | | private (cell value) | g/cm² | A:L121 | |
| `THICKNESS` | L | per run, in µm (`THICKNESS/10000.0` is a single-precision division) | cm | A:L122 | |
| `percent_active` | | 0.95 | | A:L124 | |
| derived `v_bar_AM` | v_AM | percent_active·Massloading/(ρ·L) | | A:L139 | active volume fraction |
| `molar_mass_NMC`, `density_NMC` | M, ρ | 96.46, 4.7 | g/mol, g/cm³ | A:L69 | |
| `mol_vol` | | 0.0476881609 | mol/cm³ | A:L70 | a literal, not ρ/M |
| `Qmth` | Q_th | 0.155 | Ah/g | A:L73 | |
| `xmax_c` | R_agg | 5 × 10⁻⁴ | cm | A:L135 | agglomerate radius |
| `eps_agg` | ε_agg | 0.2 | | A:L64 | agglomerate porosity |
| `diff_agg` | D_agg | private (fitted) | cm²/s | A:L128 | electrolyte diffusivity inside the agglomerate; faithful mode only |
| `xmax_xtal` | r_crystal | 200 × 10⁻⁷ | cm | A:L1724 | crystal radius inside agglomerates |
| `rxn_k` | k | private (fitted), evaluated in single precision | | A:L129 | |
| `alpha_a`, `alpha_c` | | 0.5, 0.5 | | A:L205 | |
| `Phi_1_init`, `Phi_2_init`, `cs_init` | | 4.3, 0, 10⁻⁵ | V, V, mol/cm³ | A:L77 | |
| `SEP_NODE`, `len_sep` | | 22, 25 × 10⁻⁴ cm | | A:L99–L100 | |
| `NJ`, `NJ_c` | | per run (`NUMJ`, `NUMC`, written as reals such as `68.0`) | | A:L32 | |
| `tmax` | | 72 000 | s | A:L32 | |
| `Time_mod` | | per run (`TIMEMOD`) | | A:L111 | N_steps = ⌊3600·C_rate·Time_mod⌋ |
| `C_rate` | | per run (`CRATE`, bisection values exact in single precision) | 1/h | A:L119 | |
| `R_contact` | | 0 | Ω·cm² | A:L130 | |
| current ramp | | I/50, ×1.5 per step | | A:L558–L593 | |

## Agglomerate model: public defaults (corrected mode)

The defaults of `AggParams` and of the `&agglomerate` input. Values not listed are the original template's (above), which are not cell fits: the OCP fit, the agglomerate and crystal radii, ε_agg, D₀, t₊, M, ρ, Q_th and the separator.

| Input | Default | Why |
|---|---|---|
| `L_cath_um` | 100 µm | a round, typical laboratory cathode thickness |
| `mass_loading` | 0.020 g/cm² | a round, typical laboratory loading; with 95 % active material and ε = 0.4 the solid and pore volumes add up to 0.905 |
| `eps` | 0.4 | the template's value |
| `sigma` | 0.1 S/cm | mid-range for carbon-containing NMC electrodes (the original study swept 10⁻⁴ to 8 S/cm) |
| `bruggeman` | −0.5 | electrode tortuosity τ = ε^−½ (Bruggeman) instead of a fitted τ |
| `tau_agg` | ε_agg^−½ | Bruggeman inside the agglomerates; the agglomerate pores carry the same electrolyte (D₊, D₋ from D₀ and t₊) |
| `sigma_agg` | = `sigma` | the agglomerate solid conducts like the electrode solid |
| `k_rxn` | 2.5 × 10⁻⁶ | gives i₀ ≈ 0.1 mA/cm² per crystal area at θ = 0.5 and 1 M, within the 10⁻⁵ to 10⁻³ A/cm² range usual for NMC111 |
| `dt_s` | 1 s | time step (Newton is converged at every step, so it controls only time accuracy) |
| `V_min`, `V_max` | 3.0 V, 4.4 V | NMC111 half-cell limits; the initial open-circuit voltage is 4.37 V |
| `k_Li`, `c_Li_ref` | 10⁻⁶, 10⁻³ mol/cm³ | the lithium foil, as in lfp-model |
| `kappa_bg` | 10⁻⁸ S/cm | the solvent's own ionic conductivity (neat carbonates are about 10⁻⁹ to 10⁻⁷ S/cm); keeps Φ₂ defined where the salt is exhausted (docs/model.md section 8). Both models, corrected mode |

The notebook `agglomerate_vs_uniform` shows how the results depend on σ, k and the agglomerate transport.
