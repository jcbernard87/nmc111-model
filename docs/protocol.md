# Cycling protocols

In corrected mode the cell is driven by a **protocol**: a list of steps that run in order, optionally repeated. Faithful mode ignores the protocol and runs the original constant-current discharge. Both particle models use the same protocols (the Python driver is `nmc_model/driver.py`).

The protocol is given in the `&protocol` group of the input file, read identically by the Fortran, C++ and Python implementations:

```fortran
&protocol
  steps  = 'cc C=0.5 Vmin=3.0; rest t=1800; cc C=-0.5 Vmax=4.3; cv V=4.3 Imin=0.05; rest t=1800'
  cycles = 3
/
```

Without a `&protocol` group (or with `steps = ''`) the protocol is a single step, `cc C=<C_rate> Vmin=<V_min> Vmax=<V_max>`, i.e. a constant-current discharge.

## Steps

Steps are separated by `;`. Each step is a keyword followed by `key=value` settings (case-insensitive, any order):

| Step | Settings | Ends when |
|---|---|---|
| `cc` | `C` = C-rate (1/h), **positive = discharge, negative = charge** (required); `t` = maximum duration [s]; `Vmin`, `Vmax` [V] (default: the global `V_min`, `V_max`) | the voltage reaches `Vmin` (discharge) or `Vmax` (charge), or `t` has elapsed |
| `cv` | `V` = held voltage [V] (required); `t` = maximum duration [s]; `Imin` = current magnitude, as a C-rate, below which the step ends | \|I\| ≤ `Imin`, or `t` has elapsed |
| `rest` | `t` = duration [s] (required) | `t` has elapsed |

A step with no possible end (for example `cv` with neither `t` nor `Imin`) is rejected. `cycles = n` runs the whole list n times.

**C-rate.** 1C is the current that would pass the theoretical capacity of the cathode's active material in one hour: I = Q_th × (active mass per area). For the uniform model the active mass is L_cath·ε_AM·ρ; for the agglomerate model it is `percent_active`·`mass_loading`.

## How steps are solved

- **Time step.** Every step uses a fixed Δt: `t_max/n_steps` for the uniform model (1 s by default) and `dt_s` for the agglomerate model (1 s by default). A step's last time step is shortened so the step ends exactly at `t`.
- **Newton.** Each time step is solved to convergence by Newton's method (docs/model.md section 8). If Newton fails, the sub-step is halved (down to 10⁻¹⁰ s), and after each success it doubles again up to Δt. After 200 failures within one time step the run stops. Near a full particle the voltage falls steeply, so the cutoff is approached with geometrically shrinking sub-steps.
- **Voltage cutoffs.** A sub-step that crosses a cutoff by more than 0.1 mV is repeated with half the length, so a `cc` step ends within 0.1 mV of its cutoff.
- **Constant current and rest** use the applied current directly (rest: I = 0).
- **Constant voltage.** The current is unknown: each time step finds the current I for which the cell voltage equals `V`, to |V − V_set| ≤ 10⁻⁹ V. The root is bracketed, starting from the previous step's current and expanding through I = 0, and then found by the Illinois variant of regula falsi. A current the cell cannot carry for the whole time step (Newton fails) counts as lying beyond the set voltage, which keeps the search monotone.
- **Cell voltage.** The voltage between the current collector and the lithium foil, V = Φ₁(collector) − U_Li − η_Li(I), with the foil at 0 V (docs/model.md section 7). U_Li = (RT/F)·ln(c(0)/c_Li,ref) is the foil's Nernst potential. The lithium counter electrode is a symmetric Butler–Volmer interface, η_Li = (RT/(αF))·asinh(I/(2 i₀,Li)) with α = 0.5 and i₀,Li from the electrolyte concentration at x = 0 (deviation D-12).

## Output

In corrected mode `Time_Voltage.txt` has these columns:

| Column | Header | Meaning |
|---|---|---|
| 1 | State | `D` (I > 0), `C` (I < 0) or `R` (I = 0) |
| 2 | Time (hours) | |
| 3 | Voltage (Volts) | the cell voltage against the lithium foil |
| 4 | Equivalence | electron equivalents passed per formula unit |
| 5 | Anode_Eta (mV) | −η_Li, negative on discharge |
| 6 | anode_exchange_c (mA/cm2) | i₀,Li |
| 7 | Edge_c0 (mol/cm3) | electrolyte concentration at the foil |
| 8 | Current (mA/cm2) | applied current density, positive on discharge |
| 9 | Step | 1-based index of the step in the expanded list (cycles × steps) |
| 10 | Li_Nernst (mV) | the foil's Nernst potential U_Li |
| 11 | x_front (LixNMC) | agglomerate model only: surface-crystal lithiation of the agglomerate next to the separator |
| 12 | c_collector (mol/cm3) | agglomerate model only: electrolyte concentration at the current collector |

Rows are written every `write_interval` seconds (default 18 s) and at the start and end of every step.

The run's exit reason is `end_of_protocol` when every step completed, or `nan` or `max_time` (99 h) otherwise. When a time step cannot be solved, the reason names the physical limit the cell has reached: `electrolyte_depleted` (the electrolyte in the electrode or in the agglomerate pores has fallen below 10⁻³ of c_bulk somewhere), `particles_full` or `particles_empty` (some particles or crystals within 10⁻³ of full or empty); `solver_fail` if none applies. The last row then reports the state at the start of that time step. A single-step protocol reports how that step ended (`cutoff_low`, `cutoff_high`, `duration` or `current_limit`).
