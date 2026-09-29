"""Generate the tutorial notebooks (committed without outputs; CI executes them)."""
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parents[1] / "notebooks"


def nb(cells):
    n = nbf.v4.new_notebook()
    n.cells = [nbf.v4.new_markdown_cell(c[1]) if c[0] == "md" else nbf.v4.new_code_cell(c[1]) for c in cells]
    n.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    return n


quickstart = [
("md", """# 1 · Quickstart

This notebook runs the NMC111 half-cell models from Python. Each model is a porous NMC111 cathode against a lithium-foil counter electrode, with a separator between them. There are two particle models:

- **uniform**: each electrode node holds uniform-concentration particles;
- **agglomerate**: each electrode node holds porous agglomerates of small crystals, with their own electrolyte, conduction and kinetics inside (docs/model.md section 6).

The same models are available as self-contained Fortran (`fortran/nmc.f90`) and C++ (`cpp/nmc.cpp`) programs, which read the same input files.

All parameters are generic defaults for illustration (docs/parameters.md). They are not fitted to any particular cell."""),
("code", """import numpy as np
import matplotlib.pyplot as plt

from nmc_model.inputs import load
from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import run_corrected
from nmc_model.uniform.params import Params
from nmc_model.uniform.simulate import run as run_uniform"""),
("md", """## Input files

`input/uniform.nml` and `input/agglomerate.nml` list every parameter at its default. `load` reads them into parameter objects; `python -m nmc_model input/agglomerate.nml` runs one and writes `Time_Voltage.txt`."""),
("code", """p_u, model, _ = load("../input/uniform.nml")
p_a, model_a, _ = load("../input/agglomerate.nml")
print(model, "C-rate", p_u.C_rate, "| thickness", p_u.L_cath_um, "um")
print(model_a, "C-rate", p_a.C_rate, "| thickness", p_a.thickness_um, "um, loading", p_a.mass_loading * 1e3, "mg/cm2")"""),
("md", """## A 1C discharge with each model

The agglomerate model is run here on a coarser grid than the input file's (30 electrode nodes, 12 along each agglomerate radius, 5 s steps) so that it takes a few seconds. Both runs stop at their lower voltage cutoff."""),
("code", """r_u = run_uniform(Params(mode="corrected", C_rate=1.0))
coarse = dict(nj=30, nja=12, dt_s=5.0)
r_a = run_corrected(AggParams(mode="corrected", C_rate=1.0, **coarse))
print("uniform:", r_u.exit_reason, "| agglomerate:", r_a.exit_reason)

def mAhg(r, M):
    return r.array[:, 2] * 96485.0 / (M * 3.6)          # electron equivalents -> mAh/g

plt.plot(mAhg(r_u, p_u.M), r_u.array[:, 1], label="uniform particles")
plt.plot(mAhg(r_a, p_a.M), r_a.array[:, 1], label="agglomerates")
plt.xlabel("capacity [mAh/g]"); plt.ylabel("cell voltage vs Li [V]"); plt.legend(); plt.title("1C discharge")
plt.show()"""),
("md", """The two models have different open-circuit-potential fits and electrode designs (from the original studies), so the curves are not a like-for-like comparison. Notebook 2 looks inside the agglomerate model.

## The output file

Columns written in corrected mode (docs/protocol.md): state, time [h], cell voltage against the lithium foil [V], electron equivalents, the foil overpotential [mV], its exchange current [mA/cm²], the electrolyte concentration at the foil, the current [mA/cm²], the protocol step, and the foil's Nernst potential [mV]. The agglomerate model adds the surface lithiation of the agglomerate next to the separator (`x_front`) and the electrolyte concentration at the current collector."""),
("code", """import os, tempfile
path = os.path.join(tempfile.mkdtemp(), "Time_Voltage.txt")
r_a.write(path)
print(open(path).read().splitlines()[0])
print(open(path).read().splitlines()[2])"""),
("md", """## Protocols

A protocol is a list of steps: constant current (`cc`), constant voltage (`cv`) and rest. Here: discharge at 1C, rest 10 minutes, charge at 1C to 4.3 V, hold 4.3 V until the current falls to C/20, rest."""),
("code", """cycle = "cc C=1 Vmin=3.0; rest t=600; cc C=-1 Vmax=4.3; cv V=4.3 Imin=0.05; rest t=600"
r = run_corrected(AggParams(mode="corrected", steps=cycle, **coarse))
a = r.array
fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(a[:, 0], a[:, 1]); ax[0].set_ylabel("V [V]")
ax[1].plot(a[:, 0], a[:, 6]); ax[1].set_ylabel("I [mA/cm²]"); ax[1].set_xlabel("time [h]")
plt.show()
print(r.exit_reason)"""),
]

inside = [
("md", """# 2 · Inside the agglomerate model

In the agglomerate model each cathode node holds porous spheres (agglomerates, radius 5 µm) made of small crystals (radius 200 nm). The electrolyte fills the pores between agglomerates and the pores inside them. The two scales are solved together (docs/model.md section 6):

- the **electrode** carries salt and current between the foil and the current collector;
- each **agglomerate** takes electrolyte, current and electrons in through its surface, and its crystals react with the electrolyte in its pores.

This notebook looks at the state at the end of a 2C discharge, with the default (Bruggeman) transport inside the agglomerates and with a much slower one."""),
("code", """import numpy as np
import matplotlib.pyplot as plt

from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import CorrectedModel, run_corrected

p = AggParams(mode="corrected", C_rate=2.0, nj=30, nja=12, dt_s=5.0)
slow = AggParams(mode="corrected", C_rate=2.0, nj=30, nja=12, dt_s=5.0, tau_agg=2240.0)   # 1000x Bruggeman
r, r_slow = run_corrected(p), run_corrected(slow)
m = CorrectedModel(p)
st, st_slow = r.final_state, r_slow.final_state
for name, rr in (("default", r), ("slow agglomerate transport", r_slow)):
    print(name, ":", rr.exit_reason, "after", round(rr.array[-1, 0] * 60, 1), "min")"""),
("md", """## Across the electrode

The electrolyte concentration in the macro-pores, from the foil (left) through the separator to the current collector (right)."""),
("code", """x = m.mesh.x * 1e4
plt.plot(x, m.conc(st)[0] * 1e3)
plt.axvline(p.L_sep * 1e4, color="gray", lw=0.5)
plt.xlabel("x [µm]"); plt.ylabel("c [mol/L]"); plt.title("electrolyte at the end of a 2C discharge")
plt.show()"""),
("md", """## Inside the agglomerates

The crystal lithiation x = c_s/(ρ/M) along the radius of the agglomerates at the front (next to the separator), the middle and the back (next to the current collector). The center is r = 0, the surface r = R.

With the default transport the salt reaches the agglomerate centers easily, and every crystal lithiates almost equally. With 1000 times slower transport in the agglomerate pores, the salt inside is depleted: the crystals near the surface fill first, and the discharge ends early."""),
("code", """r_um = m.xa * 1e4
fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for a_, s_, title in ((ax[0], st, "default"), (ax[1], st_slow, "slow agglomerate transport")):
    for k, label in ((0, "front"), (m.nl // 2, "middle"), (m.nl - 1, "back")):
        a_.plot(r_um, m.cs(s_)[k] / p.mol_vol, label=label)
    a_.set_xlabel("r [µm]"); a_.set_title(title)
ax[0].set_ylabel("x in LixNMC"); ax[0].legend(); plt.show()"""),
("md", """And the electrolyte inside the same agglomerates: the pores inside an agglomerate are fed only through its surface."""),
("code", """fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for a_, s_, title in ((ax[0], st, "default"), (ax[1], st_slow, "slow agglomerate transport")):
    for k, label in ((0, "front"), (m.nl // 2, "middle"), (m.nl - 1, "back")):
        a_.plot(r_um, m.conc(s_)[1][k] * 1e3, label=label)
    a_.set_xlabel("r [µm]"); a_.set_title(title)
ax[0].set_ylabel("c in the agglomerate pores [mol/L]"); ax[0].legend(); plt.show()"""),
("md", """## Conservation

The coupling conserves lithium and salt exactly: the lithium in the crystals equals the charge passed, and the salt in both pore scales is constant (tests/test_agglomerate_corrected.py)."""),
("code", """v_one = 4 / 3 * np.pi * p.R_agg ** 3
li = ((((m.cs(st) - p.cs_init) * m.dV).sum(axis=1) * (1 - p.eps_agg) / v_one * m.v_agg) * m.mesh.dx[m.nodes]).sum()
passed = r.array[-1, 0] * 3600 * 2.0 * p.i_1C / p.F
print("lithium in crystals / charge passed =", li / passed)"""),
]

sensitivity = [
("md", """# 3 · Sensitivity to conductivity, kinetics and agglomerate transport

The public defaults are generic (docs/parameters.md). This notebook shows how the rate capability depends on three of them:

- the electronic conductivity σ;
- the rate constant k;
- transport inside the agglomerates, varied through the tortuosity of their pores.

Each case is a discharge to 3.0 V at 1C, 3C and 10C on a coarse grid. With the defaults the cell is not strongly rate-limited below about 2C, so the ranges below reach into the limiting regimes."""),
("code", """import numpy as np
import matplotlib.pyplot as plt
from dataclasses import replace

from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import run_corrected

base = AggParams(mode="corrected", nj=30, nja=10, dt_s=10.0)
rates = [1.0, 3.0, 10.0]

def capacity(p):
    r = run_corrected(p)
    return r.array[-1, 2] * 96485.0 / (p.M * 3.6)      # mAh/g at the cutoff

def sweep(name, values, label):
    for v in values:
        caps = [capacity(replace(base, C_rate=c, **{name: v})) for c in rates]
        plt.plot(rates, caps, "o-", label=f"{label} = {v:g}")
    plt.xlabel("C-rate"); plt.ylabel("capacity to 3.0 V [mAh/g]"); plt.legend(); plt.show()"""),
("md", """## Electronic conductivity"""),
("code", """sweep("sigma", [0.001, 0.01, 0.1], "σ [S/cm]")"""),
("md", """## Rate constant"""),
("code", """sweep("rxn_k", [2.5e-8, 2.5e-7, 2.5e-6], "k")"""),
("md", """## Transport inside the agglomerates

The default pore tortuosity is Bruggeman's, ε_agg^−½ ≈ 2.2. Larger values slow salt transport inside the agglomerates; around 1000 times Bruggeman's value the agglomerate pores become the bottleneck, their cores run out of salt, and the discharge ends early. The model handles the exhausted cores smoothly (docs/model.md section 8)."""),
("code", """sweep("tau_agg", [2.24, 224.0, 2240.0, 22400.0], "τ_agg")"""),
]

for name, cells in (("01_quickstart", quickstart), ("02_inside_the_agglomerate", inside), ("03_sensitivity", sensitivity)):
    HERE.mkdir(exist_ok=True)
    nbf.write(nb(cells), HERE / f"{name}.ipynb")
print("wrote notebooks")
