# nmc111-model

A one-dimensional model of an NMC111 half cell (lithium foil | separator | porous NMC111 cathode | current collector) with two particle models, solved with [bandsolver](https://github.com/jcbernard87/bandsolver), an implementation of Newman's BAND method.

- **Uniform particles** (`particle_model = 'uniform'`): each electrode node holds uniform-concentration particles.
- **Agglomerates** (`particle_model = 'agglomerate'`): each electrode node holds porous spherical agglomerates of small crystals. The agglomerates have their own electrolyte, electronic conduction and reaction inside. The two scales are coupled through the agglomerate surfaces and solved together, with one Newton iteration for both scales at every step.

Both use dilute-solution transport of a binary electrolyte, electronic conduction and Butler–Volmer insertion kinetics. They are discretized with finite volumes, advanced by backward Euler with Newton's method at every step, and run constant-current, constant-voltage and rest protocols, including cycling.

The models come in three implementations that read the same input files and write the same output:

| Implementation | Where | Solver | Use it for |
|---|---|---|---|
| Fortran program | [`fortran/nmc.f90`](fortran/nmc.f90) (one file) | links bandsolver | batch runs |
| C++ program | [`cpp/nmc.cpp`](cpp/nmc.cpp) (one file) | its own embedded copy of bandsolver's C++ core; no dependencies | batch runs, embedding |
| Python package | [`python/nmc_model`](python/nmc_model), [`notebooks/`](notebooks) | bandsolver (Python) | examples, quick studies, notebooks |

This repository shares the **models only**. It contains no simulation results, no experimental data and no parameters fitted to a particular cell; its defaults are generic, documented values. Run the models to generate results.

## Quick start

**C++** (only a C++17 compiler is needed):

```sh
c++ -std=c++17 -O2 -ffp-contract=off cpp/nmc.cpp -o nmc_cpp
./nmc_cpp input/agglomerate.nml        # writes Time_Voltage.txt
./nmc_cpp input/uniform.nml
```

**Fortran and C++ with CMake** (bandsolver is downloaded automatically):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
build/fortran/nmc_f input/agglomerate.nml
```

With a local bandsolver checkout, add `-DFETCHCONTENT_SOURCE_DIR_BANDSOLVER=/path/to/bandsolver`.

**Python:**

```sh
pip install "bandsolver @ git+https://github.com/jcbernard87/bandsolver@v0.1.2"
pip install -e ".[notebooks]"
python -m nmc_model input/agglomerate.nml
```

```python
from nmc_model.agglomerate import AggParams
from nmc_model.agglomerate_corrected import run_corrected

r = run_corrected(AggParams(mode="corrected", C_rate=1.0))
t_h, volts = r.array[:, 0], r.array[:, 1]
```

## Input

One Fortran-namelist file holds every parameter: [`input/uniform.nml`](input/uniform.nml) and [`input/agglomerate.nml`](input/agglomerate.nml), each at its defaults. `&model particle_model` chooses the model. [docs/parameters.md](docs/parameters.md) lists each parameter with its unit, origin and, for the public defaults, the reason for the value. The Fortran program reads the files natively; the C++ and Python implementations parse the same format.

A cycling protocol goes in a `&protocol` group ([docs/protocol.md](docs/protocol.md)):

```fortran
&protocol
  steps  = 'cc C=1 Vmin=3.0; rest t=600; cc C=-1 Vmax=4.3; cv V=4.3 Imin=0.05; rest t=600'
  cycles = 2
/
```

The cell voltage is measured against the lithium foil.

## Two modes

The models are ports of the author's PhD research code. The port was checked for errors while it was being made:

- **`mode = 'corrected'`** (recommended, and the input files' setting) fixes every defect found. The fixes are listed with evidence in [docs/deviations.md](docs/deviations.md). For the agglomerate model this includes a consistent coupling of the two scales, the local electrolyte at the agglomerate surfaces, and a fully coupled solve (docs/model.md section 6). In both models the electrolyte can run out and particles can fill or empty smoothly, without clipping: the unknowns are the logarithm of the concentration and the log-odds of the lithiation, with exponential-fitting (Scharfetter–Gummel) fluxes (docs/model.md section 8).
- **`mode = 'faithful'`** reproduces the original programs exactly, defects included, for comparison with earlier work ([docs/validation.md](docs/validation.md)). It needs the original's fitted parameter values, which are not distributed, so it cannot run from the public inputs alone.

## Documentation

- [docs/model.md](docs/model.md): equations, boundary conditions, discretization and time stepping, for the original and the corrected models
- [docs/parameters.md](docs/parameters.md): every parameter, with value, unit and origin
- [docs/protocol.md](docs/protocol.md): protocol steps and the output columns
- [docs/deviations.md](docs/deviations.md): defects found in the original code, with evidence and fixes
- [docs/validation.md](docs/validation.md): reproduction, cross-language agreement, conservation and convergence tests, with measured numbers
- [notebooks/](notebooks): a quickstart, the agglomerate model's internal state, and sensitivity to conductivity, kinetics and agglomerate transport

## Tests

```sh
cmake --build build && python -m pytest
```

The public tests compute everything they need. They cover cross-language agreement, an independent residual, the Jacobian, conservation of lithium and salt on both scales, equilibrium, time and mesh convergence, and protocols. The tests in `private_checks/` compare against the original programs and their archived runs. They run only when `NMC_ORACLE_DIR` points at that material, which is not distributed.

## License and citation

BSD 3-Clause; see [LICENSE](LICENSE). See [CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button uses it. If you use this model, please also cite [bandsolver](https://github.com/jcbernard87/bandsolver) and Newman's method: J. Newman, *Ind. Eng. Chem. Fundam.* 7, 514 (1968); J. Newman and K. E. Thomas-Alyea, *Electrochemical Systems*, 3rd ed., Appendix C.
