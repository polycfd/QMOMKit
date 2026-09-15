# QMOMKit paper studies

This directory reproduces the figures and numbers in the QMOMKit paper.
`paper/examples/` holds the five application studies, `paper/verification/`
holds four standalone verification studies, and `paper/plotting/` holds the two
scripts that turn their output into the paper figures.

These studies complement the automated tests under `tests/` and are separate
from the user-facing tutorials under `examples/`.

[Reproduce everything](#reproduce-everything) is one copy-paste block that goes
from a fresh clone to all eight figures. The sections after it explain the
pieces, in case you want to run only part of the pipeline.

## What you need

- A C++17 compiler.
- CMake 3.16 or newer.
- Python 3 with `matplotlib`, for the figures only.

Nothing else is downloaded or required. Every command below is run **from the
repository root**. Results go to `install/` and the various `build/`
directories, which git ignores, and the figures go to `paper/plotting/`, where
they will show up as untracked files.

The whole pipeline takes a few minutes, about two of which are one verification
study (`qmomkit_paper_cqmom_affine_chain`); everything else runs in seconds.

## Reproduce everything

Copy this whole block (without the explanatory comments)
into a shell at the repository root. It builds
everything, runs all nine studies, and writes all eight figures.

```bash
# 1. Build and install QMOMKit and its Monte Carlo companion.
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build tools/monte_carlo/build --parallel
cmake --install tools/monte_carlo/build

# 2. Build and run the five application studies.
cmake -S paper/examples -B paper/examples/build -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build paper/examples/build --parallel

(cd paper/examples/build/1d_aerosols_evaporation            && ../qmomkit_paper_1d_aerosols_evaporation)
(cd paper/examples/build/2d_bubbles_cavitation              && ../qmomkit_paper_2d_bubbles_cavitation)
(cd paper/examples/build/3d_droplets_condensation-nucleation && ../qmomkit_paper_3d_droplets_condensation_nucleation)
(cd paper/examples/build/4d_particles_breakup_recoil        && ../qmomkit_paper_4d_particles_breakup_recoil)
(cd paper/examples/build/5d_particles_breakup_heat-transfer && ../qmomkit_paper_5d_particles_breakup_heat_transfer)

# 3. Build and run the four verification studies.
cmake -S paper/verification -B paper/verification/build -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build paper/verification/build --parallel

(cd paper/verification/build && ./qmomkit_paper_quadrature_1d)
(cd paper/verification/build && ./qmomkit_paper_adaptive_branch_quadrature)
(cd paper/verification/build && ./qmomkit_paper_time_convergence)
(cd paper/verification/build && ./qmomkit_paper_cqmom_affine_chain)  

# 4. Write the figures into paper/plotting/.
python3 paper/plotting/plot_examples.py
python3 paper/plotting/plot_verification.py
```

Each study is wrapped in `( ... )` so it runs in a subshell: the directory
change applies only inside the parentheses, and you stay at the repository root
throughout. That matters because every study reads its `configuration.qmom`
from the current directory while the executable itself sits one level up.

## Did it work?

Both plotting scripts print the path of every figure they write. When the run
succeeds you get these eight files in `paper/plotting/`:

```text
aerosol_evaporation_moments.pdf
cavitation_moments.pdf
condensing_droplets_moments.pdf
cqmom_spherical_breakup_quadrature_4d_moments.pdf
cqmom_spherical_breakup_moments_5d_moments.pdf
quadrature_reconstruction.pdf
adaptive_branch_quadrature.pdf
splitting_convergence.pdf
```

The verification studies also print their results as tables on standard output;
redirect them to a file if you want to keep those numbers.

## Application studies

Each study writes its results to `output/` inside its own directory under
`paper/examples/build/`. Re-running a study overwrites its previous output.

| Directory | Executable | Figure |
| --- | --- | --- |
| `1d_aerosols_evaporation` | `qmomkit_paper_1d_aerosols_evaporation` | `aerosol_evaporation_moments.pdf` |
| `2d_bubbles_cavitation` | `qmomkit_paper_2d_bubbles_cavitation` | `cavitation_moments.pdf` |
| `3d_droplets_condensation-nucleation` | `qmomkit_paper_3d_droplets_condensation_nucleation` | `condensing_droplets_moments.pdf` |
| `4d_particles_breakup_recoil` | `qmomkit_paper_4d_particles_breakup_recoil` | `cqmom_spherical_breakup_quadrature_4d_moments.pdf` |
| `5d_particles_breakup_heat-transfer` | `qmomkit_paper_5d_particles_breakup_heat_transfer` | `cqmom_spherical_breakup_moments_5d_moments.pdf` |

Note that the executable names use underscores throughout, while three of the
directory names contain a hyphen. Copy the commands rather than retyping them.

What each study does:

- `1d_aerosols_evaporation`: 1D QMOM aerosol evaporation at two, four, six,
  eight and ten nodes, compared with Monte Carlo.
- `2d_bubbles_cavitation`: 2D CQMOM bubble cavitation, compared with Monte Carlo.
- `3d_droplets_condensation-nucleation`: 3D CQMOM droplet condensation and
  nucleation, compared with Monte Carlo.
- `4d_particles_breakup_recoil`: 4D CQMOM particle breakup with recoil.
- `5d_particles_breakup_heat-transfer`: 5D CQMOM particle breakup with heat
  transfer, including generalized moments.

All five link against the Monte Carlo companion in `tools/monte_carlo`, so that
package must be installed before `paper/examples` will configure.

## Verification studies

These write CSV files into `paper/verification/build/`.

| Executable | Writes | Figure | Time |
| --- | --- | --- | --- |
| `qmomkit_paper_quadrature_1d` | `quadrature_1d.csv` | `quadrature_reconstruction.pdf` | seconds |
| `qmomkit_paper_adaptive_branch_quadrature` | `adaptive_branch_quadrature.csv` | `adaptive_branch_quadrature.pdf` | seconds |
| `qmomkit_paper_time_convergence` | `splitting_convergence.csv`, `rk_convergence.csv` | `splitting_convergence.pdf` | seconds |
| `qmomkit_paper_cqmom_affine_chain` | `cqmom_affine_chain.csv` | none | ~2 minutes |

`qmomkit_paper_cqmom_affine_chain` is by far the slowest and feeds numbers
quoted in the paper text rather than any figure, so you can skip it if you only
want the figures. `rk_convergence.csv` is likewise reported as a table rather
than plotted.

## Running only part of the pipeline

The plotting scripts default to the directories used above, so they need no
arguments after a standard run. Both accept overrides if your layout differs:

```bash
python3 paper/plotting/plot_examples.py --help
python3 paper/plotting/plot_verification.py --help
```

`plot_examples.py` reads `paper/examples/build/`, and `plot_verification.py`
reads the three CSV files above; each takes an `--output-dir` for the figures.
You can regenerate one group of figures without rerunning the other group's
studies, but every study a script reads must have been run at least once.

## Troubleshooting

**`Could not find a package configuration file provided by "QMOMKit"` or
`"QMOMKitMonteCarlo"`.** The install step was skipped, or `CMAKE_PREFIX_PATH`
does not point at the install prefix. Re-run step 1, and check that you are at
the repository root so that `$PWD/install` resolves correctly.

**A study exits with a configuration error.** It was launched from the wrong
directory. Each executable reads `configuration.qmom` from the current
directory, so run it from its own case directory as the commands above do.

**A plotting script reports a missing file.** The corresponding study has not
been run yet, or was run from the wrong directory so its output landed
elsewhere. Re-run that study and check that the expected file appeared under
`paper/examples/build/<case>/output/` or `paper/verification/build/`.

**`ModuleNotFoundError: No module named 'matplotlib'`.** Install it with
`python3 -m pip install matplotlib`. Only the figures need it; the studies
themselves do not.
