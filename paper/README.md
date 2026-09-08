# QMOMKit paper studies

`paper/examples/` contains the application studies used in the QMOMKit paper,
and `paper/verification/` contains standalone verification studies. They
complement, but do not replace, the automated tests under `tests/`; the
user-facing workflows under `examples/` are separate. Together with the
plotting scripts, these directories reproduce results used in the paper.

## Application examples

- `1d_aerosols_evaporation`: 1D QMOM aerosol evaporation, compared with Monte Carlo.
- `2d_bubbles_cavitation`: 2D CQMOM bubble cavitation, compared with Monte Carlo.
- `3d_droplets_condensation-nucleation`: 3D CQMOM droplet condensation and nucleation, compared with Monte Carlo.
- `4d_particles_breakup_recoil`: 4D CQMOM particle breakup with recoil.
- `5d_particles_breakup_heat-transfer`: 5D CQMOM particle breakup with heat transfer.

## Build

Run these commands from the repository root to install QMOMKit and its Monte
Carlo companion into `install/`, then build the paper application cases:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build

cmake -S paper/examples -B paper/examples/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build paper/examples/build
```

## Run the application examples

Run these commands from the repository root:

```bash
cd paper/examples/build/1d_aerosols_evaporation
../qmomkit_paper_1d_aerosols_evaporation
cd ../../../..

cd paper/examples/build/2d_bubbles_cavitation
../qmomkit_paper_2d_bubbles_cavitation
cd ../../../..

cd paper/examples/build/3d_droplets_condensation-nucleation
../qmomkit_paper_3d_droplets_condensation_nucleation
cd ../../../..

cd paper/examples/build/4d_particles_breakup_recoil
../qmomkit_paper_4d_particles_breakup_recoil
cd ../../../..

cd paper/examples/build/5d_particles_breakup_heat-transfer
../qmomkit_paper_5d_particles_breakup_heat_transfer
cd ../../../..
```

Results are written under the corresponding case directory in
`paper/examples/build/`, including each case's `output/` directory.

## Verification studies

Build the verification project against the installed QMOMKit package:

```bash
cmake -S paper/verification -B paper/verification/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build paper/verification/build
```

Run the studies from the repository root:

```bash
cd paper/verification/build
./qmomkit_paper_quadrature_1d
./qmomkit_paper_cqmom_affine_chain
./qmomkit_paper_adaptive_branch_quadrature
./qmomkit_paper_time_convergence
cd ../../..
```

Their CSV results are written in `paper/verification/build/`.

## Generate the figures

From the repository root, after running the application and verification
studies, run:

```bash
python3 paper/plotting/plot_examples.py
python3 paper/plotting/plot_verification.py
```

The application figures are `aerosol_evaporation_moments.pdf`,
`cavitation_moments.pdf`, `condensing_droplets_moments.pdf`,
`cqmom_spherical_breakup_quadrature_4d_moments.pdf`, and
`cqmom_spherical_breakup_moments_5d_moments.pdf`. The verification figures are
`quadrature_reconstruction.pdf`, `adaptive_branch_quadrature.pdf`, and
`splitting_convergence.pdf`; all figures are written to `paper/plotting/`.
