# QMOMKit paper studies

`paper/verification/` contains standalone, explicitly run numerical studies for
the QMOMKit paper. They complement, but do not replace, the automated tests
under `tests/`. `paper/examples/` contains standalone copies of five existing
application examples: `examples/cqmom_non-spherical-particles/`,
`examples/cqmom_condensing_droplets/`, and
`examples/cqmom_polydisperse_cavitation/`,
`examples/cqmom_spherical_breakup_quadrature_4d/`, and
`examples/cqmom_spherical_breakup_moments_5d/`.

Each paper example preserves its canonical physics, configuration, fixed-seed
Monte Carlo sample, and final time. One sampled ensemble initializes QMOMKit
and is then evolved as the MC reference. For the controlled comparison, MC
uses the same outer timestep and continuous RK23 trajectory options as QMOMKit
node advection, and both paths are output at identical physical times. This
deliberately differs from a standalone canonical MC executable if its timestep
or integration options differ.

## Paper application examples

Install QMOMKit and QMOMKitMonteCarlo, then configure the paper examples as a
standalone downstream consumer:

```bash
cmake -S paper/examples -B <paper-examples-build-dir> \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=<install-prefix>
cmake --build <paper-examples-build-dir>
```

Run each executable from its corresponding build subdirectory so it uses the
copied case configuration:

```bash
cd <paper-examples-build-dir>/non_spherical_particles
../qmomkit_paper_non_spherical_particles

cd <paper-examples-build-dir>/condensing_droplets
../qmomkit_paper_condensing_droplets

cd <paper-examples-build-dir>/polydisperse_cavitation
../qmomkit_paper_polydisperse_cavitation

cd <paper-examples-build-dir>/cqmom_spherical_breakup_quadrature_4d
../qmomkit_paper_cqmom_spherical_breakup_quadrature_4d

cd <paper-examples-build-dir>/cqmom_spherical_breakup_moments_5d
../qmomkit_paper_cqmom_spherical_breakup_moments_5d
```

Each case writes its existing example outputs plus
`output/moment_comparison.csv`, with the aligned columns `time`, `quantity`,
`qmomkit`, and `monte_carlo`.

The paper quantities are:

- Non-spherical particles `(d, chi, u, v, w)`: `M_00000`,
  `M_30000/M_20000`, `M_01000/M_00000`, `M_30100/M_00000`,
  `M_30010/M_00000`, `M_30001/M_00000`, `Theta_x`, `Theta_y`, and
  `Theta_z`, where `Theta_x = M_00200/M_00000 - (M_00100/M_00000)^2`
  with analogous `y` and `z` formulas.
- Condensing droplets `(r, w, T)`: `M_000`, `M_300/M_200`,
  `M_010/M_000`, and `M_001/M_000`. Here `r` is radius; the requested
  `M_300/M_200` is plotted exactly, without an added factor of two.
- Polydisperse cavitation `(R, U)`: `M_00`, `M_10/M_00`, `M_01/M_00`,
  `M_11/M_00`, `M_20/M_00`, and `M_02/M_00`.
- Spherical breakup with recoil `(d, u, v, w)`: raw `M[0,0,0,0]`,
  `M[3,0,0,0]`, `M[1,0,0,0]`, `M[2,0,0,0]`, `M[0,1,0,0]`,
  `M[0,2,0,0]`, `M[2,1,0,0]`, `M[2,2,0,0]`, and `M[6,0,0,0]`.
  This case demonstrates terminal-quadrature closure for nonlinear,
  diameter-dependent recoil.
- Thermal spherical breakup `(V, T, u, v, w)`: raw `M[0,0,0,0,0]`,
  `M[1,0,0,0,0]`, `M[0,1,0,0,0]`, `M[0,2,0,0,0]`,
  `M[1,1,0,0,0]`, `M[1/3,0,0,0,0]`, `M[4/3,0,0,0,0]`,
  `M[1/3,1,0,0,0]`, and `M[4/3,1,0,0,0]`. This case demonstrates
  generalized-moment evaluation of fractional-volume moments.


NumPy and Matplotlib are required to produce the figures:

```bash
cd <paper-examples-build-dir>/non_spherical_particles
python3 <source-dir>/paper/examples/non_spherical_particles/plot.py

cd <paper-examples-build-dir>/condensing_droplets
python3 <source-dir>/paper/examples/condensing_droplets/plot.py

cd <paper-examples-build-dir>/polydisperse_cavitation
python3 <source-dir>/paper/examples/polydisperse_cavitation/plot.py
```

The scripts write `non_spherical_particles_moments.pdf`,
`condensing_droplets_moments.pdf`, and
`polydisperse_cavitation_moments.pdf`, respectively, in the corresponding
build directories.

Install QMOMKit, then configure the verification project in an out-of-source
directory:

```bash
cmake -S paper/verification -B <verification-build-dir> \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=<qmomkit-install-prefix>
cmake --build <verification-build-dir>
```

Run the studies from that build directory so that generated files remain
outside the source tree:

```bash
cd <verification-build-dir>
./qmomkit_paper_quadrature_1d
./qmomkit_paper_cqmom_affine_chain
./qmomkit_paper_adaptive_branch_quadrature
./qmomkit_paper_time_convergence
```

They write `quadrature_1d.csv`, `cqmom_affine_chain.csv`,
`adaptive_branch_quadrature.csv`, `rk_convergence.csv`, and
`splitting_convergence.csv` in the current directory.

The adaptive study loads its two-dimensional `(d,u)` hierarchy, including both
node counts and conditioning order, from `configuration.qmom`. The first node
count determines the number of size branches, and the second is the maximum
conditional velocity-node count for each branch. Configurations such as
`{1,5}`, `{2,3}`, and `{2,2}` require no source changes. Its CSV uses one row
per output time and size branch with `time`, `total_mass`, `mean_velocity`,
`velocity_variance`, `branch_index`, and `active_nodes`. The relaxation,
reinjection pulse, final time, and timestep remain fixed C++ definitions because
they specify the benchmark physics. The configured splitting, node-advection,
and output options are loaded through QMOMKit.

Create the three figures from the same build directory with NumPy and
Matplotlib:

```bash
python3 <source-dir>/paper/verification/adaptive_branch_quadrature/plot.py
python3 <source-dir>/paper/verification/time_convergence/plot.py
```

The scripts write `adaptive_branch_quadrature.pdf`, `rk_convergence.pdf`, and
`splitting_convergence.pdf` in the current directory.
