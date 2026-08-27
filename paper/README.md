# QMOMKit paper studies

`paper/verification/` contains standalone, explicitly run numerical studies for
the QMOMKit paper. They complement, but do not replace, the automated tests
under `tests/`. `paper/examples/` contains frozen, standalone application cases
used for the paper; the user-facing examples under `examples/` are separate.

Each paper case preserves its configuration, physics, fixed-seed Monte Carlo
sample, and final time. A single sampled ensemble initializes QMOMKit and is
then evolved as the MC reference. MC uses the same outer timestep and
continuous RK23 trajectory options as QMOMKit node advection, and both paths
are output at identical physical times. This controlled comparison may differ
from a user-facing MC executable with different timestep or integration
options.

## Paper application examples

The paper cases are:

- `1d_aerosols_evaporation`: 1D QMOM evaporation of an aerosol diameter
  distribution, with fixed-node-count reconstructions compared with its
  fixed-seed Monte Carlo reference.
- `2d_bubbles_cavitation`: 2D CQMOM evolution of a polydisperse bubble
  population in radius and velocity, compared with Monte Carlo.
- `3d_droplets_condensation-nucleation`: 3D CQMOM condensation and nucleation
  of droplets in radius, vertical velocity, and temperature, including an
  aligned Monte Carlo reference.
- `4d_particles_breakup_recoil`: 4D CQMOM spherical-particle breakup with
  nonlinear diameter-dependent recoil. It demonstrates terminal-quadrature
  closure for the nonlinear daughter-state transformation rather than a
  physical validation of a particular breakup model.
- `5d_particles_breakup_heat-transfer`: 5D CQMOM spherical-particle breakup
  with temperature and heat transfer. It demonstrates generalized-moment
  evaluation of fractional-volume moments through terminal quadrature rather
  than a physical validation of a particular breakup model.

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
cd <paper-examples-build-dir>/1d_aerosols_evaporation
../qmomkit_paper_1d_aerosols_evaporation

cd <paper-examples-build-dir>/2d_bubbles_cavitation
../qmomkit_paper_2d_bubbles_cavitation

cd <paper-examples-build-dir>/3d_droplets_condensation-nucleation
../qmomkit_paper_3d_droplets_condensation_nucleation

cd <paper-examples-build-dir>/4d_particles_breakup_recoil
../qmomkit_paper_4d_particles_breakup_recoil

cd <paper-examples-build-dir>/5d_particles_breakup_heat-transfer
../qmomkit_paper_5d_particles_breakup_heat_transfer
```

Each CQMOM case writes its existing example outputs plus
`output/moment_comparison.csv`, with the aligned columns `time`, `quantity`,
`qmomkit`, and `monte_carlo`. The 1D aerosol case preserves its existing
per-node QMOM and Monte Carlo moment-history outputs.

The 4D and 5D breakup cases each use 10,000 initialization samples with seed
`123`, breakup seed `321`, and end time `0.1`. The 4D case uses the `(d,u,v,w)`
phase space, and the 5D case uses `(V,T,u,v,w)`. The 5D case also writes
`output/generalized_moment_history_mc.txt` for direct-ensemble histories of
`M[4/3,0,0,0,0]` and `M[4/3,1,0,0,0]`; CSV identifiers use forms such as
`M_1over3_*` and `M_4over3_*`. Native raw-moment and quadrature histories are
written according to each `configuration.qmom` without normalization.

For reproducibility, the 4D diameter marginal is the configured trimodal
lognormal mixture with component number densities `0.40`, `0.50`, and `0.10`,
median diameters `25`, `100`, and `400` micrometres, and log-space standard
deviation `0.12`; independent-product initialization normalizes those
components before applying the global number density of one.

The paper quantities are:

- Aerosol evaporation `(d)`: raw moments and quadratures for 2, 4, 6, 8, and
  10 QMOM nodes, alongside the sampled Monte Carlo moment history.
- Bubble cavitation `(R, U)`: `M_00`, `M_10/M_00`, `M_01/M_00`, `M_11/M_00`,
  `M_20/M_00`, and `M_02/M_00`.
- Droplet condensation-nucleation `(r, w, T)`: `M_000`, `M_300/M_200`,
  `M_010/M_000`, and `M_001/M_000`; `r` is radius and `M_300/M_200` is plotted
  without an added factor of two.
- Particle breakup with recoil `(d, u, v, w)`: raw `M[0,0,0,0]`,
  `M[3,0,0,0]`, `M[1,0,0,0]`, `M[2,0,0,0]`, `M[0,1,0,0]`,
  `M[0,2,0,0]`, `M[2,1,0,0]`, `M[2,2,0,0]`, and `M[6,0,0,0]`.
- Particle breakup with heat transfer `(V, T, u, v, w)`: raw
  `M[0,0,0,0,0]`, `M[1,0,0,0,0]`, `M[0,1,0,0,0]`, `M[0,2,0,0,0]`,
  `M[1,1,0,0,0]`, `M[1/3,0,0,0,0]`, `M[4/3,0,0,0,0]`,
  `M[1/3,1,0,0,0]`, and `M[4/3,1,0,0,0]`.

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
