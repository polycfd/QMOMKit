# CQMOM: 5D Non-Spherical Particles

## Purpose

This example demonstrates a five-dimensional CQMOM phase space for
non-spherical particles:

```text
xi = (d, chi, u, v, w)
```

where `d` is equivalent particle diameter, `chi` is an aspect-ratio shape
parameter, and `u`, `v`, and `w` are particle velocity components.

The example provides two standalone executables:

- `cqmom_non_spherical_particles_example`: samples the configured initial
  distribution, initializes CQMOM from the sampled raw moments, and advances the
  CQMOM transport state.
- `cqmom_non_spherical_particles_mc_example`: samples the same configured
  initial distribution and evolves the finite ensemble directly.

Both drivers use the same configuration file, sample count, initialization seed,
phase-space ordering, continuous particle RHS, and no-op breakage stage. The
current time grids differ, so the generated histories are not a validated
one-to-one quantitative comparison.

## Particle Model

The continuous particle model lives in `particle_model.hpp` and
`particle_model.cpp`. The configured phase-space order is `(d, chi, u, v, w)`,
and the velocity callback interprets the supplied coordinate vector directly in
that order.

The RHS evaluates:

```text
dd/dt   = -attrition_rate * |U_f(t) - u_p|^2
dchi/dt = (1 - chi) / aspect_ratio_relaxation_time
du/dt   = drag_rate * (U_f,x(t) - u)
dv/dt   = drag_rate * (U_f,y(t) - v)
dw/dt   = drag_rate * (U_f,z(t) - w)
```

where the drag rate uses the oblate-spheroid shape correction implemented in
`particle_model.cpp`. The current code validates a finite `gravity` parameter
but does not include gravity in `dw/dt`; the gravity term is intentionally not
active in this example.

`particle_model_parameters()` currently supplies finite model parameters in
code. The configuration file controls the CQMOM hierarchy, initialization
distribution, transport options, and output options; it does not define the
particle-model constants.

## Breakage Stage

The breakage extension point lives in `breakage_model.hpp` and
`breakage_model.cpp`.

The current CQMOM source callback returns zero moment rates. The current direct
MC breakage stage validates that samples are five-dimensional and otherwise
leaves the ensemble unchanged. The macro-step sequence is still Strang
splitting:

```text
continuous half-step
discontinuous full step
continuous half-step
```

For CQMOM, `qmomkit::advance_transport_cqmom(...)` owns this sequence. For the
direct MC executable, the two continuous half-steps are advanced with
`qmomkit::monte_carlo::evolve_ensemble(...)`, which uses RK23 for sample
evolution.

## Initialization And Configuration

Edit `configuration.qmom` to change the independent-product distribution and
CQMOM hierarchy. The phase-space sections are ordered by coordinate:

```text
0: d
1: chi
2: u
3: v
4: w
```

The default hierarchy is:

```text
node_counts = 2, 2, 2, 2, 2
conditioning_order = 0, 1, 2, 3, 4
adaptive_branch_counts = true
```

The initialization uses a fixed-seed Monte Carlo sample to estimate a dense raw
moment tensor:

```text
M[i_d, i_chi, i_u, i_v, i_w]
```

The tensor extents are controlled by:

```text
maximum_orders = 3, 3, 3, 3, 3
```

which matches the moment orders needed for the configured two-node CQMOM
hierarchy. Monte Carlo sampled moments use QMOMKit's absolute number-density
convention:

```text
M[alpha] = n0 / N * sum_p product_r xi_(r,p)^alpha_r
```

Adaptive branch counts are enabled by `phase_space.adaptive_branch_counts`.
QMOMKit automatically tracks active branch-local node counts during CQMOM
transport; source node-count evolution is not a user-configured policy in the
current API.

Optional transport settings may be provided through `[time_integration]` and
`[time_integration.node_advection]`. The CQMOM executable loads them with
`qmomkit::config::load_cqmom_transport_options(...)` and reports the transport
scheme and effective node-advection settings at startup.

## Time Controls

The current executables use different time controls:

```text
CQMOM: final_time = 100.0, time_step = 0.005
MC:    final_time = 1.0,   time_step = 0.01
```

This difference is preserved intentionally in the current example state. Use
the outputs as workflow demonstrations unless a later design decision selects a
common time horizon and macro step for quantitative comparison.

## Build

Requires:

- QMOMKit
- QMOMKitMonteCarlo

From the repository root, install both packages into the shared repository-local
`install/` prefix:

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
```

Then configure this example from its directory:

```bash
cd examples/cqmom_non-spherical-particles

cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`.

## Run

Run from the build directory so generated output stays there:

```bash
cd build
./cqmom_non_spherical_particles_example
./cqmom_non_spherical_particles_mc_example
```

The CQMOM executable writes:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

The direct MC executable writes:

- `output/moment_history_mc.txt`
