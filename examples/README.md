# QMOMKit Examples

The examples are standalone applications that show how to use QMOMKit.
They are small programs and configuration files meant for orientation. Each
example is built like a separate project that links against an installed
QMOMKit package.

## Beginner Workflow

1. `qmom_explicit_moments/` if you already have raw moments.
2. `qmom_distribution_initialization/` if you want QMOMKit to generate
   moments from an analytical distribution.
3. `qmom_1d_minimal/` for the smallest complete QMOM phase-space evolution loop.
4. `qmom_1d_aerosol_evaporation/` for a configuration-driven six-node QMOM
   aerosol evaporation workflow using node advection.
5. `qmom_1d_explicit_strang_splitting/` for application-managed adaptive QMOM
   Strang splitting with the public advection and moment-source-evolution stage API.
6. `qmom_1d_breakage_advection/` for configuration-driven QMOM phase-space evolution
   with advection and moment-source physics.
7. `cqmom_independent_product_initialization/` for the smallest
   multidimensional initialization workflow.
8. `cqmom_2d_minimal/` for the canonical CQMOM phase-space evolution workflow
   with custom velocity, source, and boundary contribution callbacks, as well as
   output.

`qmom_1d_minimal/`, `qmom_1d_aerosol_evaporation/`, `qmom_1d_explicit_strang_splitting/`,
`qmom_1d_breakage_advection/`, and
`cqmom_2d_minimal/` are the primary examples for physics callbacks. The
callback requirements are documented in the [phase-space evolution guide](../docs/user/phase_space_evolution.md). Each example
chooses the callbacks needed by its model; those combinations are illustrative,
not mandatory for every phase-space evolution workflow.

## Advanced And Reference Workflows

- `qmom_mixture_initialization/`: optional QMOM initialization extension for
  mixture distributions.
- `realizability_check/`: example with explicit opt-in diagnostics for moment realizability.
- `qmom_1d_explicit_strang_splitting/`: adaptive 1D QMOM phase-space evolution example with application-managed Strang stages.
- `cqmom_2d_adaptive_node_counts/`: 2D CQMOM phase-space evolution example for adaptive branch-local node-count reduction, persistence, and controlled recovery.
- `cqmom_2d_borderline_conditioning/`: static 2D CQMOM reconstruction sweep for conditioning-order sensitivity, borderline Vandermonde systems, solver comparisons, failed-stage interpretation, lower-node reconstructions, and human-readable plus flat diagnostic reports.
- `qmom_1d_breakage_advection_with-mc/`: 1D example with advection and moment-source physics, with initial raw moments estimated from a fixed-seed Monte Carlo sample and a companion direct Monte Carlo evolution executable.
- `qmom_1d_breakage_advection_local-rates_with-mc/`: 1D example with the same Monte Carlo initialization and physics as the previous example, but with application-managed moment-space SSP-RK3 integration through the public local moment-rate APIs.
- `cqmom_2d_breakage_advection_with-mc/`: 2D example with advection and moment-source physics, with initial raw moments estimated from a fixed-seed Monte Carlo sample and a companion direct Monte Carlo evolution executable.
- `cqmom_3d_droplets_condensation-nucleation/`: 3D CQMOM phase-space evolution example for a condensing droplet population, with matched Monte Carlo initialization and a companion direct Monte Carlo evolution executable. Builds `cqmom_3d_droplets_condensation_nucleation_example` and `cqmom_3d_droplets_condensation_nucleation_mc_example`.
- `cqmom_2d_bubbles_cavitation/`: 2D CQMOM phase-space evolution example for a polydisperse cavitation population, with matched Monte Carlo initialization and a companion direct Monte Carlo evolution executable.
- `cqmom_5d_non-spherical-particles/`: 5D CQMOM example for nonspherical particles with bounded aspect-ratio support, adaptive branch counts, moment-source fragmentation, and a companion Monte Carlo workflow.
- `cqmom_4d_particles_breakage_recoil/`: 4D spherical-particle `(diameter, velocity components)` example with constant-rate binary breakup and diameter-dependent recoil. It evaluates the nonlinear daughter-state transformation from source-stage terminal quadrature.
- `cqmom_5d_particles_breakage_heat-transfer/`: 5D spherical-particle `(volume, temperature, velocity components)` example with constant-rate binary breakup and volume-dependent daughter-temperature redistribution. It evaluates a moment-space source closed with integer and generalized moments from `CqmomSourceContext`.

## Package Dependencies

Standard examples consume the installed QMOMKit package:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(example_target PRIVATE qmomkit::qmomkit)
```

Examples that use Monte Carlo initialization or direct Monte Carlo evolution
also require the QMOMKitMonteCarlo companion package:

```cmake
find_package(QMOMKitMonteCarlo CONFIG REQUIRED)
target_link_libraries(example_target PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)
```

## Build Pattern

From the repository root, install QMOMKit into the repository-local `install/`
prefix:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build build
cmake --install build
```

For Monte Carlo examples, install the companion package into the same prefix:

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

From an individual standard QMOMKit example directory, configure against an
installed QMOMKit package:

```bash
cd examples/example_name

cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

Examples using the Monte Carlo package use the same prefix:

```bash
cd examples/example_name

cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

Run from the example build directory:

```bash
cd build
./example_executable
```

Examples with a local `configuration.qmom` copy it into `build/` during CMake
configuration, so generated output is written under the build directory.

Examples that write QMOM or CQMOM moments and quadrature use
`qmomkit::MomentWriter` and `qmomkit::QuadratureWriter` to write `.txt` files. These
files are whitespace-delimited text tables with commented header lines.
Monte Carlo particle snapshots, when present, are custom sample-level output
and are not QMOMKit moment or quadrature files.

## Configuration Files

Configuration files are kept with the executable workflows that use them. For
example, see `qmom_mixture_initialization/configuration.qmom`,
`qmom_1d_breakage_advection/configuration.qmom`, and
`cqmom_2d_minimal/configuration.qmom`. The user configuration guide also
includes inline syntax for explicit moments and runtime realizability settings.

Monte Carlo validation configuration examples live with the companion package
under `tools/monte_carlo/configuration/`:

- `qmom_distribution_lognormal.qmom`
- `cqmom_independent_product_2d.qmom`

## Examples Versus Tests

Examples are intended to be readable and short. They are built with the normal
CMake installed-package workflow but are not registered as CTest tests.

Numerical verification lives under `tests/`. See the associated
[README](../tests/README.md) for the test categories and how to run them.
Performance benchmarks live under `tests/benchmarks/`. They are timing
measurements, not examples or correctness tests.
