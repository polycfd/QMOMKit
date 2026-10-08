# QMOMKit Examples

The examples are standalone applications that show how to use QMOMKit.
They are small programs and configuration files meant for orientation. Each
example is built like a separate project that links against an installed
QMOMKit package.

## Getting started

For a short QMOM path, start with:

1. `qmom_explicit_moments/` to reconstruct a quadrature from application-owned
   raw moments.
2. `qmom_distribution_initialization/` to generate moments from an analytical
   distribution before reconstruction.
3. `qmom_1d_minimal/` for the smallest complete QMOM phase-space evolution
   workflow.

For the corresponding CQMOM path, continue with:

1. `cqmom_independent_product_initialization/` for the normal multidimensional
   configuration and initialization lifecycle.
2. `cqmom_2d_minimal/` for canonical CQMOM phase-space evolution with velocity,
   source, and boundary-contribution callbacks.

The [user example guide](../docs/user/examples.md) explains the learning path,
optional branches, and prerequisites in more detail.

## Initialization and diagnostics

- `qmom_mixture_initialization/`: generate QMOM moments from an
  absolute-density mixture when a single analytical distribution is not enough.
- `realizability_check/`: explicitly classify a raw 1D moment sequence under a
  declared support assumption before attempting reconstruction.

## QMOM evolution

- `qmom_1d_aerosol_evaporation/`: configuration-driven six-node, node-advection
  workflow without a moment source.
- `qmom_1d_breakage_advection/`: normal library-managed QMOM evolution with
  node advection, a size-dependent fragmentation source written over the
  quadrature nodes, and configured Strang splitting.
- `qmom_1d_explicit_strang_splitting/`: advanced application-owned Strang
  staging through the public advection and moment-source stage APIs; use after
  the library-managed workflow.

## CQMOM evolution

- `cqmom_2d_batch_host_data/`: many CQMOM states advanced with one
  batch-index-aware call per step, with callbacks reading time-varying host
  data indexed like the state vector.

## CQMOM numerical controls

- `cqmom_2d_adaptive_node_counts/`: adaptive conditional branch reduction,
  persistence, and controlled recovery during phase-space evolution.

## Monte Carlo companion / reference workflows

These examples require QMOMKitMonteCarlo. Their direct particle executables
are stochastic companion/reference realizations for the same sampled
initialization and physical model; they are not automatic formal validation
tests or exact solutions.

- `qmom_1d_breakage_advection_with-mc/`: sampled 1D initialization with normal
  QMOM evolution and a direct particle companion.
- `cqmom_2d_breakage_advection_with-mc/`: sampled 2D CQMOM initialization,
  evolution, and direct particle companion.
- `cqmom_2d_bubbles_cavitation/`: nonlinear 2D bubble-cavitation node advection
  with a direct sampled-population companion.
- `cqmom_3d_droplets_condensation-nucleation/`: 3D condensation and finite-time
  nucleation boundary birth with a direct droplet-population companion.

## Advanced physical and source-closure examples

- `cqmom_4d_particles_breakage_recoil/`: terminal-quadrature traversal closes a
  nonlinear daughter-state transformation with diameter-dependent recoil and a
  direct MC companion.
- `cqmom_5d_particles_breakage_heat-transfer/`: generalized fractional moments
  close a moment-space breakup source with daughter-temperature redistribution
  and a direct MC companion.
- `cqmom_5d_non-spherical-particles/`: 5D recursive CQMOM with bounded
  aspect-ratio support, adaptive branches, a multimodal velocity marginal, a
  paired bounded/unbounded reconstruction diagnostic, and a host-data pattern
  for time-varying callback inputs, plus a direct MC companion.

## Package dependencies

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

## Build pattern

Install QMOMKit first, as described in the [installation guide](../docs/user/installation.md).
From an example directory, configure against that install prefix:

```bash
cd examples/example_name
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/qmomkit
cmake --build build
cd build
./example_executable
```

Examples using the Monte Carlo companion to produce reference results
need `QMOMKitMonteCarlo` installed to the same prefix;
see its [README](../tools/monte_carlo/README.md).

Most examples with a local `configuration.qmom` copy it into `build/` during
CMake configuration, so generated output is written under the build directory.
A few instead take the configuration file as a command-line argument and are run
from the example directory; their own README says which.

Examples that write QMOM or CQMOM moments and quadrature use
`qmomkit::MomentWriter` and `qmomkit::QuadratureWriter` to write `.txt` files. These
files are whitespace-delimited text tables with commented header lines.
Monte Carlo particle snapshots, when present, are custom sample-level output
and are not QMOMKit moment or quadrature files.

## Configuration files

Configuration files are kept with the executable workflows that use them. For
example, see `qmom_mixture_initialization/configuration.qmom`,
`qmom_1d_breakage_advection/configuration.qmom`, and
`cqmom_2d_adaptive_node_counts/configuration.qmom`. The user configuration guide also
includes inline syntax for explicit moments and runtime realizability settings.

Monte Carlo companion configuration examples live with the companion package
under `tools/monte_carlo/configuration/`:

- `qmom_distribution_lognormal.qmom`
- `cqmom_independent_product_2d.qmom`

## Examples versus tests

Examples are intended to be readable and short. They are built with the normal
CMake installed-package workflow but are not registered as CTest tests.

Numerical verification lives under `tests/`. See the associated
[README](../tests/README.md) for the test categories and how to run them.
Performance benchmarks live under `tests/benchmarks/`. They are timing
measurements, not examples or correctness tests.
