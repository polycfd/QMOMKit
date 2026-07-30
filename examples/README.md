# QMOMKit Examples

The examples are standalone applications that show how to use QMOMKit.
They are small programs and configuration files meant for orientation. Each
example is built like a separate project that links against an installed
QMOMKit package.

## Beginner Workflow

1. `qmom_explicit_moments/` if you already have raw moments.
2. `qmom_distribution_initialization/` if you want QMOMKit to generate
   moments from an analytical distribution.
3. `qmom_1d_minimal/` for the smallest complete QMOM transport loop.
4. `qmom_1d_continuous_breakage/` for configuration-driven QMOM transport
   with continuous and discontinuous physics.
5. `cqmom_independent_product_initialization/` for the smallest
   multidimensional initialization workflow.
6. `cqmom_transport/` for the canonical CQMOM transport workflow
   with custom velocity, source, and boundary contribution callbacks, as well as
   output.

`qmom_1d_minimal/`, `qmom_1d_continuous_breakage/`, and
`cqmom_transport/` are the primary examples for physics callbacks. The
callback requirements are documented in the [transport guide](../docs/user/transport.md). Each example
chooses the callbacks needed by its model; those combinations are illustrative,
not mandatory for every transport workflow.

## Advanced And Reference Workflows

- `qmom_mixture_initialization/`: optional QMOM initialization extension for
  mixture distributions.
- `realizability_check/`: example with explicit opt-in diagnostics for moment realizability.
- `cqmom_2d_adaptive_node_counts/`: 2D CQMOM transport example for adaptive branch-local node-count reduction, persistence, and controlled recovery.
- `qmom_1d_continuous_breakage_mc_initialization/`: 1D example with continuous and discontinuous physics, with initial raw moments estimated from a fixed-seed Monte Carlo sample and a companion direct Monte Carlo evolution executable.
- `qmom_1d_continuous_breakage_mc_local_rates/`: 1D example with the same Monte Carlo initialization and physics as the previous example, but with application-managed moment-space SSP-RK3 integration through the public local moment-rate APIs.
- `cqmom_2d_continuous_breakage_mc_initialization/`: 2D example with continuous and discontinuous physics, with initial raw moments estimated from a fixed-seed Monte Carlo sample and a companion direct Monte Carlo evolution executable.

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

Small QMOMKit-only `.qmom` examples are provided for:

- `qmom_explicit_moments.qmom`
- `qmom_mixture_absolute_density.qmom`
- `realizability_policy_nonnegative.qmom`

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
