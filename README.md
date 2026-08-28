# QMOMKit

[![CI-Release](https://github.com/polycfd/QMOMKit/actions/workflows/release.yml/badge.svg)](https://github.com/polycfd/QMOMKit/actions/workflows/release.yml)
[![CI-Debug](https://github.com/polycfd/QMOMKit/actions/workflows/debug.yml/badge.svg)](https://github.com/polycfd/QMOMKit/actions/workflows/debug.yml)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

**QMOMKit is a C++ library implementing the quadrature method of moments to represent and evolve number-density functions in arbitrary-dimensional phase spaces**

QMOMKit reconstructs quadrature nodes and weights from raw moments and utilities for the evolution of moments and the corresponding quadratures in phase space. It is independent of any physical application or host code coupling and, hence, leaves transport of moments in the physical space, the representation and discretization of the physical space and overall simulation workflow to the host application. The QMOMKit library has **no third-party dependencies**.

## Authors

- [Fabian Denner](https://www.polymtl.ca/expertises/en/denner-fabian)
- [Bruno Blais](https://www.polymtl.ca/expertises/en/blais-bruno)

## License

QMOMKit is released under the **Apache License 2.0**. See [LICENSE](LICENSE) for details.

## Capabilities

QMOMKit provides tools to:

- reconstruct one-dimensional QMOM and multidimensional CQMOM quadratures;
- configure phase-space conditioning order, bounded support, and adaptive branch-specific CQMOM node counts;
- initialize moments from explicit moment definitions, analytical distributions, or Monte carlo sampling;
- evolve states through phase-space advection and moment-source or boundary contributions;
- use library-managed operator splitting or application-managed stage and local-rate workflows;
- perform realizability checks and report phase-space evolution diagnostics;
- write native results output;
- use the optional Monte Carlo companion code for initialization of population ensembles and to produce reference solutions.

## Repository Structure

```text
QMOMKit/
├── include/qmomkit/                  # Installed public C++ API
├── src/                              # Core implementation
├── cmake/                            # CMake package templates
├── docs/
│   ├── user/                         # User-oriented documentation
│   └── developer/                    # Developer-oriented documentation
├── examples/                         # Standalone installed-package examples
├── paper/                            # QMOMKit paper, its examples and verification
├── tests/                            # Unit, regression, and verification tests
├── tools/
│   ├── galerkin_element_coupling/    # Galerkin coupling demonstrator
│   ├── monte_carlo/                  # Monte Carlo companion code
│   └── plotting/                     # Plotting scripts
├── CMakeLists.txt
├── LICENSE
└── README.md
```


## Installation

QMOMKit requires CMake 3.16 or newer and a C++17 compiler. 

From the repository root:

```bash
cmake -S . -B build
cmake --build build
cmake --install build
```

Applications using QMOMKit can link to the installed CMake target:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE qmomkit::qmomkit)
```

Point `CMAKE_PREFIX_PATH` at the QMOMKit install prefix when configuring your application. See the [installation guide](docs/user/installation.md) and [quick start](docs/user/quick_start.md) for details.

## Examples and Tools

The standalone examples progress from basic QMOM and CQMOM reconstruction and initialization to phase-space evolution and operator-splitting workflows, adaptive quadrature and conditioning diagnostics, multidimensional physical-model demonstrations, and matched Monte Carlo comparisons. See the [examples repository index](examples/README.md) and the [user examples guide](docs/user/examples.md).

In addition, QMOMKit has several companion tools:
- `tools/monte_carlo/` provides Monte Carlo code that can be use to initialize and evolve population ensembles; see the [Monte Carlo guide](docs/user/monte_carlo.md).
- `tools/plotting/` provides optional plotting helpers for QMOMKit-native histories; see the [plotting README](tools/plotting/README.md).
- `tools/galerkin_element_coupling/` demonstrates Galerkin coupling to local QMOMKit phase-space operations; see its [README](tools/galerkin_element_coupling/README.md).

## Testing

The repository contains **over 100 individual tests**:

- unit tests for local behavior and failure paths;
- component regression tests for important interactions;
- verification tests against analytical, manufactured, convergence, or reconstruction references;
- public-header checks and installed-package consumer checks;
- separate performance benchmarks.

To enable the core test suite, use a _Debug build_ because the unit tests use
`assert()`:

```bash
cmake -S . -B build-tests \
    -DCMAKE_BUILD_TYPE=Debug \
    -DQMOMKIT_BUILD_TESTS=ON
cmake --build build-tests --parallel
ctest --test-dir build-tests --output-on-failure
```
See the [tests README](tests/README.md), [testing guide](docs/developer/testing.md),
and [benchmarks guide](docs/developer/benchmarks.md) for details.

## Documentation

- **Getting started:** [installation](docs/user/installation.md), [quick start](docs/user/quick_start.md), and the [user guide](docs/user/README.md).
- **Using QMOMKit:** [configuration](docs/user/configuration.md), [QMOM](docs/user/qmom.md), [CQMOM](docs/user/cqmom.md), [distributions](docs/user/distributions.md), [realizability](docs/user/realizability.md), and [phase-space evolution](docs/user/phase_space_evolution.md).
- **Examples:** [repository example index](examples/README.md) and [user example guide](docs/user/examples.md).
- **Optional Monte Carlo tooling:** [Monte Carlo companion guide](docs/user/monte_carlo.md).
- **Development:** [developer guide](docs/developer/README.md), [architecture](docs/developer/architecture.md), [testing](docs/developer/testing.md), and [benchmarks](docs/developer/benchmarks.md).

## Contributing

Contributions to improve and extend QMOMKit are very welcome. Before contributing, read the relevant documentation, preserve public interfaces and numerical conventions, add tests or verification cases for all new features and functionalities, and update the documentation where appropriate.

QMOMKit aims to remain a focused scientific computing library with clear interfaces, minimal dependencies, and carefully tested numerical behavior.
