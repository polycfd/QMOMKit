# QMOMKit

**QMOMKit is a C++ library implementing the quadrature method of moments to represent and evolve number-density functions in arbitrary-dimensional phase spaces**

QMOMKit reconstructs quadrature nodes and weights from raw moments and provides phase-space evolution utilities for moment and quadrature states. It leaves moment transport in physical space, the physical model, governing equations, physical-space mesh and discretization, PDE solution, and overall simulation workflow to the host application. QMOMKit is independent of any particular CFD, FVM, FEM, or DG implementation. The **core QMOMKit package** has no mandatory third-party numerical dependency; optional tools may have their own dependencies.

## Authors

- [Fabian Denner](https://www.polymtl.ca/expertises/en/denner-fabian)
- [Bruno Blais](https://www.polymtl.ca/expertises/en/blais-bruno)

## License

QMOMKit is released under the **Apache License 2.0**. See [LICENSE](LICENSE) for details.

## Capabilities

QMOMKit provides tools to:

- reconstruct one-dimensional QMOM and multidimensional recursive CQMOM quadratures;
- configure phase-space conditioning order, bounded support, and adaptive branch-specific CQMOM node counts;
- initialize moments from explicit moment definitions, analytical distributions, or Monte carlo sampling;
- evolve states through continuous phase-space node advection and discontinuous moment-source or boundary contributions;
- use library-managed operator splitting or application-managed stage and local-rate workflows;
- perform realizability checks and report transport diagnostics;
- write native results output;
- use the optional Monte Carlo companion for sampled initialization and to produce reference solutions.

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
├── tests/                            # Unit, regression, verification, and consumer checks
├── tools/
│   ├── moment_rate_element_operator/ # Downstream local-rate prototype
│   ├── monte_carlo/                  # Optional Monte Carlo companion
│   └── plotting/                     # Optional plotting scripts
├── CMakeLists.txt
├── LICENSE
└── README.md
```

The main CMake project builds and installs the core library. Core tests use the appropriate build-tree target; examples, optional tools, benchmarks, and package consumers are separate projects that use the installed package where designed to do so. Examples are not part of the test suite.

## Installation

QMOMKit requires CMake 3.16 or newer and a C++17 compiler. The core package has no mandatory third-party numerical dependency.

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

The standalone examples progress from basic QMOM and CQMOM reconstruction and initialization to transport and operator-splitting workflows, adaptive quadrature and conditioning diagnostics, multidimensional physical-model demonstrations, and matched Monte Carlo comparisons. See the [examples repository index](examples/README.md) and the [user examples guide](docs/user/examples.md).

- `tools/monte_carlo/` provides sampled initialization, ensemble evolution, and finite-sample comparison workflows; see the [Monte Carlo guide](docs/user/monte_carlo.md).
- `tools/plotting/` provides optional plotting helpers for QMOMKit-native histories; see the [plotting README](tools/plotting/README.md).
- `tools/moment_rate_element_operator/` demonstrates downstream composition of local moment-rate evaluation; see its [README](tools/moment_rate_element_operator/README.md).

## Testing

The repository contains:

- unit tests for local behavior and failure paths;
- component regression tests for important interactions;
- verification tests against analytical, manufactured, convergence, or reconstruction references;
- public-header checks and installed-package consumer checks;
- separate performance benchmarks.

Examples are not tests, and benchmarks measure performance rather than
correctness. The installed-package consumer and optional tool tests have their
own documented workflows.

To enable the core test suite:

```bash
cmake -S . -B build-tests \
    -DQMOMKIT_BUILD_TESTS=ON
cmake --build build-tests
ctest --test-dir build-tests --output-on-failure
```

See the [tests README](tests/README.md), [testing guide](docs/developer/testing.md),
and [benchmarks guide](docs/developer/benchmarks.md) for details.

## Documentation

- **Getting started:** [installation](docs/user/installation.md), [quick start](docs/user/quick_start.md), and the [user guide](docs/user/README.md).
- **Using QMOMKit:** [configuration](docs/user/configuration.md), [QMOM](docs/user/qmom.md), [CQMOM](docs/user/cqmom.md), [distributions](docs/user/distributions.md), [realizability](docs/user/realizability.md), and [transport](docs/user/transport.md).
- **Examples:** [repository example index](examples/README.md) and [user example guide](docs/user/examples.md).
- **Optional Monte Carlo tooling:** [Monte Carlo companion guide](docs/user/monte_carlo.md).
- **Development:** [developer guide](docs/developer/README.md), [architecture](docs/developer/architecture.md), [testing](docs/developer/testing.md), and [benchmarks](docs/developer/benchmarks.md).

## Contributing

Contributions to improve and extend QMOMKit are welcome. Before contributing, read the relevant documentation, preserve public interfaces and numerical conventions, add tests or verification cases at the right level, and update documentation when behavior changes.

QMOMKit aims to remain a focused scientific computing library with clear interfaces, minimal dependencies, and carefully tested numerical behavior.
