# QMOMKit

**QMOMKit is a C++ library implementing the quadrature method of moments to represent and evolve number-density functions in arbitrary-dimensional phase spaces**

QMOMKit reconstructs quadrature nodes and weights from raw moments, allowing applications to represent, evolve, and analyze number-density functions efficiently. It provides the numerical tools for moment transport and quadrature-based evaluation, while leaving the physical model, governing equations, and overall simulation workflow to the application. QMOMKit has **no** mandatory third-party numerical dependency.

## Capabilities

QMOMKit provides tools to:

- reconstruct one-dimensional quadratures using the **quadrature method of moments (QMOM)**;
- reconstruct multidimensional quadratures using the **conditional QMOM (CQMOM)**;
- evaluate derived generalized moments with finite real exponents from reconstructed QMOM or CQMOM quadrature;
- use adaptive CQMOM branch counts when different branches need different active node counts;
- initialize raw moments from supported distributions, mixtures, and independent products;
- check realizability when moments must be consistent with a chosen support;
- evolve moment-based representations through continuous node motion and source-like moment updates;
- provide user-defined physics callbacks such as `PhaseSpaceVelocity`, `MomentSource`, `CqmomMomentSource`, `MomentBoundaryContribution`, and `CqmomMomentBoundaryContribution`;
- write moment and quadrature histories with `MomentWriter` and `QuadratureWriter`;
- generate sample-based initial conditions and stochastic reference histories with the optional Monte Carlo companion tools;
- make quick plots of QMOMKit-native moment and quadrature histories;
- use tests and verification cases that check the numerical methods against
  known results.

## Repository Structure

```text
QMOMKit/
├── include/                  # Headers used by applications
├── src/                      # QMOMKit library implementation
├── docs/
│   ├── user/                 # User documentation
│   └── developer/            # Developer documentation
├── examples/                 # Example applications
├── tests/                    # Tests and verification cases
├── tools/
│   ├── monte_carlo/          # Monte Carlo companion package
│   └── plotting/             # Lightweight plotting helper scripts
├── CMakeLists.txt
└── README.md
```

The main CMake project builds and installs the QMOMKit library. Examples, tests, benchmarks, and companion tools are built separately and use the installed library, following the same workflow as external applications.

## Installation

QMOMKit uses CMake, requires a C++17 compiler and has no third-party numerical dependency.

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

The examples are small applications that show how to initialize moments, reconstruct QMOM and CQMOM quadratures, evolve moment-based representations, and write native output files. Start with the [examples guide](docs/user/examples.md).

QMOMKitMonteCarlo is an optional companion package under `tools/monte_carlo/`. It helps researchers draw finite sample ensembles, initialize QMOMKit from sample-derived moments, evolve sample ensembles, and generate stochastic reference results. See the [Monte Carlo guide](docs/user/monte_carlo.md) for more information.

The plotting helpers under `tools/plotting/` provide quick visual checks for QMOMKit-native `MomentWriter` and `QuadratureWriter` text files.

## Testing

QMOMKit includes tests that check the library at several levels:

- unit tests for local behavior and failure paths;
- component regression tests for important interactions;
- verification tests against analytical, manufactured, convergence, or reconstruction references;
- public-header checks for installed include behavior;
- an installed-package smoke test that confirms another CMake project can find,
  include, and link QMOMKit.

Optional standalone benchmarks and internal performance-instrumentation builds
support performance investigations without changing the core runtime API. See
the [testing guide](docs/developer/testing.md), [benchmarks guide](docs/developer/benchmarks.md),
and [tests README](tests/README.md) for details.

## Documentation

User documentation:

- [User guide](docs/user/README.md)
- [Installation](docs/user/installation.md)
- [Quick start](docs/user/quick_start.md)
- [Examples](docs/user/examples.md)
- [User glossary](docs/user/glossary.md)

Developer documentation:

- [Developer guide](docs/developer/README.md)
- [Architecture](docs/developer/architecture.md)
- [Testing](docs/developer/testing.md)
- [Developer glossary](docs/developer/glossary.md)

Tool documentation:

- [Monte Carlo companion](docs/user/monte_carlo.md)
- [Plotting helpers](tools/plotting/README.md)

## Typical Workflows

First define a QMOM or CQMOM case. You can load that setup from a `.qmom` file 
or construct the case directly in C++ when the surrounding application already
manages setup data.
 
A file-based QMOM reconstruction workflow is:

```text
load_qmom(...) -> QmomCase -> compute_qmom(...)
```

A file-based CQMOM reconstruction workflow is:

```text
load_cqmom(...) -> CqmomCase -> compute_cqmom(...)
```

Direct C++ setup enters the same reconstruction routines with `QmomCase` for
QMOM or `CqmomCase` for CQMOM.

After reconstruction, users may inspect quadrature, evaluate derived quantities,
write output, or optionally evolve the distribution. Transport is available for
both formulations: QMOM transport advances raw moment vectors, while normal
CQMOM transport advances a `CqmomState` created from the case.
Configuration files use the `.qmom` extension and describe phase-space
structure, initial moments or analytical initialization data, and numerical
options. They are a convenient reproducibility tool, not a requirement for using
QMOMKit. Physics enters through user-defined callbacks.

Output helpers write QMOMKit-native text tables. `MomentWriter` records native
integer-order raw moments, while `QuadratureWriter` records reconstructed nodes
and weights. Generalized moments with real exponents are evaluated separately
from valid quadrature results.

## Contributing

Contributions to improve and extend QMOMKit are welcome. Before contributing, read the relevant documentation, preserve public interfaces and numerical conventions, add tests or verification cases at the right level, and update documentation when behavior changes.

QMOMKit aims to remain a focused scientific computing library with clear interfaces, minimal dependencies, and carefully tested numerical behavior.

## Authors

- [Fabian Denner](https://www.polymtl.ca/expertises/en/denner-fabian)
- [Bruno Blais](https://www.polymtl.ca/expertises/en/blais-bruno)

## License

QMOMKit is released under the **Apache License 2.0**. See [LICENSE](LICENSE) for details.
