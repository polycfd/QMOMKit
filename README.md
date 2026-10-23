# QMOMKit

<a href="https://github.com/polycfd/QMOMKit/actions/workflows/release.yml"><img src="https://github.com/polycfd/QMOMKit/actions/workflows/release.yml/badge.svg" alt="CI-Release"></a>
<a href="https://github.com/polycfd/QMOMKit/actions/workflows/debug.yml"><img src="https://github.com/polycfd/QMOMKit/actions/workflows/debug.yml/badge.svg" alt="CI-Debug"></a>
<a href="https://github.com/polycfd/QMOMKit/actions/workflows/warnings_gcc.yml"><img src="https://github.com/polycfd/QMOMKit/actions/workflows/warnings_gcc.yml/badge.svg" alt="CI-Warnings-GCC"></a>
<a href="https://github.com/polycfd/QMOMKit/actions/workflows/clang-format.yml"><img src="https://github.com/polycfd/QMOMKit/actions/workflows/clang-format.yml/badge.svg" alt="CI-Indentation"></a>
<a href="https://github.com/polycfd/QMOMKit/actions/workflows/consumers.yml"><img src="https://github.com/polycfd/QMOMKit/actions/workflows/consumers.yml/badge.svg" alt="CI-Consumers"></a>
<a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache_2.0-blue.svg" alt="License: Apache 2.0"></a>

**QMOMKit is a C++ library implementing quadrature-based moment methods to represent and evolve number density functions in phase space.**

QMOMKit reconstructs quadratures, sets of weighted nodes, from the raw moments of a number density and evolves them in phase space: it solves the phase-space part of a population balance equation at one location. Transport in physical space, the mesh, the time loop and the coupling to a host code belong to the host application. Examples of how to couple QMOMKit to host application codes are maintained separately in [QMOMKit-Couplings](https://github.com/polycfd/QMOMKit-Couplings). The library has **no third-party dependencies**.

## Capabilities

QMOMKit provides tools to:

- reconstruct a quadrature from raw moments, with the quadrature method of moments (QMOM) [[McGraw, 1997](https://www.doi.org/10.1080/02786829708965471)] for one phase-space coordinate and the conditional quadrature method of moments (CQMOM) [[Yuan & Fox, 2011](https://www.doi.org/10.1016/j.jcp.2011.07.020)] for several;
- let the quadrature reconstruction algorithm choose the number of nodes adaptively and bound the quadrature abscissas to an interval; 
- choose the CQMOM conditioning order without changing the definition of the phase space;
- advance a state in phase space by operator splitting (Strang or Lie) of node advection and moment sources, including exchange across phase-space boundaries;
- check whether moments are realizable on a given support;
- rebuild states from moments that the host has transported in physical space;
- set up a case in C++ or from a `.qmom` configuration file, with moments given explicitly or generated from analytical distributions, mixtures and independent products;
- write moment and quadrature histories;
- compare with finite-sample Monte Carlo references from the optional companion package.

## A first example

Reconstructing a two-node quadrature from a raw moment sequence:

```cpp
#include <qmomkit/qmom.hpp>

#include <iostream>

int main()
{
    // Raw moments m0..m3 for two nodes:
    // weights {0.4, 0.6}, abscissas {-1.0, 2.0}.
    const qmomkit::MomentSet moments({1.0, 0.8, 2.8, 4.4});
    const std::size_t node_count = 2;

    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, node_count);
    if (!result.ok())
    {
        std::cerr << "QMOM failed: " << result.message << '\n';
        return 1;
    }

    for (std::size_t node = 0; node < result.weights.size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << result.weights[node]
                  << ", abscissa=" << result.abscissas[node] << '\n';
    }
}
```

[Getting started](docs/user/2_getting_started.md) builds on this to cover configuration and phase-space evolution.

## Installation

QMOMKit requires CMake 3.16 or newer and a C++17 compiler.

From the repository root:

```bash
cmake -S . -B build
cmake --build build
cmake --install build
```

Unless `CMAKE_INSTALL_PREFIX` is set, this installs into `install/` inside the repository rather than a system-wide location. The [examples](examples/README.md) assume this default.

Applications using QMOMKit can link to the installed CMake target:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE qmomkit::qmomkit)
```

Point `CMAKE_PREFIX_PATH` at the QMOMKit install prefix when configuring your application. See [Installation](docs/user/1_installation.md) and [Getting started](docs/user/2_getting_started.md) for details.

## Examples and tools

The standalone [examples](examples/README.md) progress from basic QMOM and CQMOM reconstruction and initialization to phase-space evolution and operator-splitting workflows, adaptive quadrature and conditioning diagnostics, demonstrations using multidimensional physical models, and matched Monte Carlo comparisons. See the [examples index](examples/README.md) and the [list of examples in the user guide](docs/user/8_going_further.md#examples).

Two companion tools live under `tools/`:
- `tools/monte_carlo/` samples populations, evolves the samples and estimates their moments, as a reference for QMOMKit; see its [README](tools/monte_carlo/README.md).
- `tools/plotting/` provides optional plotting helpers for QMOMKit-native histories; see the [plotting README](tools/plotting/README.md).

## Testing

The test suite consists of:

- unit tests for local behavior and failure paths;
- component regression tests for important interactions;
- verification tests against analytical, manufactured, convergence, or reconstruction references;
- public-header checks; the installed-package consumer check runs separately.

Build and run it in a _Debug build_, because the unit tests use `assert()`:

```bash
cmake -S . -B build-tests -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build-tests --parallel
ctest --test-dir build-tests --output-on-failure
```
See the [testing guide](docs/developer/testing.md) for details.

## Documentation

- **Getting started:** the [user guide](docs/user/README.md), [installation](docs/user/1_installation.md), and [getting started](docs/user/2_getting_started.md).
- **Using QMOMKit:** [concepts](docs/user/3_concepts.md), [moments and quadrature](docs/user/4_moments_and_quadrature.md), [evolving in time](docs/user/5_evolving_in_time.md), [coupling QMOMKit](docs/user/6_coupling_qmomkit.md), and the [reference](docs/user/7_reference.md).
- **Examples:** the [examples index](examples/README.md) and the [list in the user guide](docs/user/8_going_further.md#examples).
- **Optional Monte Carlo tooling:** the [Monte Carlo companion](tools/monte_carlo/README.md).
- **Development:** the [developer guide](docs/developer/README.md), [architecture](docs/developer/architecture.md), and [testing](docs/developer/testing.md).
- **Examples coupling QMOMKit to host codes:** [QMOMKit-Couplings](https://github.com/polycfd/QMOMKit-Couplings)

## Repository structure

```text
QMOMKit/
├── .github/workflows/                # Continuous integration
├── include/qmomkit/                  # Installed public C++ API
├── src/                              # Core implementation
├── cmake/                            # CMake package templates
├── docs/
│   ├── developer/                    # Developer-oriented documentation
│   ├── doxygen/                      # Doxygen documentation
│   └── user/                         # User-oriented documentation
├── examples/                         # Standalone installed-package examples
├── paper/                            # QMOMKit paper, its examples and verification
├── tests/                            # Unit, regression, and verification tests
├── tools/
│   ├── monte_carlo/                  # Monte Carlo companion code
│   ├── plotting/                     # Plotting scripts
│   └── indent-all, ...               # Formatting helpers used by CI
├── AGENTS.md                         # Instructions for AI coding agents
├── CMakeLists.txt
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── DEVELOPMENT_PHILOSOPHY.md         # How QMOMKit is developed and verified
├── LICENSE
└── README.md
```

## Development philosophy and verification

Generative AI tools have been used extensively during the development of QMOMKit. Their use has reinforced a central principle of the project: confidence in scientific software should come from evidence rather than from the origin of an implementation. Neither human-written nor AI-assisted code is inherently correct.

QMOMKit therefore relies on a deliberately rigorous approach to testing, numerical verification, and documentation, so that its behavior can be checked rather than taken on trust. See [DEVELOPMENT_PHILOSOPHY.md](DEVELOPMENT_PHILOSOPHY.md) for the full statement and for where that evidence lives.

## How to cite

If QMOMKit contributes to your academic work, please cite the QMOMKit paper:

Fabian Denner, Bruna Campos, and Bruno Blais, *QMOMKit: A C++ library for quadrature-based moment methods*, arXiv:2610.?????, 2026. doi: [10.48550/arXiv.2610.?????](https://www.doi.org/10.48550/arXiv.2610.?????)

```bibtex
@misc{qmomkit,
  author  = {Denner, Fabian and Campos, Bruna and Blais, Bruno},
  title   = {QMOMKit: A C++ library for quadrature-based moment methods},
  eprint  = {2610.?????},
  archivePrefix = {arXiv},
  year    = {2026},
  doi     = {10.48550/arXiv.2610.?????},
}
```

The studies behind the paper can be reproduced from this repository; see the [paper studies README](paper/README.md).

## Authors

- [Fabian Denner](https://www.polymtl.ca/expertises/en/denner-fabian)
- [Bruno Blais](https://www.polymtl.ca/expertises/en/blais-bruno)
- Bruna Campos

## License

Copyright 2026 The QMOMKit developers.

QMOMKit is released under the Apache License 2.0. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for
development, testing, and review.
