# Installing QMOMKit

This page describes the repository-local installation workflow for QMOMKit.
The root project builds and installs the QMOMKit CMake package. 

Normal tests
can be enabled in the root build so they validate the current source tree.
Examples, benchmarks, the installed-package consumer, and the Monte Carlo
companion package are standalone projects that use an installed package.

## Requirements

QMOMKit has **no** mandatory third-party numerical dependency.
QMOMKit only requires a C++17 compiler and CMake 3.16 or newer. CTest is needed
to build and run the test suite. 

## Build And Install QMOMKit

From the repository root, configure, build, and install the core package:

```bash
cmake -S . -B build 
cmake --build build
cmake --install build
```

This writes to `install/` inside the repository:

```text
install/
├── include/
├── lib/
└── lib/cmake/QMOMKit/
```

The installed CMake package exports the `qmomkit::qmomkit` target. Application codes
using QMOMKit should point `CMAKE_PREFIX_PATH` at this `install/` directory and 
would typically use the installed package through CMake:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE qmomkit::qmomkit)
```

For example, from a separate CMake application:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/QMOMKit/install
cmake --build build
```

## Build Tests

Normal unit, regression, verification, and public-header tests link the current
build-tree QMOMKit target. Installation is not required:

```bash
cmake -S . -B build-tests -DQMOMKIT_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-tests
ctest --test-dir build-tests --output-on-failure
```

For more detailed test selection, such as running only unit or verification
tests, see [tests/README.md](../../tests/README.md).

## Build Examples

Examples are also standalone CMake projects. Build a standard example
from its own directory after installing QMOMKit:

```bash
cd examples/qmom_explicit_moments
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_explicit_moments_example
```

Each example has its own README with the exact executable name and any local
configuration files it expects. Examples with a local `configuration.qmom`
copy it into their `build/` directory during CMake configuration.

## Optional QMOMKitMonteCarlo Package

`QMOMKitMonteCarlo` is an optional companion package under
`tools/monte_carlo/`. It provides finite ensemble sampling, sample-derived raw
moments, advection sample evolution, and moment-history comparison utilities.
It is not needed for ordinary QMOM or CQMOM use.

Build and install QMOMKit first, using the commands above. Then configure the
companion package against the installed QMOMKit package:

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

This builds the Monte Carlo library and the user-facing
`monte_carlo_sample_based_initialization_example` executable, then installs
the companion package into the same repository-local `install/` prefix. It
does not run tests.

Run Monte Carlo tests separately when developing or validating the package:

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build-tests \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DBUILD_TESTING=ON

cmake --build tools/monte_carlo/build-tests
ctest --test-dir tools/monte_carlo/build-tests --output-on-failure
```

An example or separate project that uses the companion package can point
`CMAKE_PREFIX_PATH` at the shared prefix:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="/path/to/QMOMKit/install"
```

For more detail, see [Monte Carlo Companion Package](monte_carlo.md).