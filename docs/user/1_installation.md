# Installation

QMOMKit does not have any third-party dependencies. It needs only a C++17
compiler and CMake 3.16 or newer.

## Build and install

From the repository root:

```bash
cmake -S . -B build
cmake --build build
cmake --install build
```

With a single-configuration generator this builds `Release`; it installs into
`install/` in the repository. Choose otherwise with `-DCMAKE_BUILD_TYPE=...` and
`-DCMAKE_INSTALL_PREFIX=...`.

## Use the installed package

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE qmomkit::qmomkit)
```

Configure your project with `-DCMAKE_PREFIX_PATH=/path/to/QMOMKit/install`, and
include `<qmomkit/qmomkit.hpp>` in your code.

## Build QMOMKit inside your project

Instead of installing it, you can add QMOMKit to your own build with
`add_subdirectory(...)` or `FetchContent` and link `qmomkit::qmomkit`. QMOMKit's
defaults then apply to your whole project:

- without a build type, a single-configuration build becomes `Release`;
- without an install prefix, the prefix becomes `install/` in your project's
  source directory, even if your project set a default prefix before adding
  QMOMKit;
- installing your project also installs QMOMKit.

Values from the command line or a preset are kept. Set a default prefix of
your own after adding QMOMKit:

```cmake
add_subdirectory(external/QMOMKit)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install path prefix" FORCE)
endif()
```

## Options

| Option | Default | Effect |
| --- | --- | --- |
| `QMOMKIT_BUILD_TESTS` | `OFF` | Builds the test suite; see [tests/README.md](../../tests/README.md). |
| `QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION` | `OFF` | Counts numerical work, for development; see the [developer documentation](../developer/performance_counters.md). |

## Examples and the Monte Carlo companion

Each example under [`examples/`](../../examples/) is a separate project built
against the installed package, for example:

```bash
cd examples/qmom_explicit_moments
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
```

The optional Monte Carlo companion package has its own build instructions in
[tools/monte_carlo/README.md](../../tools/monte_carlo/README.md).
