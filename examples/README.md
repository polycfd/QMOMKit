# QMOMKit Examples

The examples are standalone applications that show how to use QMOMKit.
They are small programs and configuration files meant for orientation. Each
example is built like a separate project that links against an installed
QMOMKit package.

The [user guide](../docs/user/8_going_further.md#examples) lists the
examples in order of increasing scope, with what each one shows. The examples
with a Monte Carlo reference require the QMOMKitMonteCarlo companion; their
particle executables are stochastic references, not exact solutions or tests.

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

Install QMOMKit first, as described in [Installation](../docs/user/1_installation.md).
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
`cqmom_2d_adaptive_node_counts/configuration.qmom`. The
[Reference](../docs/user/7_reference.md#configuration-file) lists every key.

## Examples versus tests

Examples are intended to be readable and short. They are built with the normal
CMake installed-package workflow but are not registered as CTest tests.

Numerical verification lives under `tests/`. See the associated
[README](../tests/README.md) for the test categories and how to run them.
