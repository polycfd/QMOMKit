# QMOMKit Tests

QMOMKit uses plain C++ test executables registered with CTest. The suite is
organized by purpose through directory placement and CTest labels.

## Testing Philosophy

```text
QMOMKit testing

                    Test purpose
                         |
        -----------------------------------
        |                 |               |
      unit        component_regression   verification
        |                 |               |
 isolated correctness   implementation   mathematical and
 local APIs             stability        workflow correctness
```

The categories are separate because they answer different questions:
- Unit tests check whether one function, type, or small component behaves as expected.
- Component regression tests check whether important implementation behavior
remains stable. 
- Verification tests check whether QMOMKit reproduces a mathematical
property or public workflow against an independent reference.

When adding a new test, use this decision rule:

1. Does it verify an isolated function, class, or local API behavior? Put it in
   `tests/unit/`.
2. Does it protect specific implementation behavior or a regression-prone
   component interaction? Put it in `tests/component_regression/`.
3. Does it verify a mathematical property or user workflow against an
   analytical, manufactured, convergence, or otherwise independent reference?
   Put it in `tests/verification/`.

## Test Categories

### Unit

Unit tests verify isolated functions, classes, and local API behavior. They are
focused, fast, and cover edge cases, preferably in ways that are independent of
implementation details where possible. They are the right place for:

- option and request validation;
- status helpers and diagnostics on local APIs;
- indexing and storage conventions;
- failure handling and state preservation;
- small component behavior that does not need a workflow-level reference.

Representative examples include indexing helpers, containers, parser helpers,
and local numerical utilities.

Unit tests live under `tests/unit/` and use the `unit` CTest label. Their CTest
names use the `unit_<area>` prefix.

### Component Regression

Component regression tests protect important implementation behavior that is
larger than an isolated unit check but narrower than mathematical verification.
They cover fragile or historically important behavior such as integrator
composition, direct component overload behavior, diagnostic/status propagation,
policy mapping, and previously fragile implementation paths.

These tests are not failed verification tests. They may use manufactured inputs
or simple analytical expectations, but their primary purpose is regression
protection rather than executable documentation of a mathematical workflow. They
are also not a replacement for unit tests.

Component regression tests live under `tests/component_regression/` and use the
`component_regression` CTest label. Their CTest names use the
`component_regression_<area>` prefix.

### Verification

Verification tests verify mathematical correctness and public workflows. They
compare against analytical references, manufactured solutions, convergence
studies, mathematical invariants, or other independently known results. They
serve as executable documentation for intended behavior.

Representative examples include QMOM/CQMOM reconstruction, phase-space evolution workflows,
physics callbacks, and output workflows.

Verification tests live under `tests/verification/` and use the `verification`
CTest label. Their CTest names use the `verification_<area>_<case>` prefix.
Each verification area has a local README when the workflow or reference
construction is non-trivial.

Verification logic must not rely on C `assert`; checks must execute in Debug
and Release builds and return a non-zero exit code on failure.

### Public Headers

Public-header checks verify at compile time that every public header is
self-contained. Each `include_*.cpp` translation unit under
`tests/public_headers/` includes one QMOMKit public header through the current
source-tree `qmomkit` target. This detects missing direct includes and
include-order dependencies before installation.

`tests/public_headers/main.cpp` is intentionally empty apart from `main()`. It
only permits the compile-check translation units to be linked into a CTest
executable. The executable performs no numerical runtime test; compilation is
the meaningful check.

### Installed Package Consumer

The installed-package smoke test under `tests/package_consumer/` is a minimal
external CMake project. It validates the installed package interface, not the
in-tree build target, by using:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(consumer PRIVATE qmomkit::qmomkit)
```

This check is intentionally not registered with the normal CTest suite because
it requires QMOMKit to be installed before the external test project is
configured. Run it separately when changing installation, export, package, or
public-header packaging behavior.

## Current Directory Structure

The current test layout is:

```text
tests/
    README.md
    unit/
        CMakeLists.txt
        conditional_moments.cpp
        configuration.cpp
        distributions.cpp
        moment_source.cpp
        node_advection.cpp
        output.cpp
        public_api.cpp
        moment_reconstruction.cpp
        realizability.cpp
        realizability_policy.cpp
        recursive_quadrature.cpp
        runge_kutta.cpp
        phase_space_evolution.cpp
        vandermonde.cpp
        wheeler.cpp
    component_regression/
        CMakeLists.txt
        README.md
        cqmom_compute_failures.cpp
        moment_source.cpp
        node_advection.cpp
        phase_space_evolution.cpp
    public_headers/
        CMakeLists.txt
        main.cpp
        include_*.cpp
    package_consumer/
        CMakeLists.txt
        README.md
        main.cpp
    verification/
        CMakeLists.txt
        README.md
        verification_helpers.hpp
        configuration/
            README.md
            configuration_failures.cpp
            *.qmom
        cqmom/
            README.md
            conditional_moments.cpp
            recursive_cqmom.cpp
            recursive_quadrature.cpp
            moment_reconstruction.cpp
            cqmom_2d_analytical.cpp
            cqmom_high_dimensional_analytical.cpp
            cqmom_public_workflows.cpp
            *.qmom
        distributions/
            README.md
            distributions_analytical.cpp
        moment_source/
            README.md
            moment_source_analytical.cpp
            *.qmom
        node_advection/
            README.md
            node_advection_analytical.cpp
            *.qmom
        output/
            README.md
            output_workflows.cpp
        qmom/
            README.md
            expected.hpp
            qmom_analytical.cpp
            qmom_compute_failures.cpp
            wheeler_analytical.cpp
            *.qmom
        realizability/
            README.md
            realizability_analytical.cpp
        runge_kutta/
            README.md
            runge_kutta_analytical.cpp
            runge_kutta_convergence.cpp
        phase-space evolution/
            README.md
            phase_space_evolution_analytical.cpp
            *.qmom
        vandermonde/
            README.md
            solver_selection.cpp
            *.qmom
```

## Building And Running Tests

The normal test suite is enabled from the repository root and links the current
build-tree QMOMKit target. Installation is not required.

```bash
cmake -S . -B build-tests \
    -DQMOMKIT_BUILD_TESTS=ON \
    -DCMAKE_BUILD_TYPE=Debug
cmake --build build-tests
```

Run the complete correctness suite:

```bash
ctest --test-dir build-tests --output-on-failure
```

Run one CTest category:

```bash
ctest --test-dir build-tests -L unit --output-on-failure
ctest --test-dir build-tests -L component_regression --output-on-failure
ctest --test-dir build-tests -L verification --output-on-failure
ctest --test-dir build-tests -L public_headers --output-on-failure
```

The labels are mutually exclusive in the current CMake configuration.

Run the installed-package smoke test:

```bash
cmake -S . -B build \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --target install

cmake -S tests/package_consumer \
    -B build-package-consumer \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build build-package-consumer
./build-package-consumer/consumer
```

This project behaves like an external user: it has no source-tree include path,
does not call `add_subdirectory`, and links only the installed
`qmomkit::qmomkit` target discovered by `find_package(QMOMKit CONFIG REQUIRED)`.

## Test Discovery And Selection

List registered tests without running them:

```bash
ctest --test-dir build-tests -N
```

Run tests whose registered CTest names match a regular expression:

```bash
ctest --test-dir build-tests -R <pattern> --output-on-failure
```

For example:

```bash
ctest --test-dir build-tests -R verification_qmom_analytical --output-on-failure
ctest --test-dir build-tests -R component_regression_node_advection --output-on-failure
ctest --test-dir build-tests -R unit_wheeler --output-on-failure
```

## Benchmarks

Benchmarks are timing tools, not correctness tests, and are not part of the
normal CTest suite. More details about benchmark build commands, workload descriptions, 
and interpretation guidance can be found in the [benchmark guide](../docs/developer/benchmarks.md).
