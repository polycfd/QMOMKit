# Testing

This page documents the current QMOMKit testing philosophy and test
architecture for maintainers, contributors, and AI-assisted development
workflows. It describes stable correctness expectations, not a complete test
inventory or temporary CI setup.

For package layout and consumer boundaries, see `architecture.md`. For moment,
tensor, and phase-space conventions that tests must protect, see
`numerical_conventions.md`. For realizability status behavior and runtime
acceptance policies, see `realizability.md`. For Monte Carlo companion package
validation guidance, see `monte_carlo.md`.

## Testing Philosophy

QMOMKit is a numerical C++ library, so correctness requires more than exercising
code paths. A passing test should usually demonstrate a mathematical invariant,
a known analytical result, a manufactured solution, or a documented failure
contract.

Unit tests alone are insufficient for moment-method software. They are good at
protecting local behavior, option validation, indexing helpers, and failure
paths, but they do not by themselves establish that an inversion, conditional
reconstruction, transport update, or realizability check matches the intended
mathematics.

The core testing principle is:

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

Analytical and manufactured solutions are especially important because they make
the expected answer independent of the implementation under test. When a
verification compares recovered quadrature, evolved moments, or reconstructed
raw moments, the reference values should come from the mathematical case being
tested, not from another call into the same production algorithm.

When adding a new test, choose the directory by purpose:

1. If it verifies an isolated function, class, or local API contract, use
   `tests/unit/`.
2. If it protects specific implementation behavior or a regression-prone
   component interaction, use `tests/component_regression/`.
3. If it verifies a mathematical property or user workflow against an
   analytical, manufactured, convergence, or otherwise independent reference,
   use `tests/verification/`.

## Test Architecture

The primary test project lives under:

```text
tests/
```

It is included by the root CMake project when `QMOMKIT_BUILD_TESTS=ON`. Normal
unit, component regression, verification, and public-header tests link the
current build-tree `qmomkit` target and compile against the current source-tree
public headers. Installation is not required for the normal correctness suite:

```bash
cmake -S . -B build-tests \
    -DQMOMKIT_BUILD_TESTS=ON \
    -DCMAKE_BUILD_TYPE=Debug
cmake --build build-tests
ctest --test-dir build-tests --output-on-failure
```

The test project currently registers CTest executables through these
subdirectories:

```text
tests/
    unit/
    component_regression/
    verification/
    public_headers/
```

Some repository tests include `src/` as a private include path for focused
internal validation. That is permitted inside the repository test project, but
it does not change the public API boundary. Downstream users should consume only
installed headers and exported package targets.

Package consumption is also tested explicitly by:

```text
tests/package_consumer/
```

That directory is a minimal external CMake project. It validates installed
package discovery, exported target linkage, and installed public-header usage
without repository source include paths or `add_subdirectory`.

## Test Categories

### Unit Tests

Unit tests verify isolated functions, classes, data structures, options, and
local API contracts. They should be focused, fast, deterministic, and specific
about the contract being protected.

Unit coverage is appropriate for:

- request and option validation;
- status values, diagnostics, and `ok()` behavior;
- indexing and storage helpers;
- local algorithm invariants;
- error handling and state preservation;
- boundary cases that should fail cleanly.

Representative unit areas include Wheeler inversion, recursive CQMOM support
structures, realizability checks and policies, configuration mapping,
distribution helpers, output helpers, node advection, moment-source evolution,
raw-moment reconstruction, and transport components.

### Component Regression Tests

Component regression tests protect important implementation behavior that is
larger than an isolated unit contract but narrower than mathematical
verification. They are appropriate for direct component composition,
diagnostics, status propagation, policy mapping, and previously fragile
implementation paths that should remain visible during maintenance.

Component regression coverage belongs under:

```text
tests/component_regression/
```

These tests may use manufactured inputs or simple analytical expectations, but
their primary purpose is regression protection. They are not failed
verification tests, not a substitute for unit tests, and not replacements for
verification tests when the behavior is a public mathematical workflow with an
independent reference.

### Verification Tests

Verification tests provide mathematical evidence that the numerical algorithms
behave as intended. They compare against analytical solutions, manufactured
solutions, reconstruction identities, convergence behavior, or other
independently known results.

Verification tests serve as executable documentation of mathematical
correctness and public workflows.

Verification coverage belongs around:

- QMOM and Wheeler reconstruction;
- CQMOM conditional moments and recursive quadrature;
- distribution-generated raw moments;
- support-aware realizability checks;
- Runge-Kutta integration;
- node advection;
- moment-source evolution;
- transport composition.

Verification tests should preserve QMOMKit conventions: raw moments are moments
of the number density function, tensor indices follow phase-space dimension
order, conditioning order does not redefine phase-space order, and quadrature
weights are not normalized unless the case explicitly defines unit total number
density.

### Public-Header Tests

Public-header tests are compile-time interface checks. Each translation unit in
`tests/public_headers/` includes one current public header from
`include/qmomkit/` through the source-tree `qmomkit` target.

These tests protect the public API boundary by verifying that headers are
self-contained and do not depend on accidental include order. The linked CTest
executable does not need numerical runtime assertions; successful compilation
and linkage are the meaningful checks.

When a public header is added, renamed, or given new direct dependencies, the
public-header coverage should be updated so downstream include behavior remains
visible.

### Package Consumer Tests

Package consumer tests validate the installed-package experience rather than the
in-tree build. The consumer project should behave like an external application:

- discover QMOMKit with `find_package(QMOMKit CONFIG REQUIRED)`;
- link only `qmomkit::qmomkit`;
- include installed public headers;
- avoid repository source paths, internal headers, and test utilities.

`tests/package_consumer/` is intentionally separate from the normal CTest unit,
verification, and public-header suite. It represents an external consumer
workflow and requires QMOMKit to be installed before the consumer project is
configured. Its purpose is to validate installed package discovery and exported
CMake target behavior, not correctness or regression behavior inside the normal
`tests/` project.

This category is most important when package exports, installation rules,
public headers, or CMake package metadata change.

## Numerical Validation Principles

Analytical solutions are preferred when the exact answer is known in closed
form. They are useful for Wheeler inversion, low-dimensional CQMOM cases,
distribution moment generation, simple ODE integration, source evolution, and
transport cases with closed-form operator behavior.

Manufactured solutions are useful when the algorithm is more complex than the
available closed-form examples. They should define moments, nodes, weights,
conditioning paths, source terms, or transport behavior from an independent
construction and then verify that QMOMKit recovers the expected mathematical
object.

Convergence studies are useful for time integration and transport. They should
test numerical order or monotonic error reduction using deterministic cases
whose exact or high-confidence reference behavior is known.

Deterministic regression cases are useful for boundary behavior, diagnostics,
and previously fragile numerical situations. They should preserve the meaningful
status, tolerance, ordering, or reconstruction property rather than simply
locking in implementation-shaped output.

High-dimensional CQMOM needs special care. Tests should separate phase-space
dimension order from conditioning order, verify terminal quadrature weights and
coordinates, and reconstruct raw tensor moments from the final quadrature when
that is the mathematical guarantee being tested.

Realizability tests should cover accepted, boundary, rejected, invalid-input,
and numerically indeterminate cases. They should also verify that policies reject
candidate moment states without implying automatic correction.

## Adding New Tests

When adding functionality, add tests at the level where the new behavior is
owned.

Add unit tests for local behavior: request validation, option defaults, status
values, diagnostics, indexing, state preservation, and explicit failure paths.

Add verification tests for mathematical guarantees: moment reconstruction,
closed-form evolution, manufactured quadrature recovery, support-aware
realizability, convergence, or invariant preservation.

Add public-header coverage when public APIs change. Installed headers should
compile when included directly and should not rely on a caller including another
header first.

Add or update package-consumer coverage when installation, export, package
configuration, or public target behavior changes.

The optional Monte Carlo companion package has its own CMake project and tests
under `tools/monte_carlo/`. Configure that project with `BUILD_TESTING=ON` when
running its tests, then run them explicitly with CTest.

Update examples only when the user-facing workflow changes. Examples are not a
substitute for automated correctness tests.

Prefer deterministic tests. Stochastic or Monte Carlo workflows are valuable for
independent validation, but normal correctness tests should not depend on random
sampling noise unless the uncertainty model is explicit and controlled.

## Relationship With Benchmarks

Tests establish correctness. Benchmarks measure performance.

A benchmark should not be the only evidence that an algorithm is correct, and a
correctness test should not depend on timing-sensitive behavior. Performance
work should preserve the same mathematical contracts checked by unit and
verification tests.

Benchmark structure and interpretation belong in `benchmarks.md`.
