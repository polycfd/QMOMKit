# Testing

This page documents the current QMOMKit testing philosophy and test
architecture.
It describes stable correctness expectations and the current root
CTest inventory, not temporary CI setup.

The bar set here is deliberately high. Today, code can be produced faster than
evidence that it is correct, so QMOMKit treats the test suite, not the
implementation, as what establishes that its numerical methods behave as
intended. The reasoning behind that position is recorded in
[DEVELOPMENT_PHILOSOPHY.md](../../DEVELOPMENT_PHILOSOPHY.md).

## Contents

- [Debug Build Requirement](#debug-build-requirement)
- [Test-suite Organization](#test-suite-organization)
- [Test Categories](#test-categories)
- [Numerical Validation Principles](#numerical-validation-principles)
- [Adding New Tests](#adding-new-tests)
- [Relationship With Benchmarks](#relationship-with-benchmarks)
- [Registered Root CTest Inventory](#registered-root-ctest-inventory)
- [Test-suite Maintenance](#test-suite-maintenance)

## Debug build requirement

The unit tests use C/C++ `assert()`. `assert()` is disabled when `NDEBUG` is
defined, and normal Release builds define `NDEBUG`; therefore a Release CTest
run does not validate unit assertions. The authoritative functional test
configuration is Debug:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DQMOMKIT_BUILD_TESTS=ON

cmake --build build --parallel

ctest --test-dir build --output-on-failure
```

Run one registered unit contract with:

```bash
ctest --test-dir build -R '^unit_wheeler_nonfinite_moment_reports_invalid_input$' --output-on-failure
```

Keep the existing assertion-based unit tests intact. Do not use a Release CTest
result as evidence that the functional unit suite passed.

## Test-suite organization

QMOMKit is a numerical C++ library, so correctness requires more than
exercising code paths. A passing test should usually demonstrate a mathematical
invariant, a known analytical result, a manufactured solution, or a documented
failure contract.

Unit tests alone are insufficient for moment-method software. They are good at
protecting local behavior, option validation, indexing helpers, and failure
paths, but they do not by themselves establish that an inversion, conditional
reconstruction, phase-space evolution update, or realizability check matches
the intended mathematics.

The core testing principle is:

```text
QMOMKit testing

                        Test purpose
                              |
        ---------------------------------------------
        |                     |                     |
      unit             component_regression     verification
        |                     |                     |
 isolated correctness    implementation      mathematical and
 local APIs                stability        workflow correctness
```

Analytical and manufactured solutions are especially important because they
make the expected answer independent of the implementation under test. When a
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

### Test project structure

The primary test project lives under:

```text
tests/
```

It is included by the root CMake project when `QMOMKIT_BUILD_TESTS=ON`. Normal
unit, component regression, verification, and public-header tests link the
current build-tree `qmomkit` target and compile against the current source-tree
public headers. Installation is not required for the normal correctness suite.

### Running tests

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

### CTest naming and current count

CTest names identify the tested contract. Several executables accept a local
selector so independently meaningful groups can be registered separately. The
individual CTest registrations are authoritative; running such an executable
without a selector is a developer convenience, not an aggregate CTest contract.
The registration counts are recorded once, under
[Registered root CTest inventory](#registered-root-ctest-inventory), and the
complete logical-test inventory is maintained in
[Test inventory](test_inventory.md); this page remains the conceptual and
how-to guide.
Obtain the current count from `ctest -N` (and confirm structured registration
with `ctest --show-only=json-v1`), rather than maintaining it by arithmetic
alone.

CTest registrations are defined from these subdirectories:

```text
tests/
    unit/
    component_regression/
    verification/
    public_headers/
```

Some repository tests include `src/` as a private include path for focused
internal validation. That is permitted inside the repository test project, but
it does not change the public API boundary. Downstream users should consume
only installed headers and exported package targets.

Package consumption is also tested explicitly by:

```text
tests/package_consumer/
```

That directory is a minimal external CMake project. It validates installed
package discovery, exported target linkage, and installed public-header usage
without repository source include paths or `add_subdirectory`.

## Test categories

### Unit tests

Unit tests verify isolated functions, classes, data structures, options, and
local API contracts. They should be focused, fast, deterministic, and specific
about the contract being protected.

They intentionally use `assert()` and therefore require a Debug build.

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
raw-moment reconstruction, and phase-space evolution components.

### Component regression tests

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

Independent numerical verification requires an analytical, manufactured,
convergence, conservation, or otherwise independent reference. Batch/scalar
parity, round trips, captured fixtures, and comparisons between related
production paths remain valuable consistency/regression evidence, but must not
be presented as independent algorithm verification.

### Verification tests

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
- phase-space evolution composition.

Verification tests should preserve QMOMKit conventions: raw moments are moments
of the number density function, tensor indices follow phase-space dimension
order, conditioning order does not redefine phase-space order, and quadrature
weights are not normalized unless the case explicitly defines unit total number
density.

### Public-header tests

Public-header tests are compile-time interface checks. Each translation unit in
`tests/public_headers/` includes one current public header from
`include/qmomkit/` through the source-tree `qmomkit` target. That header comes
first, before any standard header the translation unit also needs, so a standard
include missing from the header cannot be supplied by the test. The
repository's clang-format configuration regroups includes alphabetically, so
where a standard header would otherwise sort ahead of the public header, such as
`<cstddef>` or `<functional>`, a comment line separates them.

These tests protect the public API boundary by verifying that headers are
self-contained and do not depend on accidental include order. The linked CTest
executable does not need numerical runtime assertions; successful compilation
and linkage are the meaningful checks.

Some translation units also pin entry points. Each pin stores the address of a
function in a pointer of its exact signature, so a changed parameter or return
type fails to compile. Pins must have external linkage, for example in a named
namespace: the linker then has to find each function's definition, so a
function that is declared but never defined fails the build. In an unnamed
namespace the compiler may drop unused pins, as Clang does even without
optimization and GCC does with it, so such pins check only the signature.

When a public header is added, renamed, or given new direct dependencies, the
public-header coverage should be updated so downstream include behavior remains
visible.

### Package consumer tests

Package consumer tests validate the installed-package experience rather than
the in-tree build. The consumer project should behave like an external
application:

- discover QMOMKit with `find_package(QMOMKit CONFIG REQUIRED)`;
- link only `qmomkit::qmomkit`;
- include installed public headers;
- avoid repository source paths, internal headers, and test utilities.

`tests/package_consumer/` is intentionally separate from the normal CTest unit,
verification, and public-header suite. It represents an external consumer
workflow and requires QMOMKit to be installed before the consumer project is
configured. Its purpose is to validate installed package discovery and exported
CMake target behavior, not correctness or regression behavior inside the normal
`tests/` project. The `CI-Consumers` workflow installs QMOMKit in Debug and
builds and runs the consumer on every push to `main` and on pull requests that
change more than `docs/`, `examples/`, or `paper/`.

This category is most important when package exports, installation rules,
public headers, or CMake package metadata change.

## Numerical validation principles

Analytical solutions are preferred when the exact answer is known in closed
form. They are useful for Wheeler inversion, low-dimensional CQMOM cases,
distribution moment generation, simple ODE integration, moment-source
evolution, and phase-space evolution cases with closed-form operator behavior.

Manufactured solutions are useful when the algorithm is more complex than the
available closed-form examples. They should define moments, nodes, weights,
conditioning paths, source terms, or phase-space evolution behavior from an
independent construction and then verify that QMOMKit recovers the expected
mathematical object.

Convergence studies are useful for time integration and phase-space evolution.
They should test numerical order or monotonic error reduction using
deterministic cases whose exact or high-confidence reference behavior is known.

A convergence registration must identify its independent reference, refinement
variable and sequence, error quantity, expected order or trend, and pass
threshold. Tolerance refinement for adaptive methods is accuracy evidence; it
is not a formal fixed-step order proof.

Deterministic regression cases are useful for boundary behavior, diagnostics,
and previously fragile numerical situations. They should preserve the
meaningful status, tolerance, ordering, or reconstruction property rather than
simply locking in implementation-shaped output.

High-dimensional CQMOM needs special care. Tests should separate phase-space
dimension order from conditioning order, verify terminal quadrature weights and
coordinates, and reconstruct raw tensor moments from the final quadrature when
that is the mathematical guarantee being tested.

Realizability tests should cover accepted, boundary, rejected, invalid-input,
and numerically indeterminate cases. They should also verify that policies
reject candidate moment states without implying automatic correction.

## Adding new tests

When adding functionality, add tests at the level where the new behavior is
owned.

Add unit tests for local behavior: request validation, option defaults, status
values, diagnostics, indexing, state preservation, and explicit failure paths.

Add verification tests for mathematical guarantees: moment reconstruction,
closed-form evolution, manufactured quadrature recovery, support-aware
realizability, convergence, or invariant preservation.

Add public-header coverage when public APIs change. Installed headers should
compile when included directly and should not rely on a caller including
another header first.

Add or update package-consumer coverage when installation, export, package
configuration, or public target behavior changes.

The optional Monte Carlo companion package has its own CMake project and tests
under `tools/monte_carlo/`. Configure that project with `BUILD_TESTING=ON` when
running its tests, then run them explicitly with CTest. It follows the same
registration convention as the root suite: each logical case is registered with
CTest on its own as `monte_carlo_<subsystem>_<case>`, discovered from the
`void test_*()` definitions in each test source, so a failure names the behavior
rather than the file. Those registrations are not part of the root inventory,
which covers root CTest only. The `CI-Consumers` workflow builds these tests in
Debug against the installed QMOMKit and runs them alongside the package
consumer.

Update examples only when the user-facing workflow changes. Examples are not a
substitute for automated correctness tests.

Prefer deterministic tests. Stochastic or Monte Carlo workflows are valuable
for independent validation, but normal correctness tests should not depend on
random sampling noise unless the uncertainty model is explicit and controlled.

## Relationship with benchmarks

Tests establish correctness. Benchmarks measure performance.

A benchmark should not be the only evidence that an algorithm is correct, and a
correctness test should not depend on timing-sensitive behavior. Performance
work should preserve the same mathematical contracts checked by unit and
verification tests.

Benchmark structure and interpretation belong in
[benchmarks.md](benchmarks.md).

## Registered root CTest inventory

The root inventory contains 850 CTest registrations: 849 logical runtime tests
across the three runtime tiers (634 unit, 47 component-regression, and 168
verification), plus the `public_headers` supporting compile/link check. CTest
names use the tier prefix followed by a subsystem and a descriptive case name.

The CMake registration helpers in `tests/CMakeLists.txt` keep the target build
definition separate from the logical CTest cases. A selector passed by CTest
must execute exactly the named logical case; an executable's no-argument mode
is retained only for local developer convenience.

Use the following commands to inspect the authoritative inventory rather than
maintaining a hand-written expansion of every individual case:

```bash
ctest --test-dir build -N
ctest --test-dir build -N -L unit
ctest --test-dir build -N -L component_regression
ctest --test-dir build -N -L verification
ctest --test-dir build -N -L public_headers
ctest --test-dir build --show-only=json-v1
```

### External installed-package consumer check

`tests/package_consumer/` is outside the root CTest count. After installing
QMOMKit, configure and build it as an external project to verify
`find_package(QMOMKit CONFIG REQUIRED)`, `qmomkit::qmomkit`, and installed
headers without repository include paths.

## Test-suite maintenance

When adding, removing, renaming, or splitting a root CTest registration, update
this inventory in the same change and give every entry a distinct contract and
method. Re-run `ctest -N` and `ctest --show-only=json-v1` to record the actual
count. Keep functional validation in Debug because unit assertions are disabled
by `NDEBUG`; do not replace the assertion mechanism without an explicit design
request.
