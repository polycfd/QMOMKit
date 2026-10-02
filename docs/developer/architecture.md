# Architecture

This page describes the current QMOMKit software architecture.

## Contents

- [Project Overview](#project-overview)
- [Package Architecture](#package-architecture)
- [Repository Structure](#repository-structure)
- [Public And Internal Boundaries](#public-and-internal-boundaries)
- [Consumer Model](#consumer-model)
- [Build Philosophy](#build-philosophy)
- [Design Principles](#design-principles)

## Project overview

QMOMKit is a C++ library for quadrature-based moment methods.
It turns raw moments of a number density function into quadrature
representations that can be embedded in larger simulation, validation, or
analysis codes.

The core library provides:

- one-dimensional QMOM reconstruction using the Wheeler algorithm;
- multidimensional CQMOM reconstruction with explicit phase-space hierarchy
  and conditioning order;
- analytical moment initialization helpers;
- realizability diagnostics and optional runtime acceptance predicates;
- SSP-RK3 and adaptive embedded RK solvers;
- phase-space advection of quadrature nodes, moment-source updates,
  and operator-split phase-space evolution composition;
- moment reconstruction/evaluation utilities, output helpers, and diagnostics.

The moment reconstruction utilities reconstruct native integer-order raw
moments from quadrature for internal and output workflows. The same public
module also evaluates generalized moments with finite real exponents as derived
quantities from existing QMOM or CQMOM quadrature; generalized moments are not
reconstruction inputs or phase-space-evolved states.

Moment generation, CQMOM reconstruction, and reconstruction-time storage are
separate responsibilities. Analytical distributions can produce `MomentSet`
data for QMOM, the required CQMOM moments for a hierarchy, or a full raw moment
tensor for callers that want one. CQMOM records the hierarchy, conditioning order,
and required-moment layout. Internal recursive storage adapts those moments for
efficient traversal without making distribution objects aware of CQMOM.

QMOMKit is not a CFD solver. It does not own meshes, cells, finite-volume or
finite-element operators, boundary conditions, MPI communication,
application-specific source terms, or particle-specific physics.

QMOMKit is the project and package identity. QMOM and CQMOM are numerical
methods implemented by the library. The preferred public reconstruction APIs
are:

```text
compute_qmom(...)   one-dimensional QMOM
compute_cqmom(...)  multidimensional CQMOM
```

The recommended configured QMOM reconstruction workflow is:

```text
configuration
    -> load_qmom(...)
    -> QmomCase
    -> compute_qmom(QmomCase)
    -> QmomResult
```

The recommended configured CQMOM workflow is:

```text
configuration
    -> load_cqmom_prepared(...)
    -> CqmomPreparedCase
    -> compute_cqmom(CqmomPreparedCase)
    -> CqmomResult / CqmomQuadrature
```

After reconstruction, both QMOM and CQMOM workflows may inspect quadrature,
evaluate derived quantities, write output, or enter optional phase-space
evolution. QMOM phase-space evolution creates a `QmomState` with
`QmomCase::initial_state()` and advances it with
`advance_qmom_phase_space(QmomState, ...)`; normal CQMOM phase-space evolution
creates a `CqmomState` with `CqmomPreparedCase::initial_state()` and advances
it with `advance_cqmom_phase_space(CqmomState, ...)`.

Accepted-state batch phase-space evolution is a serial collection of these same
local operations. It does not give QMOMKit ownership of cells, meshes,
partitions, or parallel execution; applications retain spatial organization and
scheduling.

`CqmomPreparedCase` and `CqmomState` are the public objects for normal CQMOM
workflows. Normal CQMOM phase-space evolution creates a `CqmomState` from
`CqmomPreparedCase` and advances it with `advance_cqmom_phase_space(...)`. Both
case types hold the moments the hierarchy reconstructs from; a prepared case
also holds the recursive branch structure that hierarchy implies, so repeated
reconstructions reuse it instead of rebuilding it. `CqmomCase` is the CQMOM
counterpart of `QmomCase`: an open case definition, with its moments in
`CqmomMomentLayout` and `CqmomMomentData`, that a caller builds or edits in code
and then prepares. It can also be reconstructed or turned into a state
directly, but each such call rebuilds the branch structure, so the prepared case
is the normal route.

### QMOM and CQMOM terminology

Use these QMOM and CQMOM terms consistently in developer documentation and
implementation comments:

| Term | Meaning |
|------|---------|
| Case | Problem definition used by a solver workflow, including initial moments and reconstruction settings. `QmomCase` is the configured one-dimensional QMOM case. `CqmomPreparedCase` is the normal CQMOM case, holding the prepared recursive branch structure alongside the moments. `CqmomCase` is the open CQMOM case definition, mirroring `QmomCase`, from which a prepared case is built. |
| State | Evolving runtime state during time integration or phase-space evolution. QMOM uses `QmomState`; CQMOM uses `CqmomState`. |
| Raw moment tensor | Full rectangular moment representation, `RawMomentTensor`. It is an input and an output format — distributions and sample populations produce one, and `reconstruct_raw_moments(...)` returns one — not a representation a case or state holds. |
| Required moments | The moments a hierarchy reconstructs from, held by `CqmomMomentData` with its `CqmomMomentLayout`. This is what a case and a state store. |
| MomentData | Storage container for CQMOM required moments associated with a hierarchy, represented publicly by `CqmomMomentData`. |
| Request | Mathematical query/request for moments, such as `RawMomentRequest`. A request asks for moment values; it is not a QMOM or CQMOM problem definition. |
| Quadrature | Discrete representation reconstructed from moments, such as weights and abscissas in `QmomResult`, `Quadrature1D`, or `CqmomQuadrature`. |

`QmomCase` and `CqmomPreparedCase` are problem definitions, not moment
requests. `CqmomCase` is the case definition a caller builds or edits
directly; it does not replace `CqmomPreparedCase` and `CqmomState` as the normal
CQMOM workflow abstractions.

A caller holding a full raw moment tensor selects the required moments from it
with `extract_cqmom_moments(...)`, or hands the tensor to
`config::with_initial_moments(...)`, which does the same selection and reports a
malformed tensor as a configuration error.

For mathematical storage and ordering rules, see
`numerical_conventions.md`.

## Package architecture

### QMOMKit core package

The repository root CMake project builds the core library target:

```text
qmomkit
```

and exports the package target:

```text
qmomkit::qmomkit
```

The core package installs:

- the compiled library target;
- public headers from `include/qmomkit/`;
- `QMOMKitConfig.cmake`;
- `QMOMKitConfigVersion.cmake`;
- exported `QMOMKitTargets.cmake`.

Downstream CMake projects consume the installed package with:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE qmomkit::qmomkit)
```

The default root CMake build builds and installs only the core package. Normal
source-tree tests are added only when `QMOMKIT_BUILD_TESTS=ON`. Benchmarks,
examples, installed-package consumers, and Monte Carlo tooling remain separate
projects.

### QMOMKitMonteCarlo companion package

`tools/monte_carlo/` contains the optional QMOMKitMonteCarlo companion package.
It provides Monte Carlo sampling, ensemble moment extraction, evolution,
output, and comparison utilities for validation workflows.

The companion package exports:

```text
qmomkit_monte_carlo::qmomkit_monte_carlo
```

It is an installed-package consumer of QMOMKit and links the core package as a
public dependency:

```text
QMOMKitMonteCarlo
        |
      PUBLIC
        |
        v
qmomkit::qmomkit
```

The companion package configuration calls `find_dependency(QMOMKit CONFIG
REQUIRED)`, so downstream projects that consume QMOMKitMonteCarlo receive the
core dependency through the exported Monte Carlo target.

QMOMKitMonteCarlo is not part of the core library. Its headers live
under `tools/monte_carlo/include/monte_carlo/`, not under
`include/qmomkit/`.

## Repository structure

The repository is organized around a small core package and standalone
consumers:

```text
include/
src/
tests/
tests/benchmarks/
tools/monte_carlo/
examples/
docs/
```

- `include/`: public QMOMKit headers. The installed core API lives under
  `include/qmomkit/`.
- `src/`: core library implementation files and implementation-only helper
  headers.
- `tests/`: CTest unit, component regression, public-header, and verification
  executables included by the root build when `QMOMKIT_BUILD_TESTS=ON`.
- `tests/benchmarks/`: standalone optional benchmark consumer.
- `tools/monte_carlo/`: standalone QMOMKitMonteCarlo companion package and
  its tests, examples, scripts, headers, and sources.
- `examples/`: standalone example applications. Core examples consume
  QMOMKit, Monte Carlo examples consume QMOMKitMonteCarlo.
- `docs/`: user and developer documentation.

Generated build and install directories may exist inside these trees during
local development. They are not architectural components.

## Public and internal boundaries

### Public API

The supported public API is located in:

```text
include/qmomkit/
```

Installed public headers represent functionality that downstream users may
include with:

```cpp
#include <qmomkit/...>
```

The umbrella header is:

```cpp
#include <qmomkit/qmomkit.hpp>
```

Prefer narrower headers when a caller only needs a specific API family, such
as `<qmomkit/qmom.hpp>`, `<qmomkit/cqmom.hpp>`,
`<qmomkit/distributions.hpp>`, or `<qmomkit/phase_space_evolution.hpp>`.

Public headers should expose stable types, requests, results, options, and
functions. They should not expose implementation-only utilities or private
algorithm plumbing.

### Internal implementation

Core implementation files live in:

```text
src/
src/internal/
```

`src/` contains algorithm implementations and module ownership close to the
public APIs they support. `src/internal/` contains shared implementation-only
utilities and private lower-level building blocks.

Headers under `src/internal/` are not installed and are not part of the
supported downstream API. Public headers under `include/qmomkit/` should not
include `src/internal/` headers.

Current internal utilities include small dense linear-algebra and numeric
helpers used by multiple implementation modules. Internal utilities should
remain narrowly scoped and should not become a catch-all location for unrelated
algorithm code.

Some repository tests include internal headers to validate lower-level
behavior. That does not make those headers public API.

## Consumer model

The root package is consumed through its installed CMake package:

```text
QMOMKit install
    |
    +-- tests/package_consumer
    +-- tests/benchmarks
    +-- examples

QMOMKitMonteCarlo
    |
    +-- MC examples
```

Normal correctness tests are source-tree tests. They are enabled from the root
build with `QMOMKIT_BUILD_TESTS=ON` and link the local `qmomkit` target.

The standalone consumers use `find_package(...)` rather than building the
source tree as a dependency:

- `tests/package_consumer/` uses `find_package(QMOMKit CONFIG REQUIRED)`.
- `tests/benchmarks/` uses `find_package(QMOMKit CONFIG REQUIRED)`.
- core examples use `find_package(QMOMKit CONFIG REQUIRED)` or the already
  available `qmomkit::qmomkit` target.
- Monte Carlo examples use `find_package(QMOMKitMonteCarlo CONFIG REQUIRED)`.
- `tools/monte_carlo/` uses `find_package(QMOMKit CONFIG REQUIRED)` and
  exports its own package.

This model keeps package-consumer behavior visible while ensuring ordinary
correctness tests exercise the current source-tree implementation. Repository
tests may add `src/` to their include path for focused internal validation, but
downstream consumers should rely on installed headers and exported targets.

## Build philosophy

The root CMake project is the core-library/package build. It compiles the
QMOMKit library and installs the exported CMake package.

Tests, benchmarks, tools, and examples are built according to their purpose:

- correctness tests are ordinary CTest executables included by the root build
  when `QMOMKIT_BUILD_TESTS=ON`;
- benchmarks are optional performance executables in `tests/benchmarks/`;
- the Monte Carlo companion package is built and installed from
  `tools/monte_carlo/`;
- examples are built from their own directories.

The core build has no mandatory third-party numerical dependency. Benchmark
support may use Google Benchmark in the standalone benchmark project only.

### Canonical build modes

Use an explicit build type for reproducible development with
single-configuration generators such as Unix Makefiles or Ninja. When the root
QMOMKit library is configured with a single-configuration generator and
`CMAKE_BUILD_TYPE` is empty, the project defaults to `Release`. Explicit build
types should still be used for routine development, testing, profiling, and
benchmarking so the intended configuration is visible in the command history.

Debug builds are for development and debugging. They should keep assertions
enabled, leave `NDEBUG` undefined, include debug symbols, and avoid aggressive
optimization:

```sh
cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug
```

Release builds are for production-style usage and performance benchmarks. They
should enable compiler optimization, define `NDEBUG`, and disable standard C++
assertions:

```sh
cmake -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

RelWithDebInfo builds are for profiling and performance analysis. They should
keep optimization enabled, retain debug symbols, define `NDEBUG`, and preserve
a profile-friendly symbolized binary:

```sh
cmake -B build-profile -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build-profile
```

Tests, benchmarks, examples, and companion tools consume an installed QMOMKit
package. Configure those standalone consumer builds with an explicit build type
too, and use an installed QMOMKit library built with the configuration intended
for the task.

## Design principles

QMOMKit architecture follows these stable principles:

- Keep the core library CFD-code agnostic and application agnostic.
- Keep public mathematical conventions explicit and documented.
- Preserve a clear boundary between installed public headers and internal
  implementation files.
- Prefer simple runtime-flexible abstractions over hard-coded dimensions,
  hard-coded node counts, or excessive template metaprogramming.
- The core code of QMOMKit should be free of third-party dependencies.
- Use explicit result/status objects and diagnostics for numerical failures.
- Validate numerical behavior through deterministic analytical,
  manufactured, and independent Monte Carlo workflows where appropriate.
- Keep examples, tests, benchmarks, and validation tools as package consumers
  so installed-package behavior remains exercised.
- Keep user-facing APIs simple while allowing lower-level implementation
  diagnostics to support maintainers and focused tests.
