# Architecture

This page describes the current QMOMKit software architecture for maintainers,
contributors, and AI-assisted development workflows. It documents stable
package boundaries and repository responsibilities, not development history.

## Project Overview

QMOMKit is a C++ library for quadrature-based moment methods.
It turns raw moments of a number density function into quadrature
representations that can be embedded in larger simulation, validation, or
analysis codes.

The core library provides:

- one-dimensional QMOM reconstruction through Wheeler quadrature;
- multidimensional CQMOM reconstruction with explicit phase-space hierarchy
  and conditioning order;
- analytical moment initialization helpers;
- realizability diagnostics and optional runtime acceptance predicates;
- generic SSP-RK3 and adaptive embedded RK23/RK45 solvers;
- continuous quadrature-node advection, discontinuous moment-source updates,
  and staged moment transport composition;
- moment reconstruction/evaluation utilities, output helpers, and diagnostics.

The moment reconstruction utilities reconstruct native integer-order raw
moments from quadrature for internal and output workflows. The same public
module also evaluates generalized moments with finite real exponents as derived
quantities from existing QMOM or CQMOM quadrature; generalized moments are not
reconstruction inputs or transported state.

Moment generation, CQMOM reconstruction, and reconstruction-time storage are
separate responsibilities. Analytical distributions can produce `MomentSet`
data for QMOM, selected required CQMOM moments, or dense full moment tensors for
explicit dense workflows. CQMOM records the hierarchy, conditioning order, and
required-moment layout. Internal recursive storage adapts those moments for
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
    -> load_cqmom(...)
    -> CqmomCase
    -> compute_cqmom(CqmomCase)
    -> CqmomResult / CqmomQuadrature
```

After reconstruction, both QMOM and CQMOM workflows may inspect quadrature,
evaluate derived quantities, write output, or enter optional transport. QMOM
transport creates a `QmomState` with `QmomCase::initial_state()` and advances it
with `advance_transport_qmom(QmomState, ...)`; normal CQMOM transport creates a
`CqmomState` with `CqmomCase::initial_state()` and advances it with
`advance_transport_cqmom(CqmomState, ...)`.

`CqmomCase` and `CqmomState` are the neutral public objects for normal CQMOM
workflows. Normal CQMOM transport creates a `CqmomState` from `CqmomCase` and
advances it with `advance_transport_cqmom(...)`, independent of dense or lean
storage. Normal configured workflows default to lean moment storage and may
select dense storage explicitly with `cqmom.moment_storage`. `CqmomDenseCase`
with `RawMomentTensor` and `CqmomLeanCase` with `CqmomMomentLayout` and
`CqmomMomentData` remain supported advanced explicit representations for static
inversion, specialized initialization, validation, and expert workflows.

### QMOM And CQMOM Terminology

Use these QMOM and CQMOM terms consistently in developer documentation and
implementation comments:

| Term | Meaning |
|------|---------|
| Case | Problem definition used by a solver workflow, including initial moments and reconstruction settings. `QmomCase` is the configured one-dimensional QMOM case. `CqmomCase` is the normal representation-independent CQMOM case. |
| State | Evolving runtime state during time integration or transport. QMOM uses `QmomState`; CQMOM uses `CqmomState`. |
| Dense | Full tensor moment representation, normally `RawMomentTensor`. Dense is not legacy; it is an explicit full-tensor workflow for CQMOM callers that choose that representation. |
| Lean | Compact CQMOM required-moment representation containing only the moments required for a given hierarchy. Lean is not a different physical model. |
| Explicit representation workflow | Workflow where the user intentionally chooses and manages the CQMOM moment representation, either dense or lean. |
| MomentData | Storage container for CQMOM required moments associated with a hierarchy, represented publicly by `CqmomMomentData`. |
| Request | Mathematical query/request for moments, such as `RawMomentRequest` or `ConditionalMomentRequest`. A request asks for moment values or one conditional-moment reconstruction; it is not a QMOM or CQMOM problem definition. |
| Quadrature | Discrete representation reconstructed from moments, such as weights and abscissas in `QmomResult`, `Quadrature1D`, or `CqmomQuadrature`. |

`QmomCase` and `CqmomCase` are problem definitions, not moment requests.
`CqmomDenseCase` and `CqmomLeanCase` are advanced explicit representation
workflows for expert callers that need to manage dense full tensors or lean
required moment data directly. They do not replace `CqmomCase` and
`CqmomState` as the normal representation-independent CQMOM workflow
abstractions.

Dense and lean terminology is specific to CQMOM's multidimensional moment
storage choices. QMOM does not need dense or lean case terminology because its
one-dimensional setup is represented directly by `MomentSet`, `node_count`,
and `QmomOptions` inside `QmomCase`.

For mathematical storage and ordering rules, see
`numerical_conventions.md`.

## Package Architecture

### QMOMKit Core Package

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

### QMOMKitMonteCarlo Companion Package

`tools/monte_carlo/` contains the optional QMOMKitMonteCarlo companion package.
It provides Monte Carlo sampling, ensemble moment extraction, evolution, output,
and comparison utilities for validation workflows.

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

QMOMKitMonteCarlo is not part of the core library surface. Its headers live
under `tools/monte_carlo/include/monte_carlo/`, not under
`include/qmomkit/`.

## Repository Structure

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
  QMOMKit; Monte Carlo examples consume QMOMKitMonteCarlo.
- `docs/`: user and developer documentation.

Generated build and install directories may exist inside these trees during
local development. They are not architectural components.

## Public And Internal Boundaries

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
`<qmomkit/distributions.hpp>`, or `<qmomkit/transport_integrator.hpp>`.

Public headers should expose stable types, requests, results, options, and
functions. They should not expose implementation-only utilities or private
algorithm plumbing.

### Internal Implementation

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

## Consumer Model

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

## Build Philosophy

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

### Canonical Build Modes

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
keep optimization enabled, retain debug symbols, define `NDEBUG`, and preserve a
profile-friendly symbolized binary:

```sh
cmake -B build-profile -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build-profile
```

Tests, benchmarks, examples, and companion tools consume an installed QMOMKit
package. Configure those standalone consumer builds with an explicit build type
too, and use an installed QMOMKit library built with the configuration intended
for the task.

## Design Principles

QMOMKit architecture follows these stable principles:

- Keep the core library CFD-code agnostic and application agnostic.
- Keep public mathematical conventions explicit and documented.
- Preserve a clear boundary between installed public headers and internal
  implementation files.
- Prefer simple runtime-flexible abstractions over hard-coded dimensions,
  hard-coded node counts, or excessive template metaprogramming.
- Keep dependencies minimal; do not require Eigen, BLAS, LAPACK, or a testing
  framework for the core package.
- Use explicit result/status objects and diagnostics for numerical failures.
- Validate numerical behavior through deterministic analytical,
  manufactured, and independent Monte Carlo workflows where appropriate.
- Keep examples, tests, benchmarks, and validation tools as package consumers
  so installed-package behavior remains exercised.
- Keep user-facing APIs simple while allowing lower-level implementation
  diagnostics to support maintainers and focused tests.
