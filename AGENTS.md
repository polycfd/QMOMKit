# AGENTS.md

## Project and Priorities

QMOMKit is an application-agnostic C++ library for quadrature-based moment
methods, including one-dimensional QMOM, multidimensional and recursive CQMOM,
raw-moment initialization, realizability diagnostics, quadrature
reconstruction, and moment-transport utilities.

QMOMKit converts raw moments of a number density function into quadrature
representations that can be embedded in larger simulation, validation, or
analysis codes. It does not own meshes, cells, spatial discretizations, MPI
communication, physical boundary conditions, particle models, or
application-specific source terms.

Code clarity, algorithmic correctness, and computational efficiency are the
three guiding principles of QMOMKit and are of the utmost priority.

Prefer the simplest implementation that satisfies the mathematical,
performance, and API requirements. Do not introduce additional abstractions,
namespaces, directory layers, configuration mechanisms, or dependencies
without a concrete need.

## Documentation

Use the repository documentation as the source of truth before making changes.

```text
docs/user/
    Installation, configuration, examples, supported workflows, and public API
    usage.

docs/developer/
    Architecture, numerical conventions, algorithm design, transport,
    validation, testing, benchmarking, and engineering practices.

examples/
    Standalone applications consuming an installed QMOMKit package.

tests/
    Standalone unit, regression, verification, public-header,
    package-consumer, and benchmark projects.

tools/monte_carlo/
    Optional QMOMKitMonteCarlo companion package for finite-sample reference
    and validation workflows.

tools/plotting/
    Lightweight plotting helpers for QMOMKit output.
```

Do not duplicate detailed algorithm descriptions in this file. Before changing
a subsystem, read its relevant developer documentation and nearby tests.

Documentation must describe the current design rather than development history.
Do not preserve obsolete implementation stages, previous task instructions, or
superseded terminology.

## Package and Dependency Boundaries

The root CMake project builds and installs only the core library target:

```text
qmomkit::qmomkit
```

Tests, benchmarks, examples, plotting utilities, and Monte Carlo tools are
standalone consumers. They must not become dependencies of the root core
package.

The optional Monte Carlo companion package provides:

```text
qmomkit_monte_carlo::qmomkit_monte_carlo
```

It may depend on QMOMKit. The core QMOMKit package must not depend on Monte
Carlo headers, implementations, examples, scripts, or validation-only
infrastructure.

Keep mandatory dependencies minimal. Do not add a mandatory numerical backend,
test framework, application framework, or other third-party dependency without
an explicit design decision.

Google Benchmark is permitted only in the standalone benchmark project under
`tests/benchmarks/`. Benchmark infrastructure must not become a core runtime,
test, example, or Monte Carlo dependency.

## Configuration and Supported Workflows

QMOMKit configuration files use the `.qmom` extension. They describe numerical
configuration such as phase-space structure, moment initialization, node
counts, conditioning order, reconstruction options, and integration options.

Configuration files must not define application physics, own simulation loops,
or expose implementation-only storage and solver machinery.

Use the formulation appropriate to the phase-space problem:

```text
one phase-space dimension       -> QMOM
multiple phase-space dimensions -> CQMOM
```

The normal configured workflows are:

```text
load_qmom(...) -> QmomCase -> compute_qmom(...)
```

and

```text
load_cqmom(...) -> CqmomCase -> compute_cqmom(...) -> CqmomState
```

`CqmomCase` and `CqmomState` are the default representation-neutral CQMOM
interfaces. Dense or lean moment representations are advanced interfaces for
callers that intentionally manage their own moment storage or data transfer.

Do not expose parser internals, recursive storage details, or solver workspaces
through configuration or ordinary public workflows.

See `docs/developer/configuration_design.md` for configuration mapping and
extension rules.

## Public and Internal APIs

Installed public headers live under:

```text
include/qmomkit/
```

Implementation files and implementation-only helpers live under:

```text
src/
src/internal/
```

Public headers must not include headers from `src/internal/`.

Public APIs should expose stable mathematical and user-facing concepts,
including:

- requests and results;
- options and status values;
- diagnostics;
- moment containers;
- quadrature representations;
- transport state;
- documented algorithm entry points.

Do not expose parser internals, temporary solver workspaces, private
linear-algebra machinery, or recursive implementation details through public
headers.

Some mathematical or data-transfer APIs may intentionally be public for
advanced users. Such APIs should be clearly identified and documented as
advanced rather than disguised as ordinary workflow APIs.

Public API changes require deliberate review, focused tests, documentation
updates, and corresponding public-header coverage.

## Numerical Contracts

Changes to numerical conventions are high-risk. Preserve the contracts in
`docs/developer/numerical_conventions.md` unless a change is intentional,
documented, tested, and independently validated.

In particular:

- Moments are raw moments of the number density function.
- `m_0` and `M[0,...,0]` represent total number density and must not be
  normalized away implicitly.
- Phase-space dimension order and CQMOM conditioning order are distinct
  concepts.
- Moment indices, node counts, quadrature weights, abscissas, and recursive
  branch ordering must retain their documented semantics.
- Coordinate names are metadata. They must not impose physical units,
  supports, particle assumptions, velocity semantics, or application-specific
  behavior.
- Realizability assumptions, tolerances, fallback behavior, and status
  propagation must remain explicit.
- Data returned with a non-success status is diagnostic unless the public
  contract explicitly states otherwise.
- A successful numerical result must not depend on silently relaxed tolerances
  or suppressed failure conditions.

When a convention is ambiguous, make the assumption explicit and support it
with an analytical, manufactured, convergence, or independently generated
validation case.

## Transport and Time Integration

QMOMKit provides numerical operations for advancing moments and quadrature
states. Application physics enters through callbacks or caller-provided rate
evaluators.

Case and state objects carry numerical configuration and evolving numerical
state. They do not store user physical models.

Generic ODE utilities such as SSP-RK3 and adaptive RK23 must remain independent
of:

- particle size or shape;
- specific phase-space coordinate meanings;
- spatial meshes or cells;
- finite-volume or finite-element discretizations;
- application-specific source models.

QMOMKit does not implement general spatial PDE boundary conditions. Injection,
removal, nucleation, or other phase-space boundary contributions may be
represented through documented moment-rate callbacks.

Do not introduce the terms “hybrid transport,” “operator-split transport,” or
“Strang-split transport” in code or documentation.

See `docs/developer/transport_design.md` for the current transport workflow,
callback contracts, status propagation, and validation requirements.

## Implementation Principles

Use C++17 and follow the repository naming, formatting, and static-analysis
conventions.

Make changes that are small, focused, and reviewable. Do not perform unrelated
refactoring as part of a feature or bug fix.

Prefer:

- explicit control flow over hidden behavior;
- runtime-flexible algorithms over hard-coded phase-space dimensions;
- existing infrastructure over parallel abstractions;
- local implementation helpers over premature shared utilities;
- clear data ownership over mutable global state or implicit caches;
- straightforward C++ over excessive template metaprogramming.

Add a shared internal utility only when it removes meaningful duplication
across implementation modules and has a clear owner.

Public algorithm behavior should use explicit option, result, status, and
diagnostic types where these improve the caller contract. Failure paths must
remain visible and checkable.

Do not weaken mathematical checks merely to make a test or example pass.
Investigate the numerical cause and preserve useful diagnostics.

Modified C++ files must be formatted with the repository `.clang-format`. New
code must remain compatible with the repository `.clang-tidy` configuration.

Comments should explain mathematical intent, non-obvious assumptions, and
local constraints. They should not record development history or restate
obvious code.

## Testing and Validation

Every behavioral change requires tests at the level that owns the behavior.

Use:

```text
tests/unit/
    Local functions, types, options, indexing, diagnostics, validation, and
    failure paths.

tests/component_regression/
    Component interactions, status propagation, and regression-prone behavior.

tests/verification/
    Analytical, manufactured, reconstruction, convergence, and workflow
    verification.

tests/public_headers/
    Self-contained and installable public-header coverage.

tests/package_consumer/
    Installed-package discovery, exported targets, includes, and linkage.

tests/benchmarks/
    Performance measurement only.
```

Numerical algorithm changes require evidence beyond agreement with the previous
implementation. Validate them against analytical results, manufactured cases,
convergence behavior, independent implementations, or independently generated
reference data as appropriate.

Do not use benchmark timings as correctness evidence.

Do not hide discrepancies by silently loosening tolerances. Changes to
tolerances, acceptance criteria, node reduction, fallback behavior, node
ordering, or status mapping require focused tests and documented justification.

Validation changes must state:

- how the reference is constructed;
- what assumptions are made;
- which error metric is used;
- why the tolerance is appropriate;
- how the result can be reproduced.

See:

```text
docs/developer/testing.md
docs/developer/validation.md
```

## Build and Installation

The root project owns and installs the core QMOMKit package.

Tests, examples, benchmarks, package-consumer checks, and companion tools
should exercise the installed package rather than relying on private root-build
state.

Tests are plain C++ executables registered with CTest. Do not introduce a test
framework dependency unless explicitly requested and documented.

Changes to CMake, installation rules, exported targets, package metadata, or
public headers must preserve the installed downstream-consumer experience.

## Common Pitfalls

Avoid:

- normalizing raw moments by `m_0` without an explicit mathematical reason;
- confusing phase-space dimension order with conditioning order;
- assigning physical meaning based only on coordinate names;
- hard-coding dimensions, coordinate positions, or node counts into general
  algorithms;
- exposing internal implementation machinery as ordinary public API;
- coupling the core package to tests, examples, Monte Carlo tools, plotting
  scripts, or benchmarks;
- introducing CFD meshes, cells, MPI operations, spatial discretizations, or
  application physics into the core library;
- treating non-empty output from a failed result as a successful
  reconstruction;
- changing tolerances, fallback behavior, node ordering, status behavior, or
  realizability rules without focused tests and documentation;
- adding abstractions or directory structure solely for organizational
  symmetry.

## AI-Assisted Development Workflow

Before making substantial changes:

1. Read the relevant developer documentation and nearby tests.
2. Identify the owning module, package boundary, and public API implications.
3. Inspect existing patterns and reuse them where appropriate.
4. State the intended implementation approach.

While making changes:

- remain within the requested scope;
- preserve public/internal and core/consumer boundaries;
- keep mathematical assumptions explicit;
- add or update tests at the same ownership level as the behavior;
- update user or developer documentation when the supported design changes;
- run the relevant formatting, build, test, and validation commands.

When reporting completed work, include:

- files modified;
- behavior or documentation changed;
- tests and validation commands run;
- any unresolved concerns, limitations, or uncertainties.