# AGENTS.md

## Project Overview

QMOMKit is a C++ library for quadrature-based moment methods, including
one-dimensional QMOM, multidimensional CQMOM, recursive CQMOM, analytical
raw-moment initialization, realizability diagnostics, and moment-transport
components.

QMOMKit turns raw moments of a number density function into quadrature
representations that can be embedded in larger simulation, validation, or
analysis codes. It is CFD-code agnostic and application agnostic: it does not
own meshes, cells, finite-volume or finite-element operators, MPI
communication, boundary conditions, particle physics, or application-specific
source models.

Correctness, reproducibility, maintainability, and explicit mathematical
conventions take priority over premature optimization.

## Documentation Map

Use the documentation system as the source of truth before making changes.

```text
docs/user/
    User workflows, installation, examples, configuration, and public API use.

docs/developer/
    Architecture, numerical conventions, algorithm design, validation,
    testing, benchmark policy, and engineering practices.

examples/
    Standalone example applications that consume installed QMOMKit packages.

tests/
    Standalone test projects for unit, regression, verification, public-header,
    and package-consumer coverage.

tools/monte_carlo/
    Optional QMOMKitMonteCarlo companion package for finite-sample reference
    and validation workflows.

tools/plotting/
    Lightweight Python helpers for plotting QMOMKit-native output files.

AGENTS.md
    AI-assisted development rules for understanding and modifying QMOMKit
    safely.
```

Do not duplicate detailed algorithm descriptions here. When a change depends
on a specific subsystem, read the relevant `docs/developer/` page first.

## Architecture and Package Boundaries

The repository is organized around the core QMOMKit package and standalone
consumers. The root CMake project builds and installs only the core library
package:

```text
qmomkit::qmomkit
```

Tests, benchmarks, examples, and Monte Carlo tools are standalone consumers of
installed packages. They must not become dependencies of the root core build.

The optional Monte Carlo validation tools live under `tools/monte_carlo/` as a
separate companion package:

```text
qmomkit_monte_carlo::qmomkit_monte_carlo
```

The companion package may depend on QMOMKit, but the core package must not
depend on Monte Carlo headers, sources, examples, scripts, or validation-only
machinery.

The plotting helpers under `tools/plotting/` are lightweight scripts, not a
core library dependency or package surface.

For the current package layout, exported targets, consumer model, and build
philosophy, see `docs/developer/architecture.md`.

## Configuration and Workflow Principles

QMOMKit configuration files use the `.qmom` extension. They describe
phase-space structure, initial raw moments or analytical initialization data,
and numerical options. They do not define application physics, own simulation
loops, or expose implementation-only storage choices.

Formulation selection follows the phase-space problem:

```text
one phase-space dimension      -> QMOM
multiple phase-space dimensions -> CQMOM
```

The normal configured QMOM workflow is:

```text
load_qmom(...) -> QmomCase -> compute_qmom(...)
```

The normal configured CQMOM workflow is:

```text
load_cqmom(...) -> CqmomCase -> compute_cqmom(...) -> CqmomState -> transport
```

`CqmomCase` and `CqmomState` are the default representation-neutral CQMOM
objects. Dense and lean CQMOM workflows are explicit advanced representations
for callers that intentionally manage `RawMomentTensor` or required
`CqmomMomentData`. Do not add configuration fields that leak parser internals,
recursive storage machinery, or implementation-only solver workspaces.

For configuration mapping and extension rules, see
`docs/developer/configuration_design.md`.

## Public/Internal API Rules

Installed public headers live under:

```text
include/qmomkit/
```

Core implementation files and implementation-only helpers live under:

```text
src/
src/internal/
```

Public headers should expose stable user-facing concepts such as requests,
results, options, moment containers, quadrature representations, status values,
and documented algorithm entry points. Do not expose implementation-only
utilities, recursive assembly machinery, parser internals, solver workspaces,
or private dense-linear-algebra details through public headers.

Headers under `src/internal/` are not installed and are not public API. Public
headers under `include/qmomkit/` should not include internal headers.

Public API changes require deliberate review, matching tests, and documentation
updates. Update public-header tests when public headers are added, renamed, or
given new direct dependencies.

## Numerical Conventions and Non-Negotiables

Numerical convention changes are high-risk. Preserve the conventions in
`docs/developer/numerical_conventions.md` unless the change is intentional,
documented, tested, and validated.

In particular:

- QMOMKit moments are raw moments of the number density function.
- `m_0` and `M[0,...,0]` represent total number density and are not
  normalized away.
- One-dimensional moments, tensor multi-indices, phase-space dimensions,
  conditioning order, node counts, quadrature weights, and abscissas must keep
  their documented ordering semantics.
- Coordinate names are metadata and must not imply physical units, supports,
  velocity meanings, particle assumptions, or application-specific behavior.
- Realizability assumptions and status behavior must remain explicit.
- Non-success result data is diagnostic unless a public contract says
  otherwise.

When a mathematical convention is ambiguous, document the assumption, add or
update tests, and point to the relevant analytical, manufactured, or validation
case.

## Transport and Time Integration Boundaries

QMOMKit supports hybrid moment transport: continuous motion is represented by
advancing quadrature nodes, while discontinuous or source-like effects are
represented by advancing raw moments.

Runtime physics enters through callbacks such as `PhaseSpaceVelocity`,
`MomentSource`, `CqmomMomentSource`, `MomentBoundaryContribution`, and
`CqmomMomentBoundaryContribution`. Case and State objects carry numerical setup
and evolving moment representation; they do not store user physical models.

`TransportScheme::strang` is the default second-order composition. `cd` and
`dc` are supported first-order alternatives for comparisons, diagnostics, or
workflows that deliberately choose a fixed operator order.

SSP-RK3 and adaptive RK23 are generic vector ODE utilities used by transport
components. They must remain independent of particle shape, particle size,
velocity dimensions, CFD meshes, and application physics.

QMOMKit does not implement general PDE boundary conditions. Phase-space
boundary exchange, such as injection, removal, nucleation, or similar moment
rates, is represented through boundary contribution callbacks.

For transport design, status propagation, callback behavior, and verification
expectations, see `docs/developer/transport_design.md`.

## Code Modification Principles

Make small, focused, reviewable changes that fit the existing architecture.
Prefer simple runtime-flexible abstractions over hard-coded dimensions,
hard-coded node counts, hidden global state, mutable caches, or excessive
template metaprogramming.

Keep implementation code close to the module that owns the behavior. Add shared
internal utilities only when they remove real duplication across implementation
modules and do not weaken public/internal boundaries.

Use explicit request, result, option, status, and diagnostic types for public
algorithm behavior. Keep failure paths visible and checkable by callers.

Follow the repository naming, formatting, and static-analysis conventions.
Modified C++ files should be formatted with the repository `.clang-format`, and
new code should remain compatible with `.clang-tidy`.

Comments should explain mathematical intent, non-obvious assumptions, and local
constraints. Do not use comments or documentation to preserve obsolete task
history.

## Testing and Validation Expectations

Every code change should have tests at the level where the behavior is owned.
Use unit tests for local contracts, option validation, indexing, status values,
diagnostics, and failure paths.

Choose the test location by purpose:

- `tests/unit/`: isolated function, type, option, diagnostic, indexing, and
  failure-path contracts.
- `tests/component_regression/`: regression-prone component interactions,
  status propagation, and implementation behavior that should remain visible.
- `tests/verification/`: analytical, manufactured, convergence,
  reconstruction, or workflow checks against independent references.
- `tests/public_headers/`: compile-time coverage that public headers are
  self-contained and installable.
- `tests/package_consumer/`: installed-package discovery, exported target
  linkage, and downstream include behavior.

Numerical algorithm changes require verification against analytical,
manufactured, convergence, reconstruction, or independently generated
references where appropriate. A regression that only matches the previous
implementation is not enough evidence for a scientific algorithm change.

Public API changes should preserve public-header coverage. Package, export, or
installation changes should preserve installed-package consumer behavior.

Validation changes should update the validation contract: reference
construction, tolerances, error metrics, assumptions, and reproducibility data.
Do not hide discrepancies by silently relaxing tolerances.

For current testing and validation practice, see:

```text
docs/developer/testing.md
docs/developer/validation.md
```

## Documentation Expectations

Documentation should describe the current design, not development history.

Update `docs/user/` when supported user workflows, examples, configuration
semantics, installation behavior, or public API usage changes.

Update `docs/developer/` when architecture, public/internal boundaries,
numerical conventions, algorithms, transport behavior, testing strategy,
validation strategy, benchmark policy, or coding practices change.

Public interfaces should have concise comments that identify inputs, outputs,
status behavior, and important assumptions. Longer explanations belong in the
developer documentation. Local source comments should explain nearby reasoning,
not past tasks or obsolete change steps.

## Build and Dependency Principles

The core QMOMKit package keeps mandatory dependencies minimal. Do not add
mandatory third-party numerical backends or application-framework dependencies
without an explicit design decision.

The root build owns the core `qmomkit::qmomkit` package. Tests, benchmarks,
examples, and companion tools are standalone consumers of installed packages.
Preserve that model when changing CMake or package metadata.

Root builds install the core package. Tests, examples, benchmarks,
package-consumer checks, and companion tools should validate the installed
package experience rather than becoming part of the root core build.

Tests are plain C++ executables registered with CTest. Do not introduce a test
framework dependency unless explicitly requested and documented.

Google Benchmark is allowed only in the standalone benchmark project under
`tests/benchmarks/`. It must not become a core runtime, test, example, or
Monte Carlo dependency.

## Common Pitfalls

Avoid these mistakes:

- Normalizing moments by `m_0` unless the mathematical problem explicitly uses
  unit number density.
- Confusing phase-space dimension order with CQMOM conditioning order.
- Treating coordinate names as physical semantics.
- Treating internal headers, repository tests, examples, or benchmark helpers
  as public API.
- Coupling the core library to Monte Carlo validation tools.
- Adding CFD concepts, meshes, cells, boundary conditions, MPI behavior, or
  application-specific physics to the core package.
- Using benchmark timings as correctness or validation evidence.
- Treating non-empty output in a failed result as a successful numerical
  reconstruction.
- Changing tolerances, fallback behavior, node ordering, status mapping, or
  realizability behavior without focused tests and documentation.

## AI-Assisted Workflow

Before modifying code or documentation:

1. Read the relevant `docs/developer/` design page and nearby tests.
2. Identify the owning package boundary and public/internal API surface.
3. State the intended approach for substantial changes.

While making changes:

- Keep edits focused on the requested behavior.
- Preserve public/internal and core/consumer boundaries.
- Prefer existing patterns over new abstractions.
- Add or update tests and documentation at the same ownership level as the
  change.
- Ask for clarification when the mathematical convention or architecture
  boundary is genuinely ambiguous.

When reporting completed work, summarize:

- files modified;
- behavior or documentation changed;
- tests or validation commands run;
- remaining concerns, limitations, or uncertainties.
