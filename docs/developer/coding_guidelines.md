# Coding Guidelines

This page documents QMOMKit-specific development principles. It complements the
public API boundaries in [architecture.md](architecture.md), the mathematical
conventions in [numerical_conventions.md](numerical_conventions.md), and the
correctness, validation, and performance practices in [testing.md](testing.md),
[validation.md](validation.md), and [benchmarks.md](benchmarks.md).

These guidelines describe stable engineering expectations for QMOMKit. They are
not a generic C++ style guide, a contribution workflow, or an instruction file
for automated development tools.

## Project Philosophy

QMOMKit is a reusable C++ library for QMOM, CQMOM, and
phase-space evolution algorithms. The core package remains independent of CFD
frameworks, physical models, application domains, meshes, boundary conditions,
and parallel-runtime assumptions.

Development should preserve these priorities:

- numerical correctness, reproducibility, maintainability, and extensibility
  take priority over premature optimization;
- scientific abstractions should be explicit and connected to the moment-method
  mathematics they represent;
- algorithms should support arbitrary phase-space dimensions and quadrature-node
  hierarchies where the public model requires that flexibility;
- physical models and application-specific choices belong to callers, examples,
  or companion tools rather than the core library;
- package boundaries should stay explicit so QMOMKit can be consumed as an
  installed library.

Clear, conservative designs are preferred over broad abstractions that obscure
the mathematical contract.

## Public API Discipline

Public QMOMKit API lives under `include/qmomkit/`. Installed headers represent
supported functionality and should expose stable concepts such as requests,
results, options, moment containers, quadrature representations, and documented
algorithm entry points.

Implementation code belongs in `src/`. Shared private utilities belong in
`src/internal/`; those headers are not installed, are not package API, and
should not be included by public headers. Internal utilities should be shared
only when they remove real duplication between implementation modules. They
should not become a general holding area for algorithm-specific code.

Public APIs require deliberate design:

- expose mathematical inputs and outputs rather than incidental implementation
  structures;
- keep ownership, status, and failure behavior explicit;
- keep public headers self-contained and aligned with the installed package
  model;
- avoid leaking recursive assembly machinery, solver internals, parser details,
  or other private helpers into user-facing interfaces;
- update public-header tests when supported headers or public API dependencies
  change.

The QMOMKitMonteCarlo companion package is a separate consumer-oriented
package. It may depend on `qmomkit::qmomkit`, but the core package should not
depend on Monte Carlo headers or implementation details.

## Numerical Code Principles

Numerical code must preserve the conventions documented in
[numerical_conventions.md](numerical_conventions.md). In particular, QMOMKit
moments are raw moments of the number density function, with `m0` representing
total number density. Moment ordering, tensor indexing, conditioning order,
quadrature weights, and abscissa ordering should not change silently.

Developer-facing numerical expectations are:

- document mathematical assumptions, domains, indexing rules, and failure modes;
- keep status paths explicit and checkable by callers;
- prefer deterministic behavior where the algorithm allows it;
- treat fallback data returned with non-success statuses as diagnostic unless a
  public contract states otherwise;
- avoid hidden mutable caches, global state, hard-coded dimensions, and
  hard-coded node counts;
- avoid application-specific terminology in core algorithms;
- preserve realizability and conditioning diagnostics near the boundary of
  moment space;
- keep solver behavior independent of particle shape, size, velocity model, or
  other caller-owned physics.

Silent convention changes are especially risky in moment methods. A change that
alters moment interpretation, quadrature representation, conditioning order, or
realizability behavior should be reflected in developer documentation and
covered by focused verification.

## Assertions and Runtime Checks

QMOMKit distinguishes between runtime validation and assertions. Runtime checks
protect public interfaces, caller-provided data, and numerical correctness.
Assertions document internal invariants that should always hold after
validation and object construction.

Use runtime validation for conditions that can arise from valid external usage
or caller-controlled workflows, including:

- user-provided inputs and public API requests;
- configuration values and option combinations;
- inconsistent moment sets or tensor extents;
- invalid distribution parameters;
- non-finite numerical data;
- recoverable numerical failures such as non-realizability, singular systems,
  or rejected integration candidates.

Runtime checks should provide meaningful status values, diagnostics,
exceptions, or error messages according to the API contract. They should not be
removed only because they appear in numerical code; many validation branches
preserve documented failure behavior and protect callers from interpreting
diagnostic data as successful output.

Use assertions for internal invariants that should be guaranteed by earlier
validation, construction, or algorithm logic, including:

- internal indexing assumptions;
- validated relationships between object sizes, extents, and storage layouts;
- impossible states during algorithm execution;
- internal storage consistency in private implementation utilities.

Assertions should not replace checks that protect public API contracts or
conditions that can occur through supported external usage. In particular,
intermediate numerical construction can pass through states that are valid for
the algorithm even when they would be invalid as final public output.

Performance changes to validation should be justified by profiling,
benchmarking, or clear evidence that the check runs in a computational
bottleneck. Numerical safety checks should not be removed from hot paths
without considering the resulting failure behavior. Prefer maintainable
improvements such as moving validation to an appropriate boundary, avoiding
repeated validation of unchanged data, or making expensive diagnostics optional
where the public contract allows it.

## Testing and Validation Expectations

Testing, validation, and benchmarking serve different purposes:

- tests protect software correctness and regression behavior;
- validation establishes numerical and scientific credibility;
- benchmarks measure computational performance.

Changes to local behavior should be covered by unit tests. Changes to numerical
algorithms should also have verification coverage against analytical,
manufactured, or independently generated references where appropriate. Public
API changes should preserve public-header coverage. Package and export changes
should preserve installed-package consumer validation.

Scientific changes require more than matching a previous regression result.
When a method, tolerance, status rule, reconstruction path, phase-space
evolution update, or realizability decision changes the numerical behavior, the
relevant validation assumptions should be updated and documented. Numerical
discrepancies should be explained rather than hidden by undocumented tolerance
changes.

Benchmarks should remain separate from correctness tests. They are useful for
tracking representative performance, but they do not establish mathematical
correctness.

## Documentation Expectations

Documentation should describe the current design. Avoid stage labels, migration
notes, optimization history, temporary limitations, and roadmap language in
developer references.

User-facing behavior belongs in `docs/user/`. Knowledge and guidelines for
development belongs in `docs/developer/`. Local implementation details may be
documented in source comments when they help explain mathematical intent,
assumptions, or constraints.

Public interfaces should have concise comments that identify inputs, outputs,
status behavior, and important assumptions. Longer explanations should live in
developer documentation, especially for:

- mathematical conventions;
- algorithm design;
- indexing and ordering rules;
- validation methodology;
- failure handling;
- package and API boundaries.

User documentation should be updated when supported workflows, configuration
semantics, examples, or public behavior change. Developer documentation should
be updated when architecture, numerical conventions, algorithms, testing,
validation, benchmark methodology, or coding principles change.

## Dependency and Build Philosophy

The core library keeps runtime dependencies minimal. The root package builds
and installs the core `qmomkit::qmomkit` target and its public headers. Tests,
benchmarks, examples, and companion tools are standalone consumers of installed
packages rather than dependencies of the core library build.

Dependency choices should preserve the package model:

- keep the core library independent of CFD frameworks and application-specific
  runtimes;
- avoid mandatory third-party numerical backends unless the public design
  explicitly adopts them;
- keep benchmark dependencies local to `tests/benchmarks/`;
- keep test-framework choices separate from the runtime package;
- avoid coupling examples or Monte Carlo tooling back into core algorithms.

QMOMKit currently uses C++17 for the public target. Build-related changes
should preserve the installed-package workflow and exported CMake target
behavior.

## Code Quality Expectations

QMOMKit code should favor readability, simplicity, explicit behavior, and
conservative modern C++. The repository formatting and static-analysis settings
define the mechanical baseline for C++ changes, and compiler warnings should be
treated as signals to improve clarity rather than as noise to suppress.

Project-specific practices include:

- use the `qmomkit` namespace for core library code;
- follow existing naming patterns for classes, functions, variables, and private
  members;
- prefer named request, result, option, and status types over ambiguous argument
  lists for public algorithms;
- keep comments focused on mathematical intent, non-obvious assumptions, and
  constraints;
- isolate dense linear-algebra utilities behind internal implementation
  boundaries;
- prefer runtime flexibility over template-heavy designs unless compile-time
  structure is part of the public mathematical model;
- keep behavior reproducible across tests and examples where possible.

The best QMOMKit code makes the numerical contract easy to inspect: the data
representation is clear, the assumptions are named, failure behavior is visible,
and the implementation boundary matches the package boundary.
