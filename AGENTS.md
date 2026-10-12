# AGENTS.md

## Project and Priorities

QMOMKit is an application-agnostic C++ library for quadrature-based moment
methods, including one-dimensional QMOM, multidimensional CQMOM,
raw-moment initialization, realizability diagnostics, quadrature
reconstruction, and phase-space evolution utilities.

QMOMKit converts raw moments of a number density function into quadrature
representations that can be embedded in larger simulation, validation, or
analysis codes. It does not own meshes, cells, spatial discretizations, MPI
communication, physical boundary conditions, particle models, or
application-specific source terms.

QMOMKit is serial by design. Anything that pertains to the host application is
the host's responsibility, not QMOMKit's. This includes parallelism and
scheduling, the time loop and the choice of the outer time step, deciding what
to do when a QMOMKit call returns a failure status (for example, repeating the
host time step with a smaller outer time step), mapping batch positions to cells
or other host identifiers, and coupling to host fields. Step-size control inside
a single QMOMKit call, namely adaptive node-advection substeps and moment-source
retries after a rejected candidate, belongs to QMOMKit and is configured through
its options. QMOMKit returns each result with its status; what to do with it is
the host's decision. Do not add host-side machinery to QMOMKit, and do not
describe its absence as a limitation of QMOMKit.

Prefer the simplest implementation that satisfies the mathematical,
performance, and API requirements. Do not introduce additional abstractions,
namespaces, directory layers, configuration mechanisms, or dependencies
without a concrete need.

**Code clarity, algorithmic correctness, and computational efficiency are the
three guiding principles of QMOMKit and are of the utmost priority!**

## Agent Operating Rules

Before editing, inspect staged, unstaged, and untracked changes. Distinguish
pre-existing work from task changes, preserve unrelated work, and do not assume
every modified path belongs to the current task. Never reset, restore, clean,
discard, overwrite, stage, or commit unrelated changes. If unexpected unrelated
changes affect the requested scope, report them rather than silently
incorporating them. Do not stage or commit.

When multiple repositories or adjacent workspaces are present, identify their
boundaries and treat their Git states separately. Do not assume similarly named
paths belong to one repository and do not modify adjacent repositories or
manuscript material unless explicitly requested.

Verify actual paths, public APIs, options, repository state, commands,
numerical results, and validation evidence. Distinguish observed facts from
assumptions. Never invent citations, commands, test outcomes, validation
results, API names, options, paths, or repository state.

Stop and report a precise blocker when the requested behavior is unsupported by
the public API or configuration, source reality contradicts a key assumption,
the task requires an unrequested redesign, or required information is genuinely
unavailable. Do not add an expedient local workaround merely to complete a
task.

## Documentation

Read the relevant repository documentation before making changes. Use it as
design guidance, and verify public behavior against current headers, source,
configuration, and tests when necessary. Report conflicts rather than assuming
that stale prose or a prompt is authoritative. Use primary literature or
independently established evidence for scientific claims.

```text
docs/user/
    Installation, configuration, examples, supported workflows, and public API
    usage.

docs/developer/
    Architecture, numerical conventions, algorithm design, phase-space evolution,
    validation, testing, and engineering practices.

examples/
    Standalone applications consuming an installed QMOMKit package.

tests/
    Standalone unit, regression, verification, public-header, and
    package-consumer projects.

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

Update only documentation relevant to the change. Use code comments for
non-obvious implementation invariants, developer documentation for algorithms,
architecture, and extension behavior, user documentation for public API,
configuration, assumptions, units, and usage, and example README files for
physical examples. Do not broadly rewrite unrelated documentation.

## Tests and Examples

Tests and examples are strictly separate and must never be conflated.

- Tests belong under `tests/` and provide automated unit, component,
  regression, verification, and public-header coverage.
- Examples belong under `examples/` and provide clear user-facing workflows
  demonstrating how to use QMOMKit.
- Examples must follow the native-QMOMKit-first rule below, remain small and
  explicit, and must not become local frameworks.
- An example must never be registered with CTest or treated as an automated
  test, including as a smoke test.
- Example CMake projects must not call `enable_testing()` or `add_test()`.
- Important invariants demonstrated by an example must have focused coverage
  under `tests/`.
- Examples may perform concise runtime checks and return failure when their
  demonstrated workflow is invalid, but those checks do not make the example a
  test.
- Tests must not invoke example executables or depend on example output.
- Run changed examples and check their numerical behavior; compilation alone is
  not sufficient. Document relevant physical assumptions and units, and preserve
  Monte Carlo or reference companions when they are part of the example design.
  Examples are user-facing documentation and application-level evidence, not
  replacements for focused tests or verification.

## Paper and Manuscript Boundary

`paper/examples/` and `paper/verification/` are paper-specific studies,
distinct from standard examples and the test suite. Do not modify manuscript
repositories, generated manuscript data, or manuscript text unless explicitly
requested.

## Package and Dependency Boundaries

The root CMake project builds and installs only the core library target:

```text
qmomkit::qmomkit
```

Tests, examples, plotting utilities, and Monte Carlo tools are standalone
consumers. They must not become dependencies of the root core package.

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
load_qmom(...)           -> QmomCase          -> compute_qmom(...)
                                              -> initial_state() -> QmomState
                                              -> state_candidate(moments) -> QmomStateCandidate

load_cqmom_prepared(...) -> CqmomPreparedCase -> compute_cqmom(...)
                                              -> initial_state() -> CqmomState
                                              -> state_candidate(moments) -> CqmomStateCandidate
```

`compute_qmom(...)` and `compute_cqmom(...)` reconstruct a quadrature;
`initial_state()` creates the state that phase-space evolution advances.
`state_candidate(...)` creates that state from moments the host changed
itself, for example by transport in physical space, optionally after a
realizability check, and reports a failure as a status instead of throwing.

`QmomCase` and `QmomState` are the QMOM interfaces. `CqmomPreparedCase` and
`CqmomState` are the default CQMOM interfaces. QMOM has no prepared
counterpart, because a one-dimensional reconstruction has no recursive branch
structure to build. A prepared case builds the recursive branch structure of
its hierarchy once and shares it with every reconstruction, every state it
creates, and every case derived from it with
`config::with_initial_moments(...)`. `CqmomCase` is the CQMOM counterpart of
`QmomCase`: an open case definition that callers build or edit in code and then
prepare. Its own `compute_cqmom(...)` and `initial_state()` rebuild the branch
structure on every call, so they suit a single reconstruction or state rather
than repeated ones.

Do not expose parser internals, recursive storage details, or solver workspaces
through configuration or ordinary public workflows.

If requested behavior is available through QMOMKit's public API or `.qmom`
configuration, use that mechanism. Do not duplicate, override, bypass, or
reimplement configurable QMOMKit behavior in examples, paper studies, tools, or
application code. If the supported public/configuration path cannot perform the
requested task, stop and report the precise blocker instead of adding a local
workaround.

See `docs/developer/configuration_and_output.md` for configuration mapping and
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
- phase-space evolution state;
- documented algorithm entry points.

Do not expose parser internals, temporary solver workspaces, private
linear-algebra machinery, or recursive implementation details through public
headers.

Some mathematical or data-transfer APIs may intentionally be public for
advanced users. Such APIs should be clearly identified and documented as
advanced rather than disguised as ordinary workflow APIs.

Public API changes require deliberate review, focused tests, documentation
updates, and corresponding public-header coverage. When a change affects public
or configurable behavior, consider public headers, defaults, configuration
parsing and backward compatibility where applicable, observable diagnostics and
output semantics, source/state/evolution behavior, existing examples, and
installed downstream consumers where relevant.

## Numerical Contracts

Changes to numerical conventions are high-risk. Preserve the contracts in
`docs/developer/numerical_conventions.md` unless a change is intentional,
documented, tested, and independently validated.

When a convention is ambiguous, make the assumption explicit and support it
with an analytical, manufactured, convergence, or independently generated
validation case.

The following conventions must not be changed:

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

## Phase-Space Evolution and Time Integration

QMOMKit provides numerical operations for advancing moments and quadrature
states. Application physics enters through callbacks or caller-provided rate
evaluators.

Case and state objects carry numerical configuration and evolving numerical
state. They do not store user physical models.

Generic ODE utilities, such as the implemented Runge-Kutta schemes (SSP-RK3, RK23, RK45) must remain independent of:
- particle size or shape;
- specific phase-space coordinate meanings;
- spatial meshes or cells;
- finite-volume or finite-element discretizations;
- application-specific source models.

QMOMKit does not implement general spatial PDE boundary conditions. Injection,
removal, nucleation, or other phase-space boundary contributions may be
represented through documented moment-rate callbacks.

### Phase-Space Evolution Terminology

Use established numerical-method terminology accurately.

QMOMKit applies operator splitting for the phase-space advection and
moment-source operators. Strang splitting is one of the implemented
splitting schemes.

Terms such as “operator splitting,” “Lie splitting,” “Strang splitting,”
“operator-split methods,” “Strang-split integration,” “splitting scheme,” and
“splitting order” are valid and should be used where mathematically
appropriate.

Avoid misleading descriptions such as “hybrid phase-space evolution” when referring to the
standard operator-splitting workflow.

See `docs/developer/phase_space_evolution.md` for the current phase-space evolution workflow,
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

Never weaken mathematical checks merely to make a test or example pass.
Investigate the numerical cause, report suspicious behavior, and preserve useful diagnostics.

Modified C++ files must be formatted with the repository `.clang-format`. New
code must remain compatible with the repository `.clang-tidy` configuration.
Format only changed or relevant files and avoid wholesale formatting churn.
Keep temporary builds and numerical probes outside the repository unless they
are intended committed artifacts. Run `git diff --check` before final reporting.

Comments should explain mathematical intent, non-obvious assumptions, and
local constraints. They should not record development history or restate
obvious code.

## Testing and Validation

Every behavioral change requires tests at the level that owns the behavior.

A genuine user-facing QMOMKit feature must include suitable unit, regression,
and/or verification coverage; a physically meaningful example using the public
API or `.qmom` configuration, or a justified extension of an existing example;
appropriate user documentation; developer documentation where appropriate; and
code-level explanation for important non-obvious invariants. A narrow bug fix,
documentation correction, or internal-only refactor does not automatically
require a new physical example.

### Debug Functional-Test Requirement

Unit tests under `tests/unit/` intentionally use C/C++ `assert()`, and `NDEBUG`
disables those checks. Build and run QMOMKit functional tests with
`CMAKE_BUILD_TYPE=Debug`. In general, use Debug by default when building,
running, auditing, or validating the suite. A Release CTest result is not proof
that the unit suite passed. Do not propose replacing the assertions unless the
user explicitly requests a redesign of the unit-test mechanism.

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

tools/monte_carlo/tests/
    Monte Carlo companion tests. They check without assert(), so any build
    type runs them in full.
```

Numerical algorithm changes must be supported by evidence beyond agreement with the
previous implementation. Validate them against analytical results, manufactured cases,
convergence behavior, independent implementations, or independently generated
reference data as appropriate.

Passing the pre-existing test suite is necessary but does not, by itself,
verify newly introduced numerical behavior. Support new numerical functionality
with evidence proportionate to the change, such as focused unit or regression
tests, analytical or manufactured verification, convergence evidence,
independent numerical comparison, or established literature results.

Do not use timings as correctness evidence.

Do not hide discrepancies by silently loosening tolerances. Changes to
tolerances, acceptance criteria, node reduction, fallback behavior, node
ordering, or status mapping require explicit approval, focused tests, and documented justification.

Validation changes must state:

- how the reference is constructed;
- what assumptions are made;
- which error metric is used;
- why the tolerance is appropriate;
- how the result can be reproduced.

See:

```text
docs/developer/testing_and_validation.md
```

## Build and Installation

The root project owns and installs the core QMOMKit package.

The normal root CTest suite runs in a Debug build against the build-tree
`qmomkit` target. Installed-package validation is separate:
`tests/package_consumer/` is the dedicated installed-package consumer check.
Examples and companion tools may also consume an installed package where
appropriate. Run installed-package validation when changing installation,
exports, package metadata, installed public headers, or public behavior whose
downstream consumption is relevant.

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
- coupling the core package to tests, examples, Monte Carlo tools, or plotting
  scripts;
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

## Copyright and Source Code Adaptation

Do not copy or adapt source code from external projects unless explicitly instructed and the source, revision, copyright holder, and licence are recorded. For named algorithms, implement from mathematical descriptions, standards, or papers rather than from third-party source implementations. Report every external source consulted. Do not remove copyright, attribution, licence, or NOTICE information. When uncertain whether code is derived, stop and report the uncertainty rather than incorporating it.
