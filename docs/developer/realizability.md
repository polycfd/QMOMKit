# Realizability

This page documents the current QMOMKit design for realizability handling. It is
a developer reference for maintainers, contributors, and AI-assisted development
workflows.

For raw-moment, number-density, and storage conventions, see
`numerical_conventions.md`. For source updates and phase-space evolution composition, see
`phase_space_evolution_design.md`. For configuration mapping of runtime policies, see
`configuration_design.md`.

## Role Of Realizability

Realizability is the question of whether a truncated raw moment sequence can
come from a nonnegative number density function on an allowed support. QMOMKit
moments are raw moments of the number density function:

```text
m_k = integral xi^k n(xi) dxi
```

with `m_0` equal to the total number density. A vector of finite numbers is not
automatically an admissible moment sequence. For example, a one-dimensional
sequence can have nonnegative `m_0` while implying an impossible variance or a
mean outside the allowed support.

Realizability matters because QMOM and CQMOM reconstruct quadrature
representations from moments. Reconstruction algorithms assume the moment data
are mathematically consistent with an admissible distribution and with the
requested quadrature. If moment evolution, initialization, or external data
produce an inadmissible sequence, inversion may fail or a candidate source-step
state should be rejected before it becomes the next runtime state.

QMOMKit treats realizability as an explicit numerical contract. The library
provides opt-in checks and policy adapters; it does not silently repair moment
sets, renormalize moments, or infer supports from coordinate names.

## Supported Moment Domains

Current explicit realizability checks are one-dimensional. They evaluate
`MomentSet` sequences against a selected support:

```text
real_line
nonnegative
bounded_interval
```

`MomentSupport::real_line` applies when the represented coordinate may take any
real value. Typical examples are unconstrained velocity-like coordinates or
other signed scalar variables.

`MomentSupport::nonnegative` applies when the coordinate must be zero or
positive. Typical examples are particle size, mass, count-like quantities, or
other variables whose support is `[0, infinity)`.

`MomentSupport::bounded_interval` applies when the coordinate must lie within a
closed interval. The public check uses an explicit `BoundedInterval` with finite
bounds and requires `lower < upper`. This support is appropriate for
interval-constrained quantities such as bounded normalized coordinates or volume
fractions.

Support is semantic input from the caller. It is not inferred from distribution
family, phase-space dimension name, output label, or units.

## Mathematical Criteria

QMOMKit's explicit realizability checks are based on classical
one-dimensional truncated moment problems. The implementation uses these
criteria to decide whether a finite sequence

```text
m_0, m_1, ..., m_N
```

can represent a nonnegative measure on the caller-selected support. The criteria
are support-aware: the same numeric sequence may be admissible on one support
and inadmissible on another.

For `MomentSupport::real_line`, the relevant model is the Hamburger moment
problem. The coordinate may take any real value, so the primary requirement is
that the moment sequence produce nonnegative quadratic forms for polynomials on
the real line. In finite truncated form, this is expressed through positive
semidefinite moment conditions and consistency checks when the moment data lie
on the boundary of moment space.

For `MomentSupport::nonnegative`, the relevant model is the Stieltjes moment
problem. The moment sequence must satisfy the real-line positivity requirements
and the additional constraints implied by support on `[0, infinity)`. These
additional conditions reject moment sequences that may be realizable on the full
real line but require mass outside the nonnegative support.

For `MomentSupport::bounded_interval`, the relevant model is the Hausdorff
moment problem on an explicit interval `[lower, upper]`. The sequence must be
compatible with both endpoints, not merely with nonnegative density. The
bounded-interval check therefore uses interval-specific support-localizing
conditions that encode the lower and upper bounds supplied through
`BoundedInterval`.

Support-localizing criteria are the mathematical bridge between physical
support assumptions and realizability. A variable such as particle size,
volume fraction, or a bounded internal coordinate carries more information than
"there exists a nonnegative measure"; it also restricts where that measure may
live. `check_realizability(...)` applies the criteria selected by
`MomentSupport` or `BoundedInterval` and reports the resulting classification
through `RealizabilityResult`.

## Realizability Checks

The public check is:

```cpp
RealizabilityResult check_realizability(
    const MomentSet& moments,
    MomentSupport support,
    const RealizabilityOptions& options = {});

RealizabilityResult check_realizability(
    const MomentSet& moments,
    BoundedInterval interval,
    const RealizabilityOptions& options = {});
```

The input is interpreted as the inclusive one-dimensional raw moment sequence

```text
m_0, m_1, ..., m_N
```

The result classifies the sequence with `RealizabilityStatus`:

- `strictly_realizable`: accepted away from the boundary of moment space;
- `boundary`: accepted on the boundary of moment space;
- `not_realizable`: rejected for the selected support;
- `invalid_input`: the request itself is invalid;
- `numerically_indeterminate`: the sequence is too close to the numerical
  tolerance to classify confidently.

`RealizabilityResult::realizable()` is true for `strictly_realizable` and
`boundary`. `RealizabilityResult::on_boundary()` is true only for `boundary`.

Diagnostics include the highest checked order, the limiting or failed matrix
name and order, a minimum pivot/eigenvalue-like diagnostic, the tolerance used,
and a human-readable message. Callers should treat the status as authoritative
and use the diagnostics for reporting, debugging, or adaptive failure handling.
Diagnostics are not a request to continue with invalid moments.

Invalid input means the request itself is malformed or unsupported, such as an
empty moment sequence, non-finite moment values, non-positive or non-finite
tolerances, a bounded-interval request without explicit bounds, or invalid
bounded-interval limits.

`not_realizable` means the supplied numbers form a well-defined request but
cannot represent a valid nonnegative measure on the selected support. Examples
include negative `m_0`, zero total number density with nonzero higher moments,
or moment sequences that fail the support-specific realizability conditions. The
all-zero sequence is treated as a boundary-realizable zero measure.

`RealizabilityOptions` controls numerical tolerance. With relative tolerance
enabled, the effective threshold scales with the checked matrix magnitude. Near
the boundary, the implementation may report `numerically_indeterminate` rather
than accepting or rejecting with false precision.

## Realizability Policies

`RealizabilityPolicy` is an opt-in runtime policy for accepting or rejecting
generic candidate states. It is separate from the realizability check itself and
from the generic Runge-Kutta solvers.

The policy contains:

- `check_mode`;
- support or bounded interval;
- numerical check options;
- `accept_boundary`.

The supported check modes are:

- `None`: install no acceptance predicate and accept without conversion;
- `FinalOnly`: check completed candidate states;
- `EveryStage`: check intermediate SSP-RK3 stage states and final states.

Runtime policies operate through caller-provided conversion from a generic state
vector to `MomentSet`:

```text
std::vector<double> state -> MomentSet
```

This conversion is caller-owned because QMOMKit cannot infer which components of
an arbitrary state vector are a one-dimensional moment sequence. Converter
exceptions propagate through the predicate path; the policy adapters do not
hide them.

`accept_boundary` controls whether `boundary` is accepted. Strictly realizable
states are accepted. `not_realizable`, `invalid_input`, and
`numerically_indeterminate` are rejected.

## Interaction With Phase-Space Evolution

Phase-space evolution uses realizability policies through SSP-RK3 acceptance predicates for
moment-source evolution. A source update can propose candidate raw moments that
leave the admissible moment space. A policy adapter can reject that candidate so
the SSP-RK3 retry logic may reduce the timestep according to its configured
retry and minimum-timestep behavior.

The public adapters are:

```cpp
make_realizability_acceptance_predicate(...)
make_realizability_stage_acceptance_predicate(...)
```

`make_realizability_acceptance_predicate(...)` returns a final-candidate
predicate only for `FinalOnly` policies. `None` and `EveryStage` return an empty
final-candidate predicate.

`make_realizability_stage_acceptance_predicate(...)` returns:

- an empty predicate for `None`;
- a stage-aware predicate that checks only `SspRk3Stage::Final` for
  `FinalOnly`;
- a stage-aware predicate that checks every supplied SSP-RK3 stage for
  `EveryStage`.

Phase-space evolution does not install these policies automatically. Configuration parsing
can create a `RealizabilityPolicy`, but the application or driver must wire the
policy into the relevant time-integration options.

When a source update is rejected by a realizability predicate and retries cannot
produce an accepted state, the source-update layer reports predicate or
minimum-timestep failure. High-level phase-space evolution maps source-update failure to the
phase-space evolution source-failure path and records subsystem diagnostics.

## Interaction With QMOM And CQMOM

Realizability is a constraint on moment representations. It is not a replacement
for QMOM or CQMOM reconstruction status.

For one-dimensional QMOM, explicit checks can be used before calling
`compute_qmom(...)` when the workflow needs support-aware validation. Wheeler
inversion still performs its own reconstruction and status handling for the
requested node count and options.

For CQMOM, explicit realizability checks are currently one-dimensional. They can
be applied to one-dimensional moment sequences or branch-local states when the
caller has such a representation. Full multidimensional recursive CQMOM
validity is reported through `compute_cqmom(...)` status and diagnostics.

Initialization workflows should distinguish between configuration validity,
moment realizability, and reconstruction success:

- configuration parsing validates syntax and structural consistency;
- distribution helpers generate raw moments from supported analytical models;
- realizability checks classify one-dimensional moment sequences on explicit
  supports;
- QMOM and CQMOM reconstruction report whether the requested quadrature was
  accepted.

Each layer reports its own status. A successful earlier layer does not guarantee
success in a later one.

## Failure Handling

Realizability APIs use explicit result objects rather than exceptions for normal
classification failure. Callers should check `RealizabilityResult::status` or
`realizable()` before using a moment sequence as accepted input under a support
policy.

Failure handling rules are:

- `invalid_input` means the check request is malformed or unsupported;
- `not_realizable` means the moments are incompatible with the selected support;
- `numerically_indeterminate` means the sequence is too close to the configured
  tolerance boundary and should not be treated as accepted by runtime policies;
- `boundary` is mathematically accepted but may require caller policy;
- `strictly_realizable` is accepted by the check.

QMOMKit does not automatically correct failed sequences. It does not alter the
moments, adjust supports, relax tolerances, or reduce quadrature orders as part
of explicit realizability checking. Retry behavior exists only where the check
is used as an acceptance predicate inside a timestepper that has its own retry
policy.

## Public API Boundary

The public realizability API lives in:

```text
include/qmomkit/realizability.hpp
include/qmomkit/realizability_policy.hpp
```

Public concepts include:

- `check_realizability(...)`;
- `MomentSupport`;
- `BoundedInterval`;
- `RealizabilityOptions`;
- `RealizabilityStatus`;
- `RealizabilityResult`;
- `RealizabilityPolicy`;
- `RealizabilityCheckMode`;
- `accepts_realizability_policy(...)`;
- `make_realizability_acceptance_predicate(...)`;
- `make_realizability_stage_acceptance_predicate(...)`.

Implementation details include:

- support-specific matrix construction;
- eigenvalue and kernel-consistency analysis;
- tolerance scaling internals;
- internal dense linear-algebra utilities;
- helper structures used to classify limiting matrices.

These implementation details support the public contracts but are not
downstream extension points.

## Extension Principles

Realizability changes should preserve these rules:

- keep support assumptions explicit;
- keep one-dimensional explicit checks separate from full CQMOM reconstruction
  diagnostics;
- use result/status objects for classification and diagnostics;
- do not silently normalize, repair, or reinterpret moments;
- reject numerically indeterminate states in runtime acceptance policies;
- keep runtime policies opt-in and caller-wired;
- validate new checks against analytical, manufactured, or independently
  generated moment cases;
- document any new support domain, status, or policy behavior before exposing it
  through public configuration.

The goal is conservative numerical validation: realizability support should help
developers and drivers detect invalid moment states early without hiding the
source of the failure.
