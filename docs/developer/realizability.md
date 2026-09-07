# Realizability

This page documents the current QMOMKit design for realizability handling. 

For raw-moment, number-density, and storage conventions, see
[numerical_conventions.md](numerical_conventions.md). For source updates and phase-space evolution composition, see [phase_space_evolution_design.md](phase_space_evolution_design.md). For configuration mapping of runtime policies, see [configuration_design.md](configuration_design.md).

## Role Of Realizability

Realizability is the question of whether a truncated raw moment sequence can
come from a nonnegative number density function on an allowed support. QMOMKit
moments are raw moments of the number density function:
$$m_k = \int \xi^k n(\xi) \ \mathrm{d}\xi$$


with $m_0$ equal to the total number density. A vector of finite numbers is not
automatically an admissible moment sequence. For example, a one-dimensional
sequence can have nonnegative $m_0$ while implying an impossible variance or a
mean outside the allowed support.

Realizability and quadrature reconstruction are separate mathematical
operations. Realizability asks whether a moment sequence can represent a valid
distribution on an explicit support. Reconstruction asks whether QMOM or CQMOM
can construct a finite quadrature representation matching the supplied
moments. A realizability check does not validate Wheeler inversion, recursive
CQMOM reconstruction, node counts, or quadrature topology. Reconstruction
status does not replace a support-aware realizability check.

QMOMKit treats realizability as an explicit numerical contract. The library
provides opt-in checks; it does not silently repair moment sets, renormalize
moments, or infer supports from coordinate names.

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
enabled, the effective threshold, roundoff allowance, and numerical matrix
analysis scale with the checked matrix magnitude. Positive number-density
scaling therefore does not by itself change the classification within
floating-point limits. With relative tolerance disabled, `tolerance` is an
absolute threshold. Near the boundary, the implementation may report
`numerically_indeterminate` rather than accepting or rejecting with false
precision.

## Runtime Checks In Phase-Space Evolution

Runtime realizability checks apply only to transported raw moments during
moment-source SSP-RK3 evolution. They are independent of node advection:
neither QMOM nor CQMOM node trajectories, quadrature evolution, nor
reconstructed nodes and weights are realizability-check inputs.

For QMOM, `[time_integration.realizability]` configures the check. The setting
is part of `QmomPhaseSpaceEvolutionOptions` when the phase-space evolution
options are loaded from the configuration file. `FinalOnly` checks after the
complete moment-source integration step; `EveryStage` checks after every
moment-source SSP-RK3 stage. `None` performs no runtime check.

For CQMOM, realizability settings are configured separately for individual
phase-space coordinates. Each configured setting checks exactly the
one-dimensional marginal sequence for coordinate `d`:
$$M[0,\ldots,k_d,\ldots,0].$$

Only coordinates with explicit settings are checked. Mixed moments and
conditional moments are not checked. CQMOM quadrature reconstruction, adaptive
topology, and branch counts are also outside realizability checking. This is an
intentional separation between coordinate-marginal moment validation and
recursive quadrature reconstruction.

`accept_boundary` determines whether a `boundary` classification is accepted.
Strictly realizable moments are accepted; `not_realizable`, `invalid_input`,
and `numerically_indeterminate` classifications are rejected. If a source
integration attempt is rejected and the configured retry controls cannot obtain
an accepted moment state, moment-source evolution reports failure through its
normal status and diagnostics.

## Interaction With QMOM And CQMOM

Realizability is a constraint on moment representations. It is not a replacement
for QMOM or CQMOM reconstruction status.

For one-dimensional QMOM, explicit checks can be used before calling
`compute_qmom(...)` when the workflow needs support-aware validation. Wheeler
inversion remains a separate reconstruction operation with its own status
handling for the requested node count and options.

For CQMOM, realizability checks are one-dimensional coordinate-marginal checks;
they do not establish full multidimensional moment validity. Recursive CQMOM
reconstruction reports the outcome of its own quadrature construction through
`compute_cqmom(...)` status and diagnostics.

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
`realizable()` before using a moment sequence as accepted input under the
selected support settings.

Failure handling rules are:

- `invalid_input` means the check request is malformed or unsupported;
- `not_realizable` means the moments are incompatible with the selected support;
- `numerically_indeterminate` means the sequence is too close to the configured
  tolerance boundary and should not be accepted by runtime checks;
- `boundary` is mathematically accepted but may require an explicit acceptance
  setting;
- `strictly_realizable` is accepted by the check.

QMOMKit does not automatically correct failed sequences. It does not alter the
moments, adjust supports, relax tolerances, or reduce quadrature orders as part
of explicit realizability checking. Retry behavior applies only to
moment-source integration with runtime checking enabled and configured retry
controls.

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
- reject numerically indeterminate states in runtime realizability checks;
- keep runtime realizability checks opt-in and scoped to moment-source
  evolution;
- validate new checks against analytical, manufactured, or independently
  generated moment cases;
- document any new support domain, status, or runtime-check behavior before
  exposing it through public configuration.

The goal is conservative numerical validation: realizability support should help
developers and drivers detect invalid moment states early without hiding the
source of the failure.