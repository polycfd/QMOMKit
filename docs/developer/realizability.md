# Realizability

This page documents the current QMOMKit design for realizability handling.

For raw-moment, number-density, and storage conventions, see
[numerical_conventions.md](numerical_conventions.md). For source updates and
phase-space evolution composition, see
[phase_space_evolution_design.md](phase_space_evolution_design.md). For
configuration mapping of runtime policies, see
[configuration_design.md](configuration_design.md).

## Contents

- [Role Of Realizability](#role-of-realizability)
- [Supported Moment Domains](#supported-moment-domains)
- [Mathematical Criteria](#mathematical-criteria)
- [Realizability Checks](#realizability-checks)
- [Numerical Classification](#numerical-classification)
- [Runtime Checks In Phase-Space Evolution](#runtime-checks-in-phase-space-evolution)
- [Interaction With QMOM And CQMOM](#interaction-with-qmom-and-cqmom)
- [Failure Handling](#failure-handling)
- [Public API Boundary](#public-api-boundary)
- [Extension Principles](#extension-principles)

## Role of realizability

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

## Supported moment domains

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
closed interval. The public check uses an explicit `BoundedInterval` with
finite bounds and requires `lower < upper`. This support is appropriate for
interval-constrained quantities such as bounded normalized coordinates or
volume fractions.

Support is semantic input from the caller. It is not inferred from distribution
family, phase-space dimension name, output label, or units.

## Mathematical criteria

QMOMKit's explicit realizability checks are based on classical
one-dimensional truncated moment problems. The implementation uses these
criteria to decide whether a finite sequence

```text
m_0, m_1, ..., m_N
```

can represent a nonnegative measure on the caller-selected support. The
criteria are support-aware: the same numeric sequence may be admissible on one
support and inadmissible on another.

For `MomentSupport::real_line`, the relevant model is the Hamburger moment
problem. The coordinate may take any real value, so the primary requirement is
that the moment sequence produce nonnegative quadratic forms for polynomials on
the real line. In finite truncated form, this is expressed through positive
semidefinite moment conditions and consistency checks when the moment data lie
on the boundary of moment space.

For `MomentSupport::nonnegative`, the relevant model is the Stieltjes moment
problem. The moment sequence must satisfy the real-line positivity requirements
and the additional constraints implied by support on `[0, infinity)`. These
additional conditions reject moment sequences that may be realizable on the
full real line but require mass outside the nonnegative support.

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

## Realizability checks

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
- `numerically_indeterminate`: the classification depends on an eigenvalue
  within the tolerance or on digits the input moments do not carry.

`RealizabilityResult::realizable()` is true for `strictly_realizable` and
`boundary`. `RealizabilityResult::on_boundary()` is true only for `boundary`.

Diagnostics include the highest checked order, the limiting or failed matrix
name and order, a minimum pivot/eigenvalue-like diagnostic, the tolerance used,
and a human-readable message; with relative tolerance, the pivot and the
tolerance refer to the standardized matrices. Callers should treat the status
as authoritative and use the diagnostics for reporting, debugging, or adaptive
failure handling. Diagnostics are not a request to continue with invalid
moments.

Invalid input means the request itself is malformed or unsupported, such as an
empty moment sequence, non-finite moment values, non-positive or non-finite
tolerances, a bounded-interval request without explicit bounds, or invalid
bounded-interval limits.

`not_realizable` means the supplied numbers form a well-defined request but
cannot represent a valid nonnegative measure on the selected support. Examples
include negative `m_0`, zero total number density with nonzero higher moments,
or moment sequences that fail the support-specific realizability conditions.
The all-zero sequence is treated as a boundary-realizable zero measure.

`RealizabilityOptions` controls numerical tolerance; the next section
describes how it enters the classification.

## Numerical classification

### Standardized coordinate

Realizability on a support is invariant under an affine change of coordinate
$y = (x - c)/s$ with $s > 0$, which maps the support to
$\left[(a - c)/s, (b - c)/s\right]$. The moment and localizing matrices of
$y$ are congruent to those of $x$, so by Sylvester's law of inertia the signs
of their eigenvalues, and hence the classification, are unchanged. With
`relative_tolerance = true`, the check therefore tests the standardized
moments

$$\nu_k = \frac{1}{2^{e_0}\, s^k} \sum_{j=0}^{k} \binom{k}{j} (-c)^{k-j} m_j,
\qquad c = \frac{m_1}{m_0},$$

where $2^{e_0}$ is the power of two of $m_0$ and $s$ is a power of two near the
spread: $\left(|\mu_k|/m_0\right)^{1/k}$ for the first even central moment
$\mu_k$ that is resolved, otherwise $|c|$, or one if $c = 0$, so that a point
mass is not stretched. The powers of two make the scalings exact. The relative
tolerance then applies to the standardized matrices, and the classification
depends neither on the location of the distribution, nor on the unit of its
coordinate, nor on $m_0$. With `relative_tolerance = false`,
the raw moments are tested and `tolerance` is an absolute eigenvalue
threshold; this mode is not invariant.

### Input precision

Each input moment is taken as exact to its last digit, a relative error of
$\varepsilon$. Evaluating the centring sum in double precision adds its own
rounding. With round to nearest, every operation returns the exact result times
$1 + \delta$ with $|\delta| \le \mathbf{u} = \varepsilon/2$, so a value that
passes through $p$ roundings is off by a factor that differs from one by at
most $(1 + \mathbf{u})^p - 1 \le p\mathbf{u}/(1 - p\mathbf{u})$. The power
$(-c)^i$ takes $i - 1$ roundings for $i \ge 1$ and none for $i = 0$, the term
$\binom{k}{j}(-c)^{k-j} m_j$ two more (the binomial coefficients are exact
integers), and the running sum at most $k$ more, so no term passes through
more than $2k + 1$ roundings. The evaluated sum therefore differs from the
exact one by at most $(2k+1)\mathbf{u}/(1 - (2k+1)\mathbf{u})$, less than
$(k + 1)\varepsilon$, times the sum of the absolute terms; a fused
multiply-add only removes roundings. Together, the uncertainty of the standardized moment $\nu_k$ is

$$u_k = \frac{(k + 2)\,\varepsilon}{2^{e_0}\, s^k} \sum_{j=0}^{k} \binom{k}{j}
|c|^{k-j} |m_j|.$$

The Frobenius norm $\delta$ of the Hankel matrix of these uncertainties bounds
how far they can move any eigenvalue of a tested matrix. Localizing sequences
carry the uncertainties of the moments they combine. The rounding of building
them, and of the transformed support bounds, is relative to the terms involved
and lies within the roundoff allowance of the eigen solve. For a distribution
far from the origin, $u_k$ grows like $\varepsilon (|c|/\sigma)^k$ relative to
the standardized moments, which sets the resolution limit stated in the user
documentation.

### Eigenvalue partition

Every eigenvalue $\lambda$ of a tested matrix $A$ is classified with the
threshold $\tau = \tau_{\mathrm{user}} + \delta$ and the zero band
$\zeta = 1000\,\varepsilon \lVert A \rVert + \delta$:

| Class | Condition |
|---|---|
| negative | $\lambda < -\tau$ |
| zero | $|\lambda| \le \zeta$ |
| ambiguous | otherwise, $\lambda < \tau$ |
| positive | $\lambda \ge \tau$ |

A matrix with a negative eigenvalue is not positive semidefinite; otherwise
one with an ambiguous eigenvalue is `numerically_indeterminate`; otherwise one
with a zero eigenvalue is singular, which makes the sequence `boundary`;
otherwise it is positive definite. A matrix whose uncertainty $\delta$ reaches
the information scale of its sequence carries no reliable digit and is
`numerically_indeterminate`. That scale is the standardized mass $\nu_0$ for
the moments, and $\nu_0$ times the size of the coefficients of the localizing
polynomial for a localizing sequence: $(1 + |a|)\,\nu_0$ for $x - a$ and
$b - x$, and $(1 + |a + b| + |ab|)\,\nu_0$ for $(x - a)(b - x)$, with the
bounds in the standardized coordinate.

The same partition defines the numerical rank and the kernel used by the
singular-matrix tests. A kernel projection of the next column is compared
against a threshold that includes the uncertainty of the column and the
rotation of the kernel vectors, which by the Davis-Kahan theorem is at most
$\delta$ over the gap between the zero and the positive eigenvalues, bounded
below by the smallest positive eigenvalue minus $\zeta$.
The flat-extension test compares the numerical ranks of $H_n$ and its
principal submatrix $H_{n-1}$, whose eigenvalues interlace those of $H_n$. In
exact arithmetic $H_n$ is flat when the ranks agree and has one more otherwise.
The zero band of $H_n$ is the wider of the two, so its rank can also come out
lower than that of $H_{n-1}$: it then resolves the same eigenvalues less
finely, and the extension counts as flat within its precision. A higher rank
of $H_n$ rejects the sequence only if every zero of $H_{n-1}$ lies within the
eigen-solve roundoff, so that $H_{n-1}$ is singular for the data as given. If a
zero of $H_{n-1}$ holds only within the data uncertainty, $H_{n-1}$ may as well
be positive definite, which makes the sequence realizable whatever $H_n$, and
the result is `numerically_indeterminate`.

### Evidence

`shift_and_scale_preserve_classification` in
`tests/verification/realizability/realizability_analytical.cpp` checks
strictly realizable, boundary and non-realizable sequences on every support for
spreads from $10^{-9}$ to $10^9$ and locations up to the shift each resolves,
and `resolution_limit_never_reverses_classification` checks that beyond it, up
to $10^8$ spreads, a classification only gives way to `boundary` or
`numerically_indeterminate`. `exact_point_sequences_are_never_rejected` checks
50,000 sequences of up to four weighted points whose moments are exact in
double precision, many of them tens of spreads from the origin: none may be
rejected, and none on the boundary may be called strictly realizable. The unit
tests in `tests/unit/realizability.cpp` cover the partition of every
eigenvalue, the flat-extension test within the precision and within the data
uncertainty, and the loss of every reliable digit.

## Runtime checks in phase-space evolution

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

## Interaction with QMOM and CQMOM

Realizability is a constraint on moment representations. It is not a
replacement for QMOM or CQMOM reconstruction status.

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

Each layer reports its own status. A successful earlier layer does not
guarantee success in a later one.

## Failure handling

Realizability APIs use explicit result objects rather than exceptions for
normal classification failure. Callers should check
`RealizabilityResult::status` or `realizable()` before using a moment sequence
as accepted input under the selected support settings.

Failure handling rules are:

- `invalid_input` means the check request is malformed or unsupported;
- `not_realizable` means the moments are incompatible with the selected support;
- `numerically_indeterminate` means the classification depends on an eigenvalue
  within the tolerance or on digits the input does not carry, and the sequence
  should not be accepted by runtime checks;
- `boundary` is mathematically accepted but may require an explicit acceptance
  setting;
- `strictly_realizable` is accepted by the check.

QMOMKit does not automatically correct failed sequences. It does not alter the
moments, adjust supports, relax tolerances, or reduce quadrature orders as part
of explicit realizability checking. Retry behavior applies only to
moment-source integration with runtime checking enabled and configured retry
controls.

## Public API boundary

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

## Extension principles

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

The goal is conservative numerical validation: realizability support should
help developers and drivers detect invalid moment states early without hiding
the source of the failure.