# Realizability

A raw moment sequence $m_0, \ldots, m_N$ is realizable on a support when a
nonnegative measure on that support has these moments.
`check_realizability(...)` in [`src/realizability.cpp`](../../src/realizability.cpp)
decides this for one coordinate; a realizability policy applies the check
inside a workflow. A successful reconstruction is no realizability check
([QMOM](3_qmom.md)). Usage is in the
[user documentation](../user_new/4_moments_and_quadrature.md#check-realizability).

## Moment matrices

The check tests standardized moments $\nu_k$ (next section), with the support
bounds $a$, $b$ mapped to the same coordinate. For a sequence $s$, $H_n(s)$ is
the $(n+1) \times (n+1)$ Hankel matrix with entries $s_{i+j}$,
$0 \le i, j \le n$, and $H_n = H_n(\nu)$. Localizing sequences multiply the
measure by $x - a$, $b - x$ and $(x - a)(b - x)$:

```math
\ell^{a}_k = \nu_{k+1} - a\,\nu_k, \qquad
\ell^{b}_k = b\,\nu_k - \nu_{k+1}, \qquad
\ell^{ab}_k = (a + b)\,\nu_{k+1} - ab\,\nu_k - \nu_{k+2}.
```

For a realizable sequence, every tested matrix is positive semidefinite. With
$N = 2n$ or $N = 2n + 1$, the check tests, in this order:

| Support | $N = 2n$ | $N = 2n + 1$ |
| --- | --- | --- |
| `real_line` | $H_n$; if singular, flat extension over $H_{n-1}$ | $H_n$; if singular, range test |
| `nonnegative`, $a$ the image of 0 | $H_n$, $S_{n-1} = H_{n-1}(\ell^a)$; if $S_{n-1}$ is singular, range test | $H_n$, $S_n = H_n(\ell^a)$; if $H_n$ is singular, range test |
| `bounded_interval` | $H_n$, $H_{n-1}(\ell^{ab})$ | $H_n(\ell^a)$, $H_n(\ell^b)$ |

For $N = 0$ only $H_0$ is tested. `RealizabilityResult::failed_matrix` names
them `H_n`, `S_n`, `L_interval_n`, `L_lower_n` and `L_upper_n`. On a bounded
interval positive semidefiniteness is also sufficient. On the real line and the
half-line a singular matrix needs a further test: its rank must not grow with
the next moments
[[Curto & Fialkow, 1991](https://www.math.uh.edu/~hjm/v017n4/0603CURTO.pdf)].

## Standardization

`standardized_moments` shifts and scales the coordinate:

```math
\nu_k = 2^{-e_0 - k e_s} \sum_{j=0}^{k} \binom{k}{j} (-c)^{k-j} m_j,
\qquad c = \frac{m_1}{m_0}.
```

$2^{e_0}$ is the smallest power of two above $m_0$, and $2^{e_s}$ the smallest
power of two above the spread: $(|\mu_k| / m_0)^{1/k}$ for the first even
central moment $\mu_k$ that exceeds twice its uncertainty, otherwise $|c|$, or
one if $c = 0$, so that a point mass is not stretched. A bound $a$ becomes
$(a - c)/2^{e_s}$.

Realizability is invariant under an affine change of coordinate, and the
moment and localizing matrices change by a congruence, which by Sylvester's
law of inertia [[Sylvester, 1852](https://www.doi.org/10.1080/14786445208647087)] keeps the signs of their eigenvalues. The tolerance applies to
the standardized matrices, so the classification depends neither on the
location of the distribution, nor on the unit of its coordinate, nor on $m_0$.
Scaling by powers of two is exact, so only the centring sum rounds.

## Tolerance

The check has one relative tolerance $t$, `RealizabilityOptions::tolerance`,
default `1e-12`, defined in
[`include/qmomkit/realizability.hpp`](../../include/qmomkit/realizability.hpp).
It must be finite and positive. For a tested matrix $A$ it becomes
$t \lVert A \rVert_\infty$.

## Classification

### Input uncertainty

Each input moment is taken as exact to its last digit, a relative error
$\varepsilon$. Each term of the centring sum passes through at most $2k + 1$
roundings of at most $\varepsilon/2$, less than $(k + 1)\varepsilon$ of the sum
of the absolute terms. With the input error, the standardized moment $\nu_k$
carries the uncertainty

```math
u_k = (k + 2)\,\varepsilon\, 2^{-e_0 - k e_s} \sum_{j=0}^{k} \binom{k}{j} |c|^{k-j} |m_j|.
```

A localizing sequence combines the uncertainties with the absolute
coefficients of its polynomial. The Frobenius norm $\delta$ of the Hankel
matrix of the uncertainties bounds how far they can move any eigenvalue of
$A$. The rounding of the localizing sequences and of the mapped bounds is
relative to the terms involved and lies within the roundoff allowance below.

### Eigenvalue partition

`analyze_matrix` computes every eigenvalue with the Jacobi eigensolver
([Linear algebra](5_linear_algebra.md)), and `classify_eigenvalue` assigns
each one a class, tested in this order, with the zero band $\zeta$ and the
threshold $\tau$:

```math
\zeta = 1000\,\varepsilon \lVert A \rVert_\infty + \delta, \qquad
\tau = t \lVert A \rVert_\infty + \delta.
```

| Class | Condition |
| --- | --- |
| zero | $\lvert\lambda\rvert \le \zeta$ |
| negative | $\lambda < -\tau$ |
| ambiguous | $\lambda < \tau$ |
| positive | otherwise |

The factor 1000, `roundoff_allowance`, bounds the rounding of the eigensolve.
An eigenvalue inside that band has no sign, so the zero class is tested first,
and a tolerance below $1000\varepsilon$ resolves nothing finer.

A matrix with a negative eigenvalue is not positive semidefinite; otherwise,
one with an ambiguous eigenvalue is indeterminate, one with a zero eigenvalue
singular, and any other positive definite. A matrix is also
indeterminate if the eigen solve does not converge, or if $\delta$ reaches the
information scale of its sequence, so that no digit is reliable. That scale is
$\nu_0$ for the moments, and $\nu_0$ times the absolute coefficients of the
localizing polynomial otherwise: $(1 + |a|)\,\nu_0$, $(1 + |b|)\,\nu_0$ or
$(1 + |a + b| + |ab|)\,\nu_0$.

The first matrix that is not positive semidefinite makes the sequence
`not_realizable`, the first indeterminate one `numerically_indeterminate`. A
singular matrix makes it `boundary` unless its singular-matrix test decides
otherwise.

### Singular matrices

The same partition defines the numerical rank, the number of positive
eigenvalues, and the kernel, spanned by the eigenvectors of zero eigenvalues.

The **range test** (`check_range_consistency`) requires the next column of a
singular Hankel matrix to lie in its range. The largest projection of that
column on a kernel vector is compared with a threshold $\theta$: above
$\theta$ the sequence is `not_realizable`, above $\theta/10$
`numerically_indeterminate`, and otherwise `boundary`. $\theta$
(`kernel_projection_tolerance`) is the largest of $\tau$, $t$ times the
largest column entry, and ten times the uncertainty of the projection: that
of the column plus its norm times the rotation of the kernel vectors. By the
Davis–Kahan theorem [[Davis & Kahan, 1970](https://www.doi.org/10.1137/0707001)] that rotation is at most $\delta$ over the gap between the
zero and the positive eigenvalues, which the smallest positive eigenvalue
minus $\zeta$ bounds from below.

The **flat-extension test** decides a singular $H_n$ on the real line: it is
flat if it has no more positive eigenvalues than $H_{n-1}$. $H_{n-1}$ has a
narrower zero band, so $H_n$ may count fewer positive eigenvalues; that is a
coarser view of the same eigenvalues and still counts as flat. A higher
rank makes the sequence `not_realizable` only if every zero eigenvalue of
$H_{n-1}$ lies within its roundoff allowance,
$1000\,\varepsilon \lVert H_{n-1} \rVert_\infty$. If one is zero only within
the data uncertainty, $H_{n-1}$ may as well be positive definite, which makes
the sequence realizable whatever $H_n$, and the result is
`numerically_indeterminate`.

## Statuses and result

Before any matrix is built, `invalid_input` is returned for an empty sequence,
a non-finite moment, a tolerance that is not finite and positive,
`MomentSupport::bounded_interval` without a `BoundedInterval`, or an interval
whose bounds are not finite with `lower < upper`. A negative $m_0$, or a zero
$m_0$ with a nonzero higher moment, is `not_realizable`; the all-zero
sequence, the zero measure, is `boundary`. `realizable()` is true for `strictly_realizable` and
`boundary`, `on_boundary()` only for `boundary`.

`checked_order` is $N$. `failed_matrix` and `failed_matrix_order` name the
matrix that decided a failure or, on success, the matrix closest to failing:
the one with the smallest ratio of its smallest eigenvalue to its threshold.
`minimum_pivot` is that eigenvalue, or the kernel projection of a range test,
and `tolerance_used` the threshold it was compared with, both for the
standardized matrices.

## Precision limit

For a distribution centred at $c$ with width $\sigma$, $u_k$ grows relative to
$\nu_k$ like $\varepsilon (|c|/\sigma)^k$: the centring sum cancels the
leading digits. Where it cancels all of them, the check reports `boundary` or
`numerically_indeterminate` instead of the exact verdict, but never the
opposite one: a realizable sequence is never rejected, and a non-realizable or
boundary one is never strictly realizable.

## Policies

`RealizabilityPolicy`, in
[`include/qmomkit/realizability_policy.hpp`](../../include/qmomkit/realizability_policy.hpp),
holds `check_mode` (`None`, the default, `FinalOnly` or `EveryStage`),
`support` (default `real_line`), `interval` (required for `bounded_interval`),
`options` (its tolerance), and `accept_boundary` (default `true`). A policy accepts
`strictly_realizable`, and `boundary` if `accept_boundary` is set; it rejects
the other statuses. `internal::realizability_policy_error`, declared in
[`src/internal/realizability.hpp`](../../src/internal/realizability.hpp),
validates a policy with the rules `check_realizability` applies to its options
and interval; a `None` policy is always valid.

| Where | Policy | Checked moments |
| --- | --- | --- |
| `QmomCase::state_candidate(moments, policy)` | one | the candidate moments |
| `CqmomPreparedCase::state_candidate(moments, policies)` | one optional per coordinate | each coordinate's marginal moments |
| QMOM moment-source step | `QmomPhaseSpaceEvolutionOptions::realizability_policy` | the source-step candidates |
| CQMOM moment-source step, `advance_cqmom_moments_with_source` | `MomentSourceUpdateOptions::coordinate_realizability_policies` | each coordinate's marginal moments of the candidates |

`state_candidate` checks once, before the reconstruction, for any mode other
than `None`. A rejection returns `realizability_rejected`, with
`rejected_coordinate` for CQMOM, and an invalid policy `invalid_input`.
Boundary moments with $m_0$ below `QmomOptions::cutoff` are accepted even if
`accept_boundary` is false, so that they are reported `empty`, as without a
check. A QMOM sequence too short for the node count is left to the
reconstruction.

In the moment-source step a policy is an acceptance predicate; how it
interacts with the step is in
[Phase-space evolution](8_phase_space_evolution.md#moment-source). The
stateless `advance_qmom_moments_with_source` rejects coordinate policies; a
caller who wants a check installs `make_realizability_acceptance_predicate(...)` or
`make_realizability_stage_acceptance_predicate(...)`.

In CQMOM a policy checks only the marginal moments $M_{0,\ldots,k,\ldots,0}$
of its coordinate, a necessary condition; mixed and conditional moments are
never checked. A policy indexed beyond the phase-space dimensions is invalid
whatever its mode.

## Tests

[`tests/unit/realizability.cpp`](../../tests/unit/realizability.cpp) covers
the classification and
[`tests/unit/realizability_policy.cpp`](../../tests/unit/realizability_policy.cpp)
the policy predicates.
[`tests/verification/realizability/`](../../tests/verification/realizability/)
checks analytical sequences on every support, invariance under shift and
scale, the precision limit, and sequences of exact weighted points.
