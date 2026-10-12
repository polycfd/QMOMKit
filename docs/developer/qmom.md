# QMOM

QMOM turns the raw moments $m_0,\ldots,m_{2N-1}$ of a one-dimensional number
density into an $N$-node Gauss quadrature with the Wheeler algorithm
[[Wheeler, 1974](https://www.doi.org/10.1216/RMJ-1974-4-2-287)]. The same
inversion reconstructs every one-dimensional sequence of recursive
[CQMOM](cqmom.md). Public usage is in
[Moments and quadrature](../user/4_moments_and_quadrature.md#reconstruct);
moment and ordering conventions are in
[Numerical conventions](numerical_conventions.md).

## Wheeler inversion

`compute_recurrence` in [`src/wheeler.cpp`](../../src/wheeler.cpp) fills
Wheeler's $\sigma$ table from $\sigma_{0,l}=0$, $\sigma_{1,l}=m_{l-1}$,
$a_0=m_1/m_0$, $b_0=0$, for $k=2,\ldots,N$:

```math
\sigma_{k,l}=\sigma_{k-1,l+1}-a_{k-2}\,\sigma_{k-1,l}-b_{k-2}\,\sigma_{k-2,l},
\qquad
a_{k-1}=\frac{\sigma_{k,k+1}}{\sigma_{k,k}}-\frac{\sigma_{k-1,k}}{\sigma_{k-1,k-1}},
\qquad
b_{k-1}=\frac{\sigma_{k,k}}{\sigma_{k-1,k-1}} .
```

$a_k$ and $b_k$ are the three-term recurrence coefficients of the polynomials
orthogonal with respect to the moments. The Jacobi matrix $J$ has diagonal
$a_0,\ldots,a_{N-1}$ and off-diagonal $\sqrt{b_1},\ldots,\sqrt{b_{N-1}}$. With
normalized eigenvectors $Jv_i=\lambda_i v_i$,

```math
\xi_i=\lambda_i,\qquad w_i=m_0\,v_{0,i}^2 .
```

`symmetric_jacobi_eigen` ([Linear algebra](linear_algebra.md)) returns the
eigenvalues in ascending order, so abscissas come out sorted with their
weights. The moments are used as given, without centering or normalization,
and $\sum_i w_i=m_0$.

Since $\sigma_{k,k}=m_0\prod_{i=1}^{k-1}b_i$, a zero $\sigma_{k,k}$ means the
moments are carried by at most $k-1$ points, and a negative one that they are
not realizable with $k$ nodes. Both tests below act on $\sigma_{k,k}$.

## Input checks and the cutoff

`compute_quadrature_at_scale` checks, in order:

1. `rmin` and `eabs` finite and non-negative, `cutoff` finite and positive;
   $N>0$; at least $2N$ moments, the first $2N$ finite. Failure is
   `invalid_input` with empty vectors. Higher moments are not read.
2. $m_0<0$: `not_realizable`.
3. $m_0<$ `cutoff`: `not_realizable`, message `Moment sequence has m0 < cutoff.`
4. $N=1$: one node $\xi=m_1/m_0$, $w=m_0$.

`cutoff` (default `1e-30`) is an absolute threshold on $m_0$ and the only
decision of the inversion that depends on the number density. Every other test
gives the same result after scaling all moments by a constant or changing the
unit of the coordinate. State candidates report $0\le m_0<$ `cutoff` as `empty`
instead of a failure ([Realizability](realizability.md#policies)).

## Rounding bound on the recurrence

`sigma_tolerances` reruns the $\sigma$ recurrence with every term in absolute
value and takes a multiple of its diagonal as the rounding of $\sigma_{k,k}$:

```math
\bar\sigma_{1,l}=|m_{l-1}|,\qquad
\bar\sigma_{k,l}=\bar\sigma_{k-1,l+1}+|a_{k-2}|\,\bar\sigma_{k-1,l}+|b_{k-2}|\,\bar\sigma_{k-2,l},
\qquad
\tau_k=100\,\varepsilon\,\max\!\left(\bar\sigma_{k,k},\ |m_0|\,s^{2(k-1)}\right),
```

with $\varepsilon$ the double epsilon and $s$ the coordinate scale below
($s=0$ drops the floor for $k\ge2$). Both $\sigma_{k,k}$ and $\tau_k$ scale
as $m_0$ times the coordinate to the power $2(k-1)$, so the test depends on
neither unit nor number density. $\tau_k$ also grows with the cancellation of
a distribution far from the origin.

## Fixed node count

With `adaptive = false`, the first $k\in\{2,\ldots,N\}$ with
$|\sigma_{k,k}|\le\tau_k$ ends the inversion with `not_realizable` and the
message `Moment sequence supports at most k-1 nodes: …`, with $k-1$ filled in.
Without this test, rounding would decide whether a sequence on the boundary
passes, and that changes with unit and number density; an accepted rule would
also carry nodes made only of rounding. A non-finite recurrence or a
negative $b_k$ is also `not_realizable`. Otherwise the result has exactly $N$
nodes; `rmin` and `eabs` are not applied.

## Adaptive reduction

Adaptive inversion [[Marchisio & Fox, 2013](https://www.doi.org/10.1017/CBO9781139016599)]
(`adaptive = true`, the default) tries candidate counts
$n=N,N-1,\ldots,2$:

1. Scan $k=n,\ldots,2$. $\sigma_{k,k}<-\tau_k$ rejects the candidate and the
   next count is tried. $|\sigma_{k,k}|\le\tau_k$ lowers the usable count to
   $k-1$; the scan goes on, so the lowest such $k$ decides.
2. For $n_1$ from the usable count down to 2, form the $n_1$-node rule from the
   leading coefficients and return the first that passes the two criteria
   below.

If the usable count reaches 1, no rule passes, or every candidate count meets a
negative $\sigma_{k,k}$, the result is the one-node rule $(m_1/m_0,\ m_0)$
with `success`. A sequence with negative variance therefore succeeds with one
node that matches $m_0$ and $m_1$ only (`tests/unit/wheeler.cpp`).

A rule with $n_1$ nodes is accepted when

```math
\frac{\min_i w_i}{\max_i w_i}>r_{\min},
\qquad
\begin{cases}
\dfrac{\min_{i\ne j}|\xi_i-\xi_j|}{\max_{i,j}|\xi_i-\xi_j|}>e_{\mathrm{abs}}, & n_1>2,\\
\dfrac{|\xi_1-\xi_0|}{\max(|\xi_0|,|\xi_1|)}>e_{\mathrm{abs}}, & n_1=2 .
\end{cases}
```

Two nodes have a single distance, whose ratio to the span is always one, so
they use the larger abscissa magnitude instead; two nodes at zero fail. Both
criteria are independent of the unit, but the two-node test depends on the
origin: a pair whose gap is below `eabs` times its distance from zero collapses
to one node. `rmin` and `eabs` default to `1e-10` in
`QmomOptions` ([`include/qmomkit/qmom.hpp`](../../include/qmomkit/qmom.hpp))
and `WheelerOptions`
([`include/qmomkit/wheeler.hpp`](../../include/qmomkit/wheeler.hpp)); only
`adaptive` is a `.qmom` key
([Reference](../user/7_reference.md#reconstruction-options)).

## Coordinate scales

`compute_quadrature(...)` and `compute_qmom(...)` use $s=0$: given moments
carry no solve rounding. CQMOM conditional moments do, so their inversion,
`internal::compute_conditional_quadrature`, receives a positive $s$ in
adaptive mode ([CQMOM](cqmom.md#the-recursion)): `cqmom_coordinate_scales`
([`src/coordinate_scales.cpp`](../../src/coordinate_scales.cpp)) takes
$s=\sqrt{M_2/M_0}$ of the marginal of the coordinate being inverted, else
$|M_1|/M_0$, else one. The CQMOM root inverts given moments and uses $s=0$. In
fixed mode CQMOM passes $s=0$ at every level, so its fixed-count test uses the
bound without the floor.

## Abscissa bounds

`QmomCase::abscissa_bounds` and the `compute_qmom(...)` overload taking
`std::optional<AbscissaBounds>` route through `select_bounded_quadrature`
([`src/internal/bounded_quadrature.hpp`](../../src/internal/bounded_quadrature.hpp)).
The bounds are an acceptance test around ordinary Wheeler inversion. A
candidate count $c$ inverts the first $2c$ moments and is accepted when it
succeeds with exactly $c$ nodes, all in $[\ell-t,\ u+t]$ with

```math
t=10^{-12}\max(|\ell|,|u|).
```

Rounding can move a node that lies on a bound slightly outside it. That error
scales with the coordinate, so the margin is relative to the bounds and decides
the same in any unit. A fixed reconstruction
evaluates $c=N$ only; an adaptive one tries $c=N,\ldots,1$ and keeps the first
accepted, and if $c=1$ also fails, returns its status and message. Bounds must
be finite with $\ell<u$, and $2N$ finite moments are required whatever count
is accepted.

## Cases, results and statuses

`QmomCase` holds `moments`, `node_count`, `qmom_options` and
`abscissa_bounds`; `compute_qmom(case)` forwards all four. `compute_qmom(...)`
maps `QuadratureStatus` one to one onto `QmomStatus`:

| Status | Cause |
| --- | --- |
| `success` | A rule was accepted, in adaptive mode possibly with fewer nodes than requested. |
| `invalid_input` | Invalid options or bounds, no nodes, too few or non-finite moments. |
| `not_realizable` | $m_0<0$, $m_0<$ `cutoff`, a fixed count the moments do not support, negative $b_k$, non-finite recurrence or eigenpairs, an unconverged eigensolve, or nodes outside the bounds. |

`ok()` is true only for `success`; the accepted node count is `weights.size()`.
Data returned with a failure is diagnostic: the one-node rule $(m_1/m_0,\ m_0)$
where finite, $(m_1/m_0,\ |m_0|)$ for negative $m_0$, the placeholder node 0
with weight 1 below the cutoff, the rejected rule for a bounds failure (the full
rule in fixed mode, the one-node candidate in adaptive mode), and nothing for
invalid input. Callers test the status, never whether the vectors are empty.

## Where the code lives

| Function | File | Role |
| --- | --- | --- |
| `compute_qmom` | [`src/public_api.cpp`](../../src/public_api.cpp) | Public entry; bounds dispatch, status mapping. |
| `compute_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Public low-level Wheeler inversion, $s=0$. |
| `compute_conditional_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Inversion of CQMOM conditional moments with scale $s$. |
| `select_bounded_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Bounds acceptance for QMOM and CQMOM. |
| `cqmom_coordinate_scales` | [`src/coordinate_scales.cpp`](../../src/coordinate_scales.cpp) | Coordinate scales $s$. |

Coverage is in `tests/unit/wheeler.cpp`, `tests/unit/abscissa_bounds.cpp` and
`tests/verification/qmom/`.
