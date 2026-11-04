# QMOM

QMOM turns the raw moments $m_0,\ldots,m_{2N-1}$ of a one-dimensional number
density into an $N$-node Gauss quadrature with the Wheeler algorithm
[[Wheeler, 1974](https://www.doi.org/10.1216/RMJ-1974-4-2-287)]. The same
inversion reconstructs every one-dimensional sequence of recursive
[CQMOM](cqmom.md). The hyperbolic family reads $m_{2N}$ as well and returns up
to $N+1$ nodes ([Hyperbolic family](#hyperbolic-family)). Public usage is in
[Moments and quadrature](../user/4_moments_and_quadrature.md#reconstruct);
moment and ordering conventions are in
[Numerical conventions](numerical_conventions.md).

## Wheeler inversion

`compute_recurrence` in [`src/wheeler.cpp`](../../src/wheeler.cpp) computes
Wheeler's $\sigma$ table row by row from $\sigma_{-1,l}=0$, $\sigma_{0,l}=m_l$,
$a_0=m_1/m_0$, $b_0=0$, for $k=1,\ldots,N-1$:

```math
\sigma_{k,l}=\sigma_{k-1,l+1}-a_{k-1}\,\sigma_{k-1,l}-b_{k-1}\,\sigma_{k-2,l},
\qquad
a_k=\frac{\sigma_{k,k+1}}{\sigma_{k,k}}-\frac{\sigma_{k-1,k}}{\sigma_{k-1,k-1}},
\qquad
b_k=\frac{\sigma_{k,k}}{\sigma_{k-1,k-1}} .
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

Since $\sigma_{k,k}=m_0\prod_{i=1}^{k}b_i$, a zero $\sigma_{k,k}$ means the
moments are carried by at most $k$ points, and a negative one after positive
ones that $m_0,\ldots,m_{2k}$ are the moments of no non-negative density.
Both tests below act on $\sigma_{k,k}$.

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

`compute_recurrence` also runs the $\sigma$ recurrence with every term in
absolute value and takes a multiple of its diagonal as the rounding of
$\sigma_{k,k}$:

```math
\bar\sigma_{0,l}=|m_l|,\qquad
\bar\sigma_{k,l}=\bar\sigma_{k-1,l+1}+|a_{k-1}|\,\bar\sigma_{k-1,l}+|b_{k-1}|\,\bar\sigma_{k-2,l},
\qquad
\tau_k=100\,\varepsilon\,\max\!\left(\bar\sigma_{k,k},\ |m_0|\,s^{2k}\right),
```

with $\varepsilon$ the double epsilon and $s$ the coordinate scale below
($s=0$ drops the floor). Both $\sigma_{k,k}$ and $\tau_k$ scale as $m_0$
times the coordinate to the power $2k$, so the test depends on neither unit
nor number density. $\tau_k$ also grows with the cancellation of
a distribution far from the origin.

`compute_wheeler_recurrence` returns $a_k$, $b_k$, $\sigma_{k,k}$ and $\tau_k$
for $k=0,\ldots,N-1$ as `compute_recurrence` computes them for the
reconstruction, with $s=0$; entry 0 holds $\sigma_{0,0}=m_0$ and $\tau_0=0$,
since no test is made at $k=0$. The test of
[Fixed node count](#fixed-node-count) and step 1 of
[Adaptive reduction](#adaptive-reduction) act on these values. Its input
checks are those of `compute_quadrature_at_scale`, without the option and
`cutoff` tests, and $m_0=0$ is `not_realizable`.

## Fixed node count

With `adaptive = false`, the first $k\in\{1,\ldots,N-1\}$ with
$|\sigma_{k,k}|\le\tau_k$ ends the inversion with `not_realizable` and the
message `Moment sequence supports at most k nodes: …`, with $k$ filled in.
Without this test, rounding would decide whether a sequence on the boundary
passes, and that changes with unit and number density; an accepted rule would
also carry nodes made only of rounding. A non-finite recurrence or a
negative $b_k$ is also `not_realizable`. Otherwise the result has exactly $N$
nodes; `rmin` and `eabs` are not applied.

## Adaptive reduction

Adaptive inversion [[Marchisio & Fox, 2013](https://www.doi.org/10.1017/CBO9781139016599)]
(`adaptive = true`, the default) selects its rule in two steps:

1. The usable count is $n=k$ for the first $k\in\{1,\ldots,N-1\}$ with
   $\sigma_{k,k}\le\tau_k$, zero within its rounding or negative, and $n=N$
   if there is none.
2. For $n_1=n,n-1,\ldots,2$, form the $n_1$-node rule from the leading
   coefficients and return the first that passes the two criteria below.

If $n=1$ or no rule passes, the result is the one-node rule
$(m_1/m_0,\ m_0)$ with `success`. A sequence with negative variance therefore
succeeds with one node that matches $m_0$ and $m_1$ only.

A rule with $n_1$ nodes is accepted when

```math
\frac{\min_i w_i}{\max_i w_i}>r_{\min},
\qquad
\begin{cases}
\dfrac{\min_{i\ne j}|\xi_i-\xi_j|}{\max_{i,j}|\xi_i-\xi_j|}>e_{\mathrm{abs}}, & n_1>2,\\
\dfrac{|\xi_1-\xi_0|}{\max(|\xi_0|,|\xi_1|)}>e_{\mathrm{abs}}, & n_1=2 .
\end{cases}
```

The weight ratio and, for $n_1>2$, the separation ratio are the criteria of
[[Yuan & Fox, 2011](https://www.doi.org/10.1016/j.jcp.2011.07.020)]. The
abscissas are sorted, so the smallest distance is the smallest gap between
neighbors and the largest is the span $\xi_{n_1-1}-\xi_0$. Two nodes have a
single distance, whose ratio to the span is always one. Yuan & Fox compare that
distance itself with `eabs`, which depends on the unit; QMOMKit divides it by
the larger abscissa magnitude instead, and two nodes at zero fail. Both
criteria are then independent of the unit, but the two-node test depends on the
origin: a pair whose gap is below `eabs` times its distance from zero collapses
to one node. `rmin` and `eabs` default to `1e-10` in
`QmomOptions` ([`include/qmomkit/qmom.hpp`](../../include/qmomkit/qmom.hpp))
and `WheelerOptions`
([`include/qmomkit/wheeler.hpp`](../../include/qmomkit/wheeler.hpp)); neither
has a `.qmom` key
([Reference](../user/8_reference.md#reconstruction-options)).

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
the same in any unit. A fixed reconstruction evaluates $c=N$ only; an adaptive
one tries $c=N,\ldots,1$ and keeps the first accepted, and if $c=1$ also fails,
returns its status and message. In evolution, an adaptive reconstruction starts
from the state's active count instead of $N$, and the post-source one may then
try one node more ([Phase-space evolution](phase_space_evolution.md#node-counts)).
Bounds must be finite with $\ell<u$, and $2N$ finite moments are required
whatever count is accepted.

## Hyperbolic family

HyQMOM [[Fox & Laurent, 2022](https://www.doi.org/10.1137/21M1406143)] closes
$m_{2n+1}$ from $m_0,\ldots,m_{2n}$ on the real line so that the moment
system of free transport is globally hyperbolic, as shown for $n\le9$ (Fox &
Laurent, 2022, Theorems 4.6 and 4.8). The internal
`compute_hyperbolic_quadrature`
([`src/internal/hyperbolic_quadrature.hpp`](../../src/internal/hyperbolic_quadrature.hpp))
returns the $(n+1)$-node rule of that closure, which reproduces
$m_0,\ldots,m_{2n}$ and the closed $m_{2n+1}$. It uses $s=0$, as
`compute_quadrature(...)` does.

In `.qmom` files, `[phase_space] qmom_family` sets the family
([Configuration and output](configuration_and_output.md#mapping-to-objects)).
`compute_qmom(...)` dispatches on `QmomOptions::family` in
`compute_qmom_from_moments`
([`src/public_api.cpp`](../../src/public_api.cpp)): `classical` takes the
paths above, and `hyperbolic` calls `compute_hyperbolic_quadrature` with
`node_count` as $n$. HyQMOM is defined on the real line, so abscissa bounds
with `hyperbolic` are `invalid_input`, with the message
`The hyperbolic QMOM family takes no abscissa bounds.`, before any
reconstruction. CQMOM refuses the family ([CQMOM](cqmom.md)), and states
and their evolution follow the count rule of
[Phase-space evolution](phase_space_evolution.md#node-counts). The moments a
family reads, $2n$ or $2n+1$, have one definition,
`internal::qmom_required_moment_count`
([`src/internal/qmom_family.hpp`](../../src/internal/qmom_family.hpp)).

`compute_recurrence` runs the table of [Wheeler inversion](#wheeler-inversion)
on $2n+1$ moments: one more column in every row $k<n$, and row $n$, whose only
entry is $\sigma_{n,n}$, with $b_n=\sigma_{n,n}/\sigma_{n-1,n-1}$ and $\tau_n$.
The entries for $k<n$ are those of the $2n$-moment table. $a_n$ would need
$m_{2n+1}$; the closure (Fox & Laurent, 2022, eq. 4.16, or eq. 4.13 in
standardized moments) sets it:

```math
a_n=\frac{1}{n}\sum_{k=0}^{n-1}a_k .
```

The rule is the Gauss rule of $J_{n+1}$, with diagonal $a_0,\ldots,a_n$ and
off-diagonal $\sqrt{b_1},\ldots,\sqrt{b_n}$, the matrix associated with the
closure (Fox & Laurent, 2022, Section 7). For moments symmetric about their
mean, every $a_k$ is the mean, the closure is exact, and the rule is the
classical $(n+1)$-node Gauss rule of $m_0,\ldots,m_{2n}$ and the symmetric
$m_{2n+1}$. As $b_n\to0$, $J_{n+1}$ splits into $J_n$ and the entry $a_n$: the
rule tends to the classical $n$-node rule and a node at $a_n$ of vanishing
weight. On the boundary the moments determine the distribution, a sum of at
most $n$ weighted points (Fox & Laurent, 2022, Section 2.4), and the classical
rule is that distribution.

The input checks are those of
[Input checks and the cutoff](#input-checks-and-the-cutoff) with $2n+1$
moments, messages naming `node_count`, the message
`Moment sequence has fewer than 2*n+1 entries.`, and no one-node shortcut:
$n=1$ gives the two nodes of $J_2$. The rounding test runs to $k=n$:

| Moments | Adaptive | Fixed |
| --- | --- | --- |
| Interior: $\sigma_{k,k}>\tau_k$ for $k=1,\ldots,n$ | Rule of $J_{n+1}$, without `rmin` and `eabs` | Rule of $J_{n+1}$ |
| On or below the boundary: first $k$ with $\sigma_{k,k}\le\tau_k$ | The classical selection from $k$, exactly what `compute_quadrature(moments, n, options)` returns | `supports at most k nodes` or the negative $b$ |

Before the rule of $J_{n+1}$, both modes check $a_0,\ldots,a_n$ and
$b_1,\ldots,b_n$ as in [Fixed node count](#fixed-node-count), and a failed
eigensolve is the result. If $\sigma_{n,n}$ overflows, so does $\tau_n$.
Adaptive mode then takes the classical selection, since
$\sigma_{n,n}\le\tau_n$ holds; fixed mode skips the non-finite $\tau_n$ and
reports a non-finite recurrence, even where
`compute_quadrature(moments, n, options)`, which does not read $m_{2n}$, would
succeed.

Adaptive mode applies `rmin` and `eabs` only to the classical rules. Applied to
the rule of $J_{n+1}$, they would reduce interior moments to a rule outside the
closure: the classical $n$-node rule does not reproduce $m_{2n}$, and its
closure is QMOM's, which is only weakly hyperbolic (Fox & Laurent, 2022,
Theorem 4.3). They are not needed either. The rounding test already excludes
nodes made of rounding, as in fixed mode, and near the boundary the node of
small weight tends to $a_n$, a mean of well-determined coefficients, with a
weight that vanishes continuously. So every interior result is the HyQMOM rule,
and every boundary result its limit. A reduced result is always a classical
rule of at most $n$ nodes, and only on or below the boundary.

### Characteristic speeds

The moment system $\partial_t m_k+\partial_x m_{k+1}=0$, $k=0,\ldots,2n$,
closed by HyQMOM has the characteristic polynomial $P_{2n+1}=Q_nR_{n+1}$, with
$Q_n$ the $n$-th monic orthogonal polynomial of the recurrence and

```math
R_{n+1}=(x-a_n)\,Q_n-\beta_nQ_{n-1},\qquad \beta_n=\frac{2n+1}{n}\,b_n
```

(Fox & Laurent, 2022, Theorem 4.6, eqs. 4.11, 4.13 and 4.16), shown for
$n\le9$ and postulated beyond; the verification confirms it to $n=10$
([QMOM verification](../../tests/verification/qmom/README.md#hyqmom)). Its
$2n+1$ roots, the characteristic speeds, are the eigenvalues of $J_n$, diagonal
$a_0,\ldots,a_{n-1}$ and off-diagonal $\sqrt{b_1},\ldots,\sqrt{b_{n-1}}$, and
of $K_{n+1}$, which is $J_{n+1}$ with $\sqrt{\beta_n}$ in place of $\sqrt{b_n}$
(Fox & Laurent, 2022, eqs. 4.17 and 4.18); $\beta_n$ enters only the speeds,
never the rule. For $\beta_n>0$ the roots of $R_{n+1}$ strictly separate those
of $Q_n$ (Fox & Laurent, 2022, Theorem 4.8), so the speeds are distinct.

The public `compute_hyqmom_characteristic_speeds(state)`
([`src/public_api.cpp`](../../src/public_api.cpp)) follows the rule the state
holds, as the evolution does
([Phase-space evolution](phase_space_evolution.md#node-counts)), and
reconstructs nothing:

| State | Speeds |
| --- | --- |
| $n+1$ nodes | The eigenvalues of $J_n$ and $K_{n+1}$ for the stored moments $m_0,\ldots,m_{2n}$, the moments of the nodes, with the closure: $2n+1$ values |
| A classical rule of $k\le n$ nodes | Its $k$ abscissas, the roots of $Q_k$, each twice, and its closure $a_k$, the mean of $a_0,\ldots,a_{k-1}$, which is the trace of $J_k$ over $k$ and so the mean of the abscissas: $2k+1$ values |

The second row is the $k$-node system on its own boundary, $\beta_k=0$, where
$R_{k+1}=(x-a_k)\,Q_k$ (Fox & Laurent, 2022, Remark 4.9). At $k=n$ it is the
first row with $\beta_n=0$, so the speeds are continuous as the moments
approach the boundary. Fox & Laurent (2022) define no speeds on the boundary
for $k<n$; returning those of the $k$-node system is QMOMKit's choice.

The internal `compute_hyperbolic_characteristic_speeds` computes the first row
from the nodes $x_i$ and weights $w_i$. The first $2n+1$ moments of $n+1$
weighted points determine $a_0,\ldots,a_{n-1}$ and $b_1,\ldots,b_n$, the
recurrence of the points, which Lanczos on $\operatorname{diag}(x)$ from the
unit vector with components $\sqrt{w_i/m_0}$ computes by orthogonal
transformations. The moment table would lose them to cancellation near the
boundary: the moments of an accepted rule of $n+1$ nodes differ from those it
was built from by rounding, and can fail the rounding test. Lanczos removes
the projections onto all earlier vectors twice, which keeps them orthogonal in
floating point; once is not enough where they cancel most of the product with
$\operatorname{diag}(x)$, as for nodes of small spread far from the origin. As
nodes come together, the speeds tend continuously to the boundary values of
the second row, which coinciding nodes give to rounding. A $b_k$ that vanishes,
as for nodes all at zero or coordinates that underflow, or that is not finite
ends the recurrence, and the result is `not_realizable`. The closure is
`hyperbolic_closure`, shared with `compute_hyperbolic_quadrature`.
The options of the state do not enter, and no performance counter records the
speeds.

## Cases, results and statuses

`QmomCase` holds `moments`, `node_count`, `qmom_options` and
`abscissa_bounds`; `compute_qmom(case)` forwards all four. `compute_qmom(...)`
maps `QuadratureStatus` one to one onto `QmomStatus`:

| Status | Cause |
| --- | --- |
| `success` | A rule was accepted, in adaptive mode possibly with fewer nodes than requested. |
| `invalid_input` | Invalid options or bounds, bounds with the hyperbolic family, no nodes, too few or non-finite moments. |
| `not_realizable` | $m_0<0$, $m_0<$ `cutoff`, a fixed count the moments do not support, negative $b_k$, non-finite recurrence or eigenpairs, an unconverged eigensolve, or nodes outside the bounds. |

`ok()` is true only for `success`; the accepted node count is `weights.size()`.
Data returned with a failure is diagnostic, such as a one-node rule or the
rejected candidate, and empty for invalid input. Callers test the status, never
whether the vectors are empty.

## Where the code lives

| Function | File | Role |
| --- | --- | --- |
| `compute_qmom` | [`src/public_api.cpp`](../../src/public_api.cpp) | Public entry; family and bounds dispatch, status mapping. |
| `compute_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Public low-level Wheeler inversion, $s=0$. |
| `compute_wheeler_recurrence` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Public Wheeler recurrence and rounding bounds, $s=0$. |
| `compute_conditional_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Inversion of CQMOM conditional moments with scale $s$. |
| `compute_hyperbolic_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Internal HyQMOM reconstruction, $s=0$. |
| `compute_hyqmom_characteristic_speeds` | [`src/public_api.cpp`](../../src/public_api.cpp) | Public HyQMOM characteristic speeds of a state; family and count dispatch. |
| `compute_hyperbolic_characteristic_speeds` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Internal speeds of $n+1$ nodes, the eigenvalues of $J_n$ and $K_{n+1}$ from the Lanczos recurrence of the nodes. |
| `select_bounded_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Bounds acceptance for QMOM and CQMOM. |
| `cqmom_coordinate_scales` | [`src/coordinate_scales.cpp`](../../src/coordinate_scales.cpp) | Coordinate scales $s$. |
