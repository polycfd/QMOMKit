# QMOM verification

## Purpose

These cases verify the inversion of raw moments into a one-dimensional
quadrature, through `load_qmom` and `compute_qmom` and through the low-level
Wheeler inversion `compute_quadrature`, whose recurrence
`compute_wheeler_recurrence` reports, the internal HyQMOM reconstruction
`compute_hyperbolic_quadrature`, and the HyQMOM characteristic speeds. A rule
of $n$ weighted points is determined by its first $2n$ moments, and a HyQMOM
rule of $n + 1$ nodes and its speeds by $m_0, \ldots, m_{2n}$ and its closure,
so each reference has one correct answer.
They also check the contracts of [QMOM](../../../docs/developer/qmom.md): $m_0$
is the number density and is not normalized, a fixed node count is delivered or
refused, adaptive reduction returns only accepted rules, and abscissa bounds
reject a rule without clipping it.

## Cases

### Classical and manufactured rules

`verification_qmom_analytical_<case>` (`qmom_analytical.cpp`, `expected.hpp`,
`<case>.qmom`). Each fixture holds the first $2n$ raw moments of a rule listed
independently in `expected.hpp`:

| Cases | Reference rule |
| --- | --- |
| `single_dirac`, `two_node_discrete`, `three_node_discrete` | chosen weighted points; `single_dirac` has $m_0 = 2.5$ |
| `uniform_n2`, `uniform_n3` | Gauss-Legendre mapped to $[0, 1]$ |
| `exponential_n2`, `exponential_n3` | Gauss-Laguerre, $\alpha = 0$ |
| `gaussian_n2`, `gaussian_n3` | Gauss-Hermite for the weight $e^{-\xi^2/2}/\sqrt{2\pi}$ |
| `hyperbolic_gaussian_n1`, `hyperbolic_gaussian_n2`, `hyperbolic_normal_distribution_n2` | the 2- and 3-point Gauss-Hermite rules, from `qmom_family = hyperbolic` with node count 1 and 2 |

The sorted nodes and the weights must match the rule; the weights must be
positive, sum to $m_0$ and rebuild every moment. Uniform nodes must lie in
$[0, 1]$, exponential nodes must be positive, uniform and Gaussian rules must
be symmetric, and Gaussian odd moments must vanish. The hyperbolic fixtures
give the standard normal moments $m_0, \ldots, m_{2n}$, explicitly or, in
`hyperbolic_normal_distribution_n2`, through `[initial]` without
`maximum_order`, so that an inferred order below $2n$ would fail; the unit tests
pin the inferred order itself. For moments symmetric about their mean the
closure is exact, and the HyQMOM rule with node count $n$ is the $(n+1)$-point
Gauss rule.

### Wheeler inversion

`verification_wheeler_analytical_<case>` (`wheeler_analytical.cpp`). Rules of
one to six and of eight chosen weighted points must be recovered from their
moments, including the eight-point rule at number density 2500, three points
scaled by every fourth power of ten from $10^{-20}$ to $10^{20}$, and droplet
diameters of 20 to 200 micrometers given in meters
(`droplet_diameters_in_metres_are_recovered`).

### Wheeler recurrence

`verification_wheeler_recurrence_exact_<case>` (`wheeler_recurrence_exact.cpp`).
`compute_wheeler_recurrence` must return the recurrence of
[QMOM](../../../docs/developer/qmom.md#wheeler-inversion) as exact rational
arithmetic gives it, for rules of two to five weighted points
(`two_points` to `five_points`), a symmetric rule with $a_k = 0$
(`symmetric_three_points`) and a rule at number density 1024
(`three_points_at_number_density_1024`). Nodes and weights are short binary
fractions, so the moments are exact in double. The references for $a_k$, $b_k$
and $\sigma_{k,k}$ were computed from those moments with Python's `fractions`,
by the script in the test source, and are stored to 17 significant digits; the
coefficients of the two-point and symmetric rules are exact in double, those of
the other rules are not. Entry 0 must hold $a_0 = m_1/m_0$, $b_0 = 0$,
$\sigma_{0,0} = m_0$ and tolerance 0, and every $\tau_k$ must be finite, positive
and below $\sigma_{k,k}$, so that all $n$ points are supported.

### HyQMOM

`verification_hyperbolic_quadrature_exact_<case>`
(`hyperbolic_quadrature_exact.cpp`). `compute_hyperbolic_quadrature` must return
the Gauss rule of $J_{n+1}$ with the closure of
[QMOM families](../../../docs/developer/qmom_families.md#hyperbolic), for
$n = 1$ to 10 and three shapes: a symmetric binomial rule (`symmetric`), a
skewed geometric one at number density 2500 (`skewed`), and two binomial
clusters of masses 3/4 and 1/4 shifted off the origin (`bimodal`). Each shape
has more than $n + 1$ weighted points, so its moments are interior; the script
in the test source asserts that every $\sigma_{k,k}$ exceeds $\tau_k$ at least a
thousandfold. The moments $m_0, \ldots, m_{2n}$ are rounded to double, and the
references are computed from those doubles: the recurrence and the closure with
Python's `fractions`, the eigenpairs of $J_{n+1}$ with mpmath at 60 digits.
Fixed mode must return $n + 1$ nodes, and adaptive mode the identical rule.

`verification_hyperbolic_characteristic_speeds_exact_<case>`
(`hyperbolic_characteristic_speeds_exact.cpp`). The characteristic speeds of the
closed system, from
[QMOM families](../../../docs/developer/qmom_families.md#characteristic-speeds),
of the given moments, the moments the references are computed from. A state's
speeds are those of its stored moments, which the unit tests check bit for bit.
Three references:

- `normal`: normal moments with dyadic mean $\mu$ and standard deviation
  $\sigma$, $n = 1$ to 3, whose speeds are $\mu + \sigma s$ for the roots $s$
  of $He_n$ and of $R_{n+1}$ of the standard normal distribution.
- `boundary`: the moments of $k \le n$ dyadic weighted points, whose speeds are
  the points, each twice, and their mean.
- `skewed`, `bimodal`: the moments of the HyQMOM cases above,
  $n = 1$ to 10, against the eigenvalues of the Jacobian of the closed
  system. The script in the test source closes $m_{2n+1}$ in exact rational
  arithmetic with the reverse step of Fox & Laurent (2022, Algorithm 2),
  differentiates it exactly with dual numbers, and finds the roots of the
  characteristic polynomial with mpmath at 60 digits. It asserts that they
  agree with the eigenvalues of $J_n$ and $K_{n+1}$ and interlace; the test
  asserts that the speeds strictly increase.

### Moment-space boundary

`verification_qmom_moment_space_boundary_<case>` (`moment_space_boundary.cpp`)
and `verification_adaptive_qmom_adaptive_public_qmom_reconstruction`
(`adaptive_qmom.cpp`). Three-point rules approach a two-point rule as two nodes
coalesce (separation down to $10^{-10}$) or the middle weight vanishes (down to
$10^{-14}$); the six moments of a two-point rule lie on the boundary for three
nodes.

- `coalescing_nodes_preserve_supplied_moments`: a fixed three-node request is
  accepted or `not_realizable`. At separation $d = 10^{-2}$ it returns the
  three chosen points, positions within $100\varepsilon/d^2$ and weights within
  $100\varepsilon/d^3$, the rounding the close pair's conditioning allows.
- `adaptive_reduction_preserves_represented_moments`,
  `accepted_reduction_respects_weight_ratio`: adaptive inversion accepts both
  sweeps; in the weight sweep, rules of two or more nodes have
  $\min w / \max w$ above the option `rmin`.
- `boundary_crossing_matches_realizability`: raising $m_4$ gives three nodes and
  `strictly_realizable`; lowering it gives `not_realizable` from both
  `compute_qmom` and `check_realizability`.
- `fixed_count_does_not_reduce_on_the_boundary`: in five units and three number
  densities, a fixed request reports "supports at most 2 nodes"; adaptive
  inversion accepts fewer nodes.
- `adaptive_public_qmom_reconstruction`: adaptive inversion returns the two
  chosen points.

Every accepted rule of $n$ nodes must have non-negative weights summing to
$m_0$, sorted abscissas, and rebuild the first $2n$ moments. Where reduction
starts depends on conditioning, not on a contract, and is not asserted.

### Scaling and support

`verification_qmom_scaling_and_support_<case>` (`scaling_and_support.cpp`).
Seeded random rules of one to eight clustered weighted points, with weights
spanning six decades, are inverted with up to two nodes more than they have, so
many sequences lie on the boundary; both fixed-count refusals and adaptive
reductions occur.

- `power_of_two_scaling_is_exact`: scaling the coordinate by $2^p$ and the
  number density by $2^q$ scales every moment exactly, so the status, the node
  count, the abscissas times $2^p$ and the weights times $2^q$ must match bit
  for bit, in fixed and adaptive inversion. No decision may depend on the unit
  or the number density, except the absolute cutoff on $m_0$.
- `hyperbolic_power_of_two_scaling_is_exact`: the same for the hyperbolic
  family, whose recurrence coefficients, closure (their mean) and rounding
  bounds scale exactly. Rules of two to eight points are given a node count
  $n$ from 1 to the point count less one and $m_0, \ldots, m_{2n}$, so that
  most lie in the interior; rules of $n + 1$ nodes and classical fallbacks must
  both occur, and no rule may have more nodes than points.
- `no_more_nodes_than_weighted_points`: a fixed request for more nodes than
  points is `not_realizable`, and adaptive inversion returns at most as many
  nodes as points. No node may be made of rounding.

### Abscissa bounds

`verification_qmom_quadrature_bounds_<case>` (`quadrature_bounds.cpp`), with
bounds $[\ell, u]$. Nodes inside the bounds give the unbounded result
(`fixed_bounded_qmom_accepts_nodes_inside_interval`); a node outside gives
`not_realizable` with that node unclipped
(`fixed_bounded_qmom_rejects_nodes_outside_interval_without_clipping`);
adaptive inversion keeps the largest admissible count, here one node
(`adaptive_bounded_qmom_uses_highest_admissible_count`). At scales $10^{-6}$ to
$10^{6}$, a node half the margin $10^{-12}\max(|\ell|, |u|)$ beyond a bound is
accepted and one twice the margin beyond is rejected
(`bounded_qmom_acceptance_does_not_depend_on_the_unit`).

### Failure statuses

`verification_qmom_compute_failure_<case>` (`qmom_compute_failures.cpp`). Each
request must return the listed status and a message containing the fragment;
`invalid_input` returns no nodes. Requests are fixed-count unless marked
adaptive, because adaptive inversion would answer the negative variance with
one node.

| Case | Request | Status | Fragment |
| --- | --- | --- | --- |
| `negative_zeroth_moment` | $(-1, 0, 1, 0)$, 2 nodes | `not_realizable` | `m0 < 0` |
| `negative_variance` | $(1, 0, -1, 0)$, 2 nodes | `not_realizable` | `negative recurrence` |
| `insufficient_moments` | 3 moments, 2 nodes | `invalid_input` | `fewer than 2*n` |
| `insufficient_moments_bounded_adaptive` | the same, bounded and adaptive | `invalid_input` | `fewer than 2*n` |
| `non_finite_moment_bounded_adaptive` | $(1, 0.5, 0.5, \mathrm{NaN})$, 2 nodes, bounded and adaptive | `invalid_input` | `non-finite` |
| `zero_node_count` | 0 nodes | `invalid_input` | `positive` |

## Tolerances

| Cases | Tolerance | Metric | Why |
| --- | --- | --- | --- |
| Gauss and manufactured rules, the hyperbolic fixtures included | $10^{-12}$ | both absolute and relative to $\max(\lvert\text{reference}\rvert, 1)$ | observed $10^{-16}$ to $10^{-14}$ |
| Wheeler, up to three points | $10^{-14}$ to $10^{-11}$ per case | absolute | observed $10^{-15}$ to $10^{-13}$ |
| Wheeler, from four points and the scaled points | $10^{-12}$; $10^{-10}$ for eight points, $2 \times 10^{-12}$ for the droplets | relative to the largest $\lvert\xi\rvert$ and the total weight | conditioning worsens with $n$ and spread; observed $10^{-15}$ at four points to $10^{-12}$ at eight |
| Wheeler recurrence | $10^{-13}$ | $b_k$ and $\sigma_{k,k}$ relative to the reference; $a_k$ relative to the largest $\lvert\xi\rvert$, since $a_k$ lies between the smallest and largest node and vanishes for the symmetric rule | observed at most 13 $\varepsilon$ ($2.9 \times 10^{-15}$, five points); the rules with exact coefficients are reproduced exactly |
| HyQMOM | $10^{-12}$ up to $n = 7$, $10^{-11}$ for $n = 8$, $5 \times 10^{-10}$ for $n = 9$ and 10 | relative to the largest $\lvert\xi\rvert$ and the total weight | conditioning worsens with $n$ |
| HyQMOM characteristic speeds | as for HyQMOM; $10^{-12}$ for `normal` and `boundary` | relative to the largest $\lvert\text{speed}\rvert$ of the reference | conditioning worsens with $n$; a safety bound for `normal` and `boundary` |
| Boundary, adaptive and bounds | $10^{-12}$ | absolute; the boundary sweeps relative to $\lvert m_k\rvert$ where it exceeds one | a safety bound over observed errors of $10^{-15}$ or zero |
| Scaling and support | none | exact equality of status, node count, abscissas and weights | power-of-two scalings are exact; scalings that would leave the normal double range are skipped |
| Failure statuses | none | exact status and message fragment | |

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_(qmom_(analytical|compute_failure|moment_space_boundary|quadrature_bounds|scaling_and_support)|wheeler_analytical|wheeler_recurrence_exact|hyperbolic_quadrature_exact|hyperbolic_characteristic_speeds_exact|adaptive_qmom)_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing.md#building-and-running).
