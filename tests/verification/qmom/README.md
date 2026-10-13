# QMOM verification

## Purpose

These cases verify the inversion of raw moments into a one-dimensional
quadrature, through `load_qmom` and `compute_qmom` and through the low-level
Wheeler inversion `compute_quadrature`. A rule of $n$ weighted points is
determined by its first $2n$ moments, so each reference has one correct answer.
They also check the contracts of [QMOM](../../../docs/developer/qmom.md): $m_0$
is the number density and is not normalized, a fixed node count is delivered or
refused, adaptive reduction returns only accepted rules, and abscissa bounds
reject a rule without clipping it.

## Cases

**Classical and manufactured rules**, `verification_qmom_analytical_<case>`
(`qmom_analytical.cpp`, `expected.hpp`, `<case>.qmom`). Each fixture holds the
first $2N$ raw moments of a rule listed independently in `expected.hpp`:

| Cases | Reference rule |
| --- | --- |
| `single_dirac`, `two_node_discrete`, `three_node_discrete` | chosen weighted points; `single_dirac` has $m_0 = 2.5$ |
| `uniform_n2`, `uniform_n3` | Gauss-Legendre mapped to $[0, 1]$ |
| `exponential_n2`, `exponential_n3` | Gauss-Laguerre, $\alpha = 0$ |
| `gaussian_n2`, `gaussian_n3` | Gauss-Hermite for the weight $e^{-\xi^2/2}/\sqrt{2\pi}$ |

The sorted nodes and the weights must match the rule; the weights must be
positive, sum to $m_0$ and rebuild every moment. Uniform nodes must lie in
$[0, 1]$, exponential nodes must be positive, uniform and Gaussian rules must
be symmetric, and Gaussian odd moments must vanish.

**Wheeler inversion**, `verification_wheeler_analytical_<case>`
(`wheeler_analytical.cpp`). Rules of one to six and of eight chosen weighted
points must be recovered from their moments, including the eight-point rule at
number density 2500, three points scaled by every fourth power of ten from
$10^{-20}$ to $10^{20}$, and droplet diameters of 20 to 200 micrometers given in
meters (`droplet_diameters_in_metres_are_recovered`). Reordering the input
points must change neither the moments nor the sorted output.

**Moment-space boundary**, `verification_qmom_moment_space_boundary_<case>`
(`moment_space_boundary.cpp`) and
`verification_adaptive_qmom_adaptive_public_qmom_reconstruction`
(`adaptive_qmom.cpp`). Three-point rules approach a two-point rule as two nodes
coalesce (separation down to $10^{-10}$) or the middle weight vanishes (down to
$10^{-14}$); the six moments of a two-point rule lie on the boundary for three
nodes.

- `coalescing_nodes_preserve_supplied_moments`: a fixed three-node request is
  accepted or `not_realizable`, and returns three nodes at separation
  $10^{-2}$.
- `adaptive_reduction_preserves_represented_moments`,
  `accepted_reduction_respects_weight_ratio`: adaptive inversion accepts both
  sweeps; in the weight sweep, rules of two or more nodes have
  $\min w / \max w > r_{\min}$.
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

**Scaling and support**, `verification_qmom_scaling_and_support_<case>`
(`scaling_and_support.cpp`). Seeded random rules of one to eight clustered
weighted points, with weights spanning six decades, are inverted with up to two
nodes more than they have, so many sequences lie on the boundary; both
fixed-count refusals and adaptive reductions occur.

- `power_of_two_scaling_is_exact`: scaling the coordinate by $2^p$ and the
  number density by $2^q$ scales every moment exactly, so the status, the node
  count, the abscissas times $2^p$ and the weights times $2^q$ must match bit
  for bit, in fixed and adaptive inversion. No decision may depend on the unit
  or the number density.
- `no_more_nodes_than_weighted_points`: a fixed request for more nodes than
  points is `not_realizable`, and adaptive inversion returns at most as many
  nodes as points. No node may be made of rounding.

**Abscissa bounds**, `verification_qmom_quadrature_bounds_<case>`
(`quadrature_bounds.cpp`). Nodes inside the bounds give the unbounded result
(`fixed_bounded_qmom_accepts_nodes_inside_interval`); a node outside gives
`not_realizable` with that node unclipped
(`fixed_bounded_qmom_rejects_nodes_outside_interval_without_clipping`);
adaptive inversion keeps the largest admissible count, here one node
(`adaptive_bounded_qmom_uses_highest_admissible_count`). At scales $10^{-6}$ to
$10^{6}$, a node half the margin $10^{-12}\max(|\ell|, |u|)$ beyond a bound is
accepted and one twice the margin beyond is rejected
(`bounded_qmom_acceptance_does_not_depend_on_the_unit`).

**Failure statuses**, `verification_qmom_compute_failure_<case>`
(`qmom_compute_failures.cpp`). Fixed-count requests, except the bounded
adaptive ones, since adaptive inversion answers a negative variance with one
node, must return the status and a message containing the fragment;
`invalid_input` returns no nodes.

| Case | Request | Status | Fragment |
| --- | --- | --- | --- |
| `negative_zeroth_moment` | $(-1, 0, 1, 0)$, 2 nodes | `not_realizable` | `m0 < 0` |
| `negative_variance` | $(1, 0, -1, 0)$, 2 nodes | `not_realizable` | `negative recurrence` |
| `insufficient_moments` | 3 moments, 2 nodes | `invalid_input` | `fewer than 2*n` |
| `insufficient_moments_bounded_adaptive` | the same, bounded and adaptive | `invalid_input` | `fewer than 2*n` |
| `non_finite_moment_bounded_adaptive` | $(1, 0.5, 0.5, \mathrm{NaN})$, 2 nodes, bounded and adaptive | `invalid_input` | `non-finite` |
| `zero_node_count` | 0 nodes | `invalid_input` | `positive` |

## Tolerances

- Classical rules: $10^{-12}$, absolute and relative to
  $\max(|\text{reference}|, 1)$; observed $10^{-16}$ to $10^{-14}$.
- Wheeler up to three points: absolute $10^{-14}$ to $10^{-11}$ per case,
  observed $10^{-15}$ to $10^{-13}$. From four points, relative to the largest
  $|\xi|$ and the total weight: $10^{-12}$, but $10^{-10}$ for eight points and
  $2 \times 10^{-12}$ for the droplets, because conditioning worsens with $n$
  and spread; observed $10^{-15}$ at four points to $10^{-12}$ at eight.
- Boundary, adaptive and bounds: $10^{-12}$, a safety bound over observed
  errors of $10^{-15}$ or zero.
- Scaling and support: exact equality and node counts, no tolerance. Scalings
  that would take a moment out of the normal double range are skipped, because
  they are not exact.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_(qmom_(analytical|compute_failure|moment_space_boundary|quadrature_bounds|scaling_and_support)|wheeler_analytical|adaptive_qmom)_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing.md#building-and-running).
