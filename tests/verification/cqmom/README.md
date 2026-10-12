# CQMOM verification

## Purpose

These cases verify that CQMOM recovers a quadrature known in advance:
positions, weights, conditioning paths and branch node counts, stored in
phase-space order for any conditioning order, in fixed and adaptive mode, with
rebuilt moments equal to the given ones. Further cases cover moment
reconstruction, abscissa bounds and branch growth. See
[CQMOM](../../../docs/developer/cqmom.md) for the recursion and
[QMOM](../../../docs/developer/qmom.md) for the inversion on each branch.

Every reference is a set of weighted points whose raw moments
$M_k = \sum_q W_q \prod_d \xi_{q,d}^{k_d}$ are summed by the test or written
into the `.qmom` fixture, never computed by QMOMKit; the expected result is the
set of points itself. `recursive_quadrature.cpp` and `recursive_cqmom.cpp` call
the internal recursion and inspect its traversal tree; all other cases use the
public API. Failures of `compute_cqmom(...)` are covered in
[`cqmom_compute_failures.cpp`](../../component_regression/cqmom_compute_failures.cpp).

## Cases

### Two dimensions

`verification_cqmom_{independent_product,deterministic_affine,branch_dependent}`
(`cqmom_2d_analytical.cpp`): node counts {2, 2}, conditioning order {0, 1}.

| Case | Rule | Adaptive | Branch counts |
| --- | --- | --- | --- |
| `independent_product` | $\xi_0 \in \{1, 3\}$, weights {0.4, 0.6}; independently $\xi_1 \in \{-1, 2\}$, weights {0.25, 0.75} | no | {2, 2} |
| `deterministic_affine` | same $\xi_0$; $\xi_1 = 2\xi_0 + 1$ | yes | {1, 1} |
| `branch_dependent` | $\xi_0 \in \{1, 3\}$, weights {0.35, 0.65}; $\xi_1 = -0.5$, or $\{1, 4\}$ with weights {0.3, 0.7} | yes | {1, 2} |

Checked: parent and conditional positions and weights, terminal paths, positive
weights summing to $M_{0,0} = 1$, node counts in the diagnostics, and the
rebuilt required moments.

### Three to six dimensions

`verification_cqmom_{affine_chain_5d,independent_product_[3-6]d,mixed_affine_independent_4d}`
(`cqmom_high_dimensional_analytical.cpp`): two nodes per coordinate. The
configured moments are first checked against the reference rule.

- `independent_product_3d` to `_6d`, fixed: products of two-point rules, with
  conditioning orders {1, 0, 2}, {2, 0, 3, 1}, {4, 1, 3, 0, 2} and
  {5, 2, 0, 4, 1, 3}. Every branch must recover the same rule for its
  coordinate.
- `affine_chain_5d`, adaptive: $\xi_0 \in \{-0.5, 1.5\}$, each later coordinate
  affine in the one before, so every branch below the root reduces to one node.
- `mixed_affine_independent_4d`, adaptive, order {0, 2, 1, 3}: independent
  two-point $\xi_0$ and $\xi_2$, with $\xi_1 = 2\xi_0 + 1$ and
  $\xi_3 = -\xi_1 + 0.5\xi_2$; only the affine levels reduce.

Checked: terminal positions, global and local weights, branch counts, the affine
relations, positive weights summing to the reference number density, and the
full rebuilt tensor.

### Coordinate units

`verification_cqmom_{unit_invariance,disparate_units}` (`unit_invariance.cpp`)
use a ragged adaptive five-point rule in three dimensions with order {2, 0, 1},
with coordinate $d$ in units scaled by $s_d$: $s = (1, 1, 1)$ and
$(10^6, 10^{-6}, 10^3)$ from the two `unit_invariance_*_3d.qmom` fixtures, and
$(10^8, 10^{-8}, 10^4)$ summed in the test, with moments from about $10^{-24}$
to $10^{37}$. Checked: paths, positions divided by $s$, weights and the required
moments, plus branch counts for the fixtures. The adaptive rounding bound scales
with the coordinate
([QMOM](../../../docs/developer/qmom.md#rounding-bound-on-the-recurrence)), so
units must not change the result.

### Public workflows

`verification_cqmom_public_workflows` (`cqmom_public_workflows.cpp`) loads
`independent_product.qmom` as the default `CqmomPreparedCase` and as the open
`CqmomCase`. Both results are compared with the reference rule, since agreement
alone would pass if both were wrong, and with each other. The `CqmomCase` must
hold exactly the ten required moments in their documented order.

### Moment reconstruction

`verification_cqmom_moment_reconstruction_*` (`moment_reconstruction.cpp`):
`qmom_all_orders` for a three-point rule; `cqmom_selected`, where the request
{2,0}, {0,0}, {1,1}, {0,2}, {2,0}, {3,1} keeps its order and duplicate;
`cqmom_required_layout`, aligned with `layout.required_indices`; and
`failure_paths`, where an unsuccessful 1D quadrature, a wrong-rank request and a
mismatched hierarchy each throw.

### Recursive inversion

`verification_recursive_quadrature_manufactured_*` (`recursive_quadrature.cpp`)
covers one to three dimensions, both 2D orders, node counts 2 and 3, and
adaptive branches with fewer points than requested, mostly through the rebuilt
tensor. In one dimension the result must equal Wheeler. The
`..._velocity_collapse_...` case uses non-realizable moments {1, 0, -1, 0} for
$\xi_1$ instead of a rule and checks only its reduction to one node in
either conditioning order.

`verification_cqmom_recursive_*` (`recursive_cqmom.cpp`) compares positions and
weights in six fixed cases, `1d` to `4d_adversarial_non_identity` (clustered
positions, a weight of 0.09, a check of path-to-storage mapping), and in the
eight-point ragged `4d_adversarial_adaptive_non_identity`. `dimension_sweep_D2`
to `D6` compare moments only.

### Abscissa bounds

`verification_cqmom_quadrature_bounds_*` (`quadrature_bounds.cpp`): a rule
inside its bounds is returned unchanged; one with a point outside fails as
`qmom_failed`, the point not moved onto the bound; and in 2D with order {1, 0}
each inversion uses the bound of its own coordinate.

### Branch growth

`verification_cqmom_adaptive_branch_growth_adaptive_branch_reinjection`
(`adaptive_branch_growth.cpp`): a constant source carries one $\xi_1$ per
$\xi_0$ branch to two on the branch $\xi_0 = 1$ over unit time. Checked: branch
counts, paths, the recorded growth, positions, weights and moments
([active branch counts](../../../docs/developer/cqmom.md#active-branch-counts)).

## Tolerances

Errors are absolute unless relative, $|a - e| / \max(|e|, 1)$.

| Cases | Tolerance | Metric | Why |
| --- | --- | --- | --- |
| Two dimensions | $10^{-12}$ | positions, weights; relative moments | observed $10^{-15}$; moments reach $10^3$ |
| Three to six dimensions | $10^{-11}$ | positions, weights; relative moments | error grows from $10^{-15}$ in 3D to $10^{-13}$ in 6D; the bound sits one order above |
| Coordinate units | $2 \times 10^{-10}$; `disparate_units` $10^{-10}$ | relative positions divided by $s$; weights; relative moments | safety bounds, observed $10^{-16}$; moments below 1 are checked only absolutely |
| Workflows, reconstruction, bounds | $10^{-12}$ | positions, weights, moments | observed $10^{-15}$; moments below $10^2$ |
| Recursive inversion | $10^{-12}$ to $10^{-7}$ per case | positions, weights, moments | safety bounds; the loosest cover 3D moments with errors of $10^{-11}$ |
| Branch growth | $5 \times 10^{-11}$ | positions, weights, moments | safety bound: SSP-RK3 integrates a constant rate exactly; observed $10^{-16}$ |

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_(cqmom_|recursive_quadrature_manufactured_)' --output-on-failure
```

See [Building and running](../../../docs/developer/testing_and_validation.md#building-and-running).
