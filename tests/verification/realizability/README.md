# Realizability verification

## Purpose

These cases verify `check_realizability`, which classifies a truncated
raw-moment sequence on the real line (Hamburger), on $\xi \ge 0$ (Stieltjes)
or on an interval $[a, b]$ (Hausdorff). They check the verdict and the
diagnostics that explain it, its invariance under shift, scale and number
density, and that rounding of the input may weaken a verdict to `boundary` or
`numerically_indeterminate` but never reverse it
([Realizability](../../../docs/developer/realizability.md)).

## Cases

All are `verification_realizability_analytical_<case>`
(`realizability_analytical.cpp`). Realizable references are analytical
densities (normal, gamma, lognormal, beta, uniform) or chosen weighted points;
non-realizable ones violate a stated condition.

- `diagnostic_success_cases`, `diagnostic_failure_cases`,
  `boundary_diagnostics`: checked order, pivot and message for a density on
  each support, and the limiting matrix for the failing matrices `m0`, `H_1`,
  `S_0` and `L_upper_0` and for boundary cases; the zero measure is `boundary`
  with `zero_measure`.
- `hamburger`, `stieltjes`, `hausdorff`, `hausdorff_on_other_intervals`:
  densities, weighted points and violations on each support, including
  odd-order range conditions, on $[0, 1]$, $[-1, 2]$ and $[2, 5]$. A variance
  0.25 above the Bhatia-Davis bound $(b - \mu)(\mu - a)$
  [[Bhatia & Davis, 2000](https://www.doi.org/10.1080/00029890.2000.12005203)]
  is rejected, one 0.25 below it accepted.
- `even_order_boundaries`: $n$ interior points are strictly realizable to order
  $2n - 1$ and on the boundary from $2n$; a point at an end of the support puts
  the sequence on the boundary once that end's localizing matrix is larger than
  the number of other points.
- `shift_and_scale_preserve_classification`,
  `resolution_limit_never_reverses_classification`: 21 sequences built from a
  location $c$ and a spread $s$ keep their verdict for $s$ from $10^{-9}$ to
  $10^{9}$ and $c/s$ up to a per-case shift of 2 to $10^3$, below the first at
  which double-precision input changes it; beyond, up to $c/s = 10^8$, it may
  only become `boundary` or `numerically_indeterminate`.
- `exact_point_sequences_are_never_rejected`: 50,000 seeded sequences of up to
  four points on a lattice of $1/4$ below 16 in magnitude, with weights in
  multiples of $1/8$, whose moments up to order 8 are exact in double precision.
  None may be rejected on the real line, and none of order $2n$ or more for $n$
  points, over 20,000 of them, may be called strictly realizable.
- `near_boundary_and_scaling`: the tolerance band below; strict classification
  for number densities from $10^{-200}$ to $10^{50}$, boundary and
  non-realizable classification at $10^{-50}$, and rejection at 0.6 times the
  largest double, where a matrix row sum overflows.

## Tolerances

Verdicts, and in the diagnostic cases the checked order and, for rejections
and boundary cases, the limiting matrix, must match exactly. There a strict verdict needs a pivot above the tolerance
used and a rejection one below its negative, so neither rests on rounding. The
classifier runs at its default relative tolerance $10^{-12}$; the near-boundary
cases set $10^{-8}$, so that moving $m_4$ of two points at $\pm 1$ by $10^{-10}$
is `numerically_indeterminate` and by $10^{-6}$ is decided, as are variances of
$\pm 10^{-10}$.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_realizability_analytical_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing_and_validation.md#building-and-running).
