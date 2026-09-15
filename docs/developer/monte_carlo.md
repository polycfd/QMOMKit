# Monte Carlo companion package

This page documents the development and validation conventions for the optional
QMOMKitMonteCarlo companion package under `tools/monte_carlo/`.

For user-facing build and run instructions, see
[docs/user/monte_carlo.md](../user/monte_carlo.md). For general numerical
validation policy, see [validation.md](validation.md). For test category
guidance, see [testing.md](testing.md).

## Package role

QMOMKitMonteCarlo provides finite ensemble sampling, ensemble moment
estimation, advection sample evolution, text output, and comparison utilities.
It is a companion package for sample-based analysis and stochastic validation
references.

It is not part of the core `qmomkit::qmomkit` library. It must continue to use
an installed QMOMKit package rather than depending on core build-tree details.
The companion target is:

```text
qmomkit_monte_carlo::qmomkit_monte_carlo
```

The package should remain scientifically aligned with QMOMKit conventions:
absolute raw moments of the number density function, phase-space dimension
order for sample coordinates and tensor indices, explicit status values, and
deterministic output ordering.

## Reference construction

Monte Carlo validation compares raw moment histories:

```text
QMOMKit deterministic moments
           |
        versus
           |
Monte Carlo estimated moments
```

It does not validate by comparing quadrature nodes, internal CQMOM conditioning
paths, or branch-local structures. Those quantities may be useful diagnostics,
but raw moments are the primary comparison quantities.

For a direct QMOMKit-versus-Monte-Carlo validation case, prefer starting both
paths from the same finite realization:

```text
sampled ensemble
    -> Monte Carlo sample evolution

sample-derived raw moments
    -> QMOM or CQMOM evolution
    -> reconstructed raw moments
```

This separates representation and evolution differences from initial
finite-sample discrepancy. If the deterministic path starts from analytical
moments while the Monte Carlo path starts from a finite sample, the comparison
also includes sampling error at the initial time.

## Moment conventions

For a multidimensional multi-index $\alpha$:
$$M_{\alpha} = \int \prod_{d} \xi_{d}^{\alpha_{d}} n(\xi) \ \mathrm{d} \xi$$

For deterministic quadrature:
$$M_{\alpha} = \sum_{i} W_i \prod_d \xi_{i,d}^{\alpha_d}$$

For a Monte Carlo ensemble with $N$ samples and total number density $n_0$:
$$M_\alpha = \frac{n_0}{N} \sum_p \prod_d \xi_{p,d}^{\alpha_d}$$

These are unnormalized raw moments of the number density function. They are not
centered, standardized, divided by $M[0,...,0]$, or converted to probability
moments unless a case explicitly defines a unit-number-density problem.

Sample coordinates use QMOMKit phase-space dimension order:
$$\xi = \left(\xi_0, \xi_1, ..., \xi_{D-1}\right)$$

Coordinate names are metadata and must not change sampling, storage, evolution,
moment-index, or output order.

### Population birth/injection

`inject_population(...)` is a local equally weighted birth operation. The
caller supplies the integrated number-density amount $\Delta n$ for its own
chosen interval and a callback sampling the normalized birth distribution. The
companion package does not integrate a physical birth-rate law or own timing or
operator splitting.

For a valid existing ensemble with $N$ samples and positive number density
$n_0$, its representative weight is $w = n_0 / N$. The operation uses one local
`std::mt19937_64` seeded by the public seed argument to draw
$K \sim \mathrm{Poisson}(\Delta n / w)$ and then to sample all $K$ birth coordinates.
It returns $n_{0,\mathrm{new}} = w (N + K)$. Consequently the realized post-birth
estimator is:
$$M_{\alpha,\mathrm{new}} = w \left[\sum_p \xi_{p}^{\alpha} + \sum_r Y_{r}^{\alpha}\right]$$

and its expected increment equals the caller's integrated birth contribution.
Samples remain unweighted coordinate vectors; variable per-sample weights are
not introduced. A deterministic callback may ignore the engine. A zero birth
amount is a successful no-op. Positive birth into an empty or zero-density
ensemble is deliberately unsupported because no positive representative weight
can be inferred.

### Selective multidimensional estimation

`compute_requested_ensemble_moments(...)` estimates only the indices in a
`qmomkit::RawMomentRequest`. It preserves request order and duplicate entries,
and returns aligned `indices`, `values`, and `standard_errors` in
`SelectedEnsembleMomentEstimate`. An empty request is a successful empty result
after the ensemble and number-density inputs pass validation.

The implementation evaluates each requested monomial directly with work
proportional to the number of requested indices, samples, and phase-space
dimensions. It must not construct maximum orders, dense extents, a dense
tensor, or a shared power table. Its normalization and standard-error
convention are identical to the dense ensemble moment estimator.

Validation is ordered so ensemble structure and sample finiteness are checked
first, followed by requested-index rank, number density, and exponent overflow.
Failures return the existing `EnsembleMomentStatus` values with an empty
payload; a successful result is the only result whose values are usable as
estimated moments. The selective API is an additive Monte Carlo-package API and
does not alter the existing dense overloads.

## Reproducibility

A Monte Carlo validation case should record enough information to reproduce and
interpret the comparison:

- random seed;
- sample count;
- number density;
- distribution or explicit quadrature used to construct the reference;
- requested moment orders or tensor extents;
- phase-space coordinate order;
- time interval and output times, when evolution is used;
- RK23 options, when sample evolution is used;
- tolerance or uncertainty rule used by the validation check.

Sampling currently uses an explicit seed and sample count. Fixed seeds make
repository tests repeatable on a given implementation, but standard
distribution adapters are not guaranteed to produce bit-for-bit identical
samples on every platform. Tests should use mathematical tolerances or
uncertainty checks rather than exact sampled values.

## Comparison quantities

Use explicit reference and computed histories. For each compared moment:

```text
signed_error = computed - reference
absolute_error = abs(signed_error)
relative_error = absolute_error / scale
z_score = signed_error / standard_error
```

For C++ comparison utilities, the first `compare_moments(...)` argument is the
reference and the second argument is the computed Monte Carlo estimate. For the
Python history-comparison script, the first input file is the reference table
and the second input file is the computed table.

The relative-error scale should include a small positive floor for near-zero
reference values. When Monte Carlo standard errors are available, they may also
be used as part of the scale so nearly zero references with visible sampling
uncertainty remain interpretable.

If the standard error is zero, the z-score is zero for exact agreement and
signed infinity for nonzero signed error. Large z-scores are diagnostic
signals; they are not automatic proof that either calculation is wrong.

## Output format

Normal Monte Carlo output follows the QMOMKit text-table convention:

- `.txt` filenames;
- whitespace-delimited columns;
- one commented data header line beginning with `#`;
- optional metadata comments before the data header;
- 17-digit floating-point text precision for reproducible scientific checks.

Example moment-history table:

```text
# time sample_count M0 M1 M2
0.0 100000 1.0 2.0 5.0
```

One-dimensional moment labels use `M0`, `M1`, and so on. Tensor moment labels
use compact multi-index labels such as `M00`, `M01`, and `M10` in QMOMKit
phase-space order. Comparison tables should include both `signed_error` and
`absolute_error`; a signed quantity must not be named `absolute_error`.

## Tolerances and uncertainty

Deterministic validation cases should use exact or manufactured expected values
and fixed tolerances. Monte Carlo validation should account for sampling
uncertainty through standard errors, z-scores, sample-count studies, or another
documented statistical rule.

Tolerance choice should consider:

- sample count;
- moment order;
- number density and moment scale;
- distribution tails;
- time integration tolerances;
- finite-node capacity in the deterministic representation.

Do not relax tolerances simply to match a previous output file. If a tolerance
changes, document the reason in the relevant test or validation description.

## Adding a Monte Carlo validation case

Add a Monte Carlo validation case when a stochastic ensemble reference gives
useful evidence that deterministic QMOM or CQMOM moment evolution behaves as
intended.

Before adding the case, define:

- what capability is being checked;
- how the reference ensemble is sampled or constructed;
- which raw moments are compared;
- whether both paths start from the same finite realization;
- the sample count and seed;
- the time grid and evolution options;
- the comparison metrics and tolerance rule;
- which diagnostics are reported but not used as pass/fail criteria.

Use deterministic analytical or manufactured validation instead when an exact
reference is available. Monte Carlo evidence complements deterministic
validation; it does not replace it.