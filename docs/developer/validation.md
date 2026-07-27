# Validation

This page documents the QMOMKit numerical validation strategy for maintainers,
contributors, and AI-assisted development workflows. It answers how QMOMKit
establishes that its numerical methods are correct and trustworthy.

For software test architecture, see `testing.md`. For raw-moment, tensor, node,
weight, and phase-space conventions, see `numerical_conventions.md`. For
algorithm-specific contracts, see `qmom_wheeler.md` and `cqmom_recursive.md`.
User-facing Monte Carlo workflow guidance lives in `docs/user/monte_carlo.md`.
Developer guidance for the companion package lives in `monte_carlo.md`.

## Validation Philosophy

Numerical correctness requires more than passing unit tests. Unit tests protect
local contracts, failure handling, and regression behavior, but validation asks
whether the numerical method reproduces behavior that is independently known or
scientifically meaningful.

The distinction is:

```text
testing:
    software correctness

validation:
    scientific and numerical credibility
```

Validation compares QMOMKit behavior against known expectations:

- exact analytical solutions;
- manufactured moment or quadrature cases;
- deterministic convergence behavior;
- independent stochastic Monte Carlo references.

Deterministic validation cases are essential. They provide repeatable evidence
that core algorithms preserve the mathematical conventions they claim to
implement. Stochastic validation is useful when an independent ensemble
reference is valuable, but it must be interpreted statistically rather than as
exact pass/fail equality.

## Validation Hierarchy

### Algorithm Verification

Algorithm verification validates an individual numerical algorithm against a
known mathematical object.

Examples include:

- reconstructing every supplied raw moment from a recovered Wheeler quadrature;
- comparing one-dimensional QMOM nodes and weights against known Gaussian
  quadrature rules or manufactured discrete rules;
- verifying CQMOM conditional moments from an independently defined weighted
  Vandermonde problem;
- reconstructing complete raw moment tensors from terminal CQMOM quadrature;
- checking realizability classifications against known truncated moment
  domains.

Algorithm verification should isolate the algorithm as much as practical. If a
case verifies Wheeler inversion, the expected nodes, weights, and reconstructed
moments should be defined independently of the Wheeler implementation. If a case
verifies recursive CQMOM, it should make terminal weights, phase-space
coordinates, conditioning paths, and raw tensor reconstruction visible.

### Numerical Method Verification

Numerical method verification validates a composition of numerical methods,
typically with time evolution or operator splitting.

Examples include:

- RK23 behavior against simple ODEs with known exact solutions;
- SSP-RK3 source evolution against closed-form source updates;
- node advection with constant, affine, or linear phase-space velocities;
- moment-source integration with known moment-space dynamics;
- transport composition against commuting and noncommuting analytical cases;
- convergence studies that verify expected error reduction.

These cases should define the exact or manufactured reference independently of
the production integration path. They should report the error metric, tolerance,
and expected qualitative behavior, such as conservation of `m_0`, unchanged node
weights during node advection, or second-order behavior for Strang splitting in
a noncommuting case.

### Application-Level Validation

Application-level validation checks whether a complete QMOMKit workflow behaves
credibly for a physically meaningful or workflow-level setup. QMOMKit remains a
numerical C++ library, so these cases should not introduce CFD meshes, boundary
conditions, MPI behavior, or application-specific source models into the core
library.

Examples include:

- distribution-based initialization followed by QMOM or CQMOM reconstruction;
- evolution of moments generated from meaningful analytical distributions;
- sample-based comparisons between deterministic QMOMKit evolution and an
  independently evolved Monte Carlo ensemble;
- high-dimensional CQMOM workflows that exercise conditioning order,
  branch-dependent quadrature, and raw tensor reconstruction.

Application-level validation may combine multiple features, but each case should
state what capability is being validated and which observed quantities are
primary evidence rather than diagnostics.

## Analytical And Manufactured Solutions

Analytical references are preferred when an exact solution is available. They
are useful because the expected behavior can be derived without calling the
production algorithm under test.

For QMOM, analytical references include exact discrete quadratures and classical
one-dimensional quadrature rules for known moment families. Useful checks
include node and weight values, positivity, total mass preservation, and
reconstruction of all configured raw moments.

For CQMOM, analytical and manufactured references should make the conditional
structure explicit. Independent-product cases verify marginal and terminal
weight products. Deterministic affine cases verify adaptive branch reduction and
conditioning-path behavior. Higher-dimensional cases verify phase-space order,
conditioning order, terminal coordinates, global weights, and full raw tensor
reconstruction.

For transport, analytical references should isolate the numerical behavior being
validated. Constant or affine node advection tests exact phase-space evolution.
Closed-form source equations test moment-source integration. Commuting and
noncommuting operator cases test transport composition and convergence behavior.

Manufactured solutions are appropriate when the exact physical model is less
important than a controlled mathematical reference. A manufactured case may
start from exact nodes, weights, source equations, conditioning paths, or
terminal quadrature and then generate the raw moments that QMOMKit should
recover.

## Monte Carlo Validation

`tools/monte_carlo/` contains the optional QMOMKitMonteCarlo companion package.
It provides an independent stochastic ensemble reference for validation
workflows without adding headers to the core `include/qmomkit/` API.

Monte Carlo validation compares raw moment histories:

```text
QMOMKit deterministic moments
        versus
Monte Carlo estimated moments
```

It does not validate by comparing quadrature nodes, CQMOM conditioning paths, or
internal branch structures. Those quantities are diagnostics. The primary
comparison is the moment representation.

Monte Carlo references are stochastic sample-based references that validate
moment representations and their evolution. They are not production physical
models, analytical solutions, source-term models, or replacements for
QMOM/CQMOM transport. A Monte Carlo case supplies a sampled ensemble reference
whose uncertainty depends on sample count, seed, moment order, distribution
tails, and the selected dynamics.

For direct QMOMKit-versus-Monte-Carlo comparisons, initialize both paths from
the same finite realization when practical:

```text
sampled ensemble
    -> Monte Carlo evolution

sample-derived raw moments
    -> QMOM or CQMOM evolution
    -> reconstructed raw moments
```

This avoids mixing initial finite-sample discrepancy with representation or
evolution error. If QMOMKit starts from analytical moments while Monte Carlo
starts from a finite ensemble, the comparison includes sampling discrepancy at
the initial time.

Monte Carlo comparisons should record the sample count, seed, number density,
moment selection, time grid, evolution options, and uncertainty estimates.
Standard errors and z-scores are interpretation aids. A large z-score is a
diagnostic signal, not automatically a QMOMKit failure.

For the detailed Monte Carlo reference-data and output conventions, see
`monte_carlo.md`.

## Validation Conventions

Validation cases use the raw moment conventions in `numerical_conventions.md`.
For a multidimensional multi-index `alpha`:

```text
M_alpha = integral product_d xi_d^alpha_d n(xi) dxi
```

For deterministic quadrature:

```text
M_alpha = sum_i W_i product_d xi_(i,d)^alpha_d
```

For a Monte Carlo ensemble with `N` samples and total number density `n0`:

```text
M_alpha = n0 / N * sum_p product_d xi_(p,d)^alpha_d
```

These are unnormalized raw moments of the number density function. They are not
centered, standardized, divided by `M[0,...,0]`, or converted to probability
moments unless a validation case explicitly defines a unit-number-density
problem.

Relative errors should be interpreted with scale in mind. Near-zero references
need an absolute tolerance or a denominator floor. High-order moments may be
sensitive to tail behavior, roundoff, sampling noise, and finite-node capacity.

Deterministic references should normally use exact or manufactured expected
values and fixed tolerances. Stochastic references should include uncertainty
information and should distinguish deterministic numerical error from sampling
variation.

Validation output should preserve enough information to reproduce the
comparison:

- raw moment labels or tensor extents;
- phase-space coordinate order;
- total number density;
- node counts and conditioning order;
- adaptive settings, when used;
- reference construction;
- time interval and output times, when evolution is used;
- numerical tolerances and error metrics;
- Monte Carlo seed, sample count, and standard errors, when applicable.

## Adding Validation Cases

When adding a new numerical method or extending an existing one, define the
validation contract before relying on the result.

The contract should specify:

- the capability being validated;
- the known reference and how it is constructed;
- the raw moments or tensor moments being compared;
- expected accuracy and tolerance rationale;
- deterministic regression coverage for the core behavior;
- diagnostics that should be reported but are not primary pass criteria;
- assumptions about support, conditioning order, node counts, time integration,
  or adaptive behavior.

Add deterministic validation whenever a known analytical or manufactured
reference is available. Add Monte Carlo validation when an independent sampled
reference gives useful evidence for a representation or evolution workflow that
is not fully covered by deterministic cases.

Do not use benchmark timings as validation evidence. Benchmarks measure
performance; validation establishes numerical credibility.
