# Distribution Design

This page documents the current QMOMKit design for analytical distribution
handling. It is a developer reference for how distribution specifications are
represented, validated, and converted into raw moments for QMOM and CQMOM
initialization.

For raw moments, number density function (NDF) conventions, tensor storage, and
phase-space ordering, see [numerical_conventions.md](numerical_conventions.md).
For configuration-file mapping, see
[configuration_design.md](configuration_design.md). User-facing examples and
workflow guidance live in `docs/user/distributions.md`.

## Role Of Distributions

Distribution support is an initialization layer. A distribution describes a
compact analytical source for moments; QMOM and CQMOM algorithms consume the
generated moments, not the distribution object itself.

The design flow is:

```text
analytical distribution specification
      |
      v
raw NDF moments or raw NDF moment tensor
      |
      v
MomentSet, CqmomMomentData, or RawMomentTensor
      |
      v
compute_qmom(...) or compute_cqmom(...)
```

After moment generation, the reconstruction algorithms have no memory of whether
the input moments came from an analytical distribution, an explicit moment list,
a sample-derived ensemble, or another upstream calculation. This keeps
distribution handling separate from quadrature reconstruction, phase-space evolution,
realizability checks, and application physics.

## Distribution Abstraction

The public distribution API lives in:

```text
include/qmomkit/distributions.hpp
```

One-dimensional analytical families are represented by strongly typed parameter
structures. The supported public variant is:

```cpp
DistributionSpecification
```

which can hold:

- `NormalDistribution`;
- `LogNormalDistribution`;
- `BetaDistribution`;
- `GammaDistribution`;
- `ChiDistribution`;
- `ChiSquaredDistribution`;
- `WeibullDistribution`.

The design intentionally uses typed C++ specifications instead of string-based
dispatch or untyped parameter maps. String names and text fields are a
configuration concern; the numerical API receives typed values.

The one-dimensional moment-generation entry point is:

```cpp
RawMomentResult compute_raw_moments(
    const DistributionSpecification& distribution,
    std::size_t maximum_order,
    double total_weight = 1.0);
```

`maximum_order` is inclusive. A request with `maximum_order = 3` returns
`m_0, m_1, m_2, m_3` on success.

Distributions also support a general selective raw-moment request:

```cpp
RawMomentRequest
SelectedRawMomentResult
```

`RawMomentRequest` stores requested phase-space-order moment multi-indices.
Request order is significant, duplicate indices are preserved, and result
values are aligned one-to-one with the request. This abstraction is generic: it
does not name QMOM or CQMOM and it is not a sparse tensor type. It allows a
caller to ask for exactly the raw moments a workflow needs.

For one-dimensional analytical distributions, each requested index has rank
one and `index[0]` is the moment order. For multidimensional independent
products, each requested index has one component per phase-space dimension.

For a single analytical distribution, the generated sequence follows the NDF
moment convention:

```text
n(xi) = n0 p(xi)
m_k = n0 E[xi^k]
m_0 = n0
```

where `p(xi)` is the normalized analytical shape and `n0` is the supplied
`total_weight`. The result is not normalized away before being returned.

## Analytical Distribution Families

Analytical families provide closed-form raw moments for normalized
one-dimensional shapes. Parameter validation happens before moment values are
used.

The current family conventions are:

- Normal: finite mean and positive finite standard deviation.
- Log-normal: finite underlying-normal mean and positive finite
  underlying-normal standard deviation.
- Beta: standard `[0, 1]` support with positive finite `alpha` and `beta`.
- Gamma: positive finite shape and scale.
- Chi: positive finite degrees of freedom.
- Chi-squared: positive finite degrees of freedom.
- Weibull: positive finite shape and scale.

Normal moments are generated through a raw-moment recurrence. Positive-support
families and beta moments use logarithmic formulas with `std::lgamma` where
needed, so large gamma or beta-function values are not formed directly before
scaling.

This page does not serve as a user catalog. The stable developer contract is
that each family maps typed parameters to raw analytical moments, validates
invalid parameters as `RawMomentStatus::invalid_input`, and reports numerical
overflow, underflow, or non-finite results through `RawMomentResult`.

## Mixture Distributions

A one-dimensional mixture is represented by:

```cpp
DistributionMixture
```

which contains `DistributionComponent` values. Each component has:

- a `DistributionSpecification`;
- an absolute component number density.

Mixture moments are accumulated as:

```text
m_k = sum_r n_r E_r[xi^k]
m_0 = sum_r n_r
```

The component values `n_r` are absolute number-density contributions. They are
not fractions, are not normalized, and are not required to sum to one. A
one-component mixture reproduces the corresponding single-distribution result
with `total_weight = n_r`.

Validation and failure behavior are part of the developer contract:

- a mixture must contain at least one component;
- at least one component must have positive number density;
- each component number density must be finite and non-negative;
- zero-density components contribute exactly zero and are skipped before their
  distribution parameters are evaluated;
- a positive-density component failure is propagated with
  `RawMomentResult::failed_component`;
- accumulation failures report the failing moment order and, where applicable,
  the component index.

Mixture accumulation uses compensated summation independently for each moment
order. This improves robustness for heterogeneous component scales, but it does
not change the mathematical convention or guarantee that all high-order moments
remain representable in `double`.

## Independent Product Distributions

Multidimensional analytical initialization is represented by:

```cpp
IndependentProductDistribution
```

It contains:

- one normalized one-dimensional marginal distribution per phase-space
  dimension;
- one global absolute number density.

The tensor moment generation entry point is:

```cpp
TensorMomentResult compute_raw_moments(
    const IndependentProductDistribution& distribution,
    const std::vector<std::size_t>& maximum_orders);
```

The product convention is:

```text
n(xi_0, ..., xi_(d-1)) = n0 product_r p_r(xi_r)
integral p_r(xi_r) dxi_r = 1
m_alpha = n0 product_r E_r[xi_r^alpha_r]
m_[0, ..., 0] = n0
```

Marginals are ordered by phase-space dimension, not by CQMOM conditioning order.
`maximum_orders[r]` is inclusive and defines
`RawMomentTensor::extents[r] = maximum_orders[r] + 1`.

The global number density is applied exactly once. Independent products do not
support per-dimension number densities, correlated distributions, copulas,
marginal mixtures, or mixtures of multidimensional products. Those concepts
would require additional public representation and validation rules rather than
being inferred from current distribution names or configuration fields.

The implementation generates each marginal moment sequence with unit total
number density, then forms tensor entries by direct multiplication with the
single global number density. Direct multiplication is intentional because
normal marginals may have negative or exactly zero odd moments; a purely
logarithmic product would not represent all supported families.

The selective independent-product implementation follows the same convention
without constructing the full tensor. It validates request rank, groups unique
requested marginal orders per phase-space dimension, evaluates each marginal
through the one-dimensional selective API, and reconstructs each requested
multidimensional moment as:

```text
m_alpha = n0 product_r E_r[xi_r^alpha_r]
```

The output preserves request order and duplicate requests.

## Numerical Behavior

Distribution moment generation uses explicit result objects:

```cpp
RawMomentResult
TensorMomentResult
RawMomentStatus
```

Status values are:

- `success`;
- `invalid_input`;
- `overflow`;
- `underflow`;
- `non_finite`.

Callers must check `ok()` or `status` before using generated moments. Non-empty
result data is not a success criterion.

Stable numerical rules include:

- total weights and global number densities must be finite and non-negative;
- zero total weight for a single distribution returns a zero moment sequence;
- zero global number density for an independent product returns a zero tensor
  with the requested valid extents and does not evaluate marginal parameters;
- positive-support moments are checked before exponentiation where logarithmic
  formulas are used;
- exactly positive moments that underflow to zero are reported as `underflow`;
- infinite results are reported as `overflow`;
- other non-finite intermediate or final values are reported as `non_finite`;
- `failed_order`, `failed_component`, `failed_dimension`, and `failed_order`
  diagnostics identify the location of failures when available.

Distribution generation validates parameters and numerical representability. It
does not perform support-aware realizability analysis, run QMOM or CQMOM
inversion, or repair generated moments. Downstream algorithms still report their
own statuses.

## Initialization Interaction

Distribution-generated moments enter QMOMKit through the same public case
objects as explicit moments. One-dimensional QMOM initialization uses
`MomentSet`:

```text
DistributionSpecification
      -> RawMomentResult
      -> MomentSet
      -> QmomCase

DistributionMixture
      -> RawMomentResult
      -> MomentSet
      -> QmomCase
```

Normal configured CQMOM initialization for independent products uses selected
required moment data:

```text
IndependentProductDistribution
      -> SelectedRawMomentResult
      -> CqmomMomentData
      -> CqmomCase
```

The explicit dense CQMOM initialization route remains available when callers
intentionally use full rectangular tensor storage:

```text
IndependentProductDistribution
      -> TensorMomentResult
      -> RawMomentTensor
      -> CqmomDenseCase
```

For CQMOM workflows that can benefit from selective generation, the data flow is:

```text
CQMOM hierarchy
      -> CqmomMomentLayout
      -> RawMomentRequest

IndependentProductDistribution
      -> SelectedRawMomentResult
      -> CqmomMomentData
      -> CqmomCase

or, for explicit lean representation workflows:

SelectedRawMomentResult
      -> CqmomMomentData
      -> CqmomLeanCase
```

The CQMOM layer defines `CqmomMomentLayout`, CQMOM moment data, and
reconstruction. The distribution layer only sees the generic `RawMomentRequest`
and returns generic selected raw moments. Conversion from selected raw moments
into CQMOM moment data occurs outside the distribution layer, so distributions
remain independent of QMOM/CQMOM reconstruction details.

The dense `CqmomDenseCase` workflow remains supported. Public CQMOM
reconstruction internally extracts the required CQMOM moment data from a dense
input tensor before recursive reconstruction.

Configuration parsing uses the same public generation functions. Explicit
moment initialization has precedence in configuration: explicit moments are
treated as the source of truth, and distribution fields that would otherwise
generate moments are not active in that path.

Generated moments must still be sufficient for the requested reconstruction. For
one-dimensional QMOM, an `N`-node inversion needs at least `2 * N` moments. For
dense CQMOM compatibility, the generated tensor extents must be large enough for
the configured phase-space hierarchy and node counts. For lean CQMOM
initialization, sufficiency is encoded by `CqmomMomentLayout`. Those
reconstruction requirements belong to QMOM case and CQMOM case validation and
algorithm status handling, not to the distribution object itself.

## API Boundary

Public distribution concepts include:

- `NormalDistribution`;
- `LogNormalDistribution`;
- `BetaDistribution`;
- `GammaDistribution`;
- `ChiDistribution`;
- `ChiSquaredDistribution`;
- `WeibullDistribution`;
- `DistributionSpecification`;
- `DistributionComponent`;
- `DistributionMixture`;
- `IndependentProductDistribution`;
- `RawMomentResult`;
- `RawMomentRequest`;
- `SelectedRawMomentResult`;
- `TensorMomentResult`;
- `RawMomentStatus`;
- `compute_raw_moments(...)`.

Implementation details include:

- family-specific moment kernels;
- recurrence and logarithmic evaluation helpers;
- compensated summation machinery;
- checked tensor-product assembly;
- internal failure-construction helpers.

These implementation details support the public API but are not downstream
extension points.

## Extension Principles

Adding a distribution capability should preserve the separation between
initialization and reconstruction:

- add a typed public representation for new supported analytical families;
- define parameter validation and moment formulas explicitly;
- preserve raw NDF moment conventions and `m_0`/`m_[0,...,0]` semantics;
- define finite-value, overflow, underflow, and diagnostic behavior;
- keep configuration string mapping separate from the numerical API;
- add unit and verification coverage with independent analytical references;
- avoid inferring physical units, supports, or realizability policies from
  distribution names.

Distribution support should remain a deterministic moment-generation layer. It
should not become an application-level solver, a sampling library, a PDF/CDF evaluation
package, or a hidden control path for QMOM/CQMOM reconstruction behavior.
