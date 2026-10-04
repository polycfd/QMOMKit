# Distribution-based moment initialization

QMOMKit can generate raw moments from supported analytical distributions. This
is useful when an initial population is easier to describe through familiar
distribution parameters than through a manually assembled moment table. QMOM
and CQMOM still consume moments; the distribution helpers simply create those
moments with QMOMKit's raw number-density convention.

The workflow is:

```text
distribution
      ↓
raw moments
      ↓
MomentSet / required CQMOM moment data / RawMomentTensor
      ↓
QMOM/CQMOM reconstruction
```

QMOMKit does not carry a probability-density object through the rest of the
workflow. Distribution helpers are initialization tools: they produce raw
moments, and the usual QMOM or CQMOM APIs consume those moments. A full tensor
of raw moments is one explicit CQMOM input; configured CQMOM initialization
generates only the moments the hierarchy requires.

## Raw moment convention

Generated moments are unnormalized raw moments of the number density function.
The distribution parameters describe the shape of the population, and
`number_density` supplies the total amount of distribution mass carried by the
moment sequence.
For a one-dimensional distribution:

```text
m_k = number_density * E[xi^k]
m_0 = number_density
```

`number_density` is an absolute total weight. It is not a probability, fraction,
or normalized mixture coefficient unless you intentionally choose values that
make it so.

For a multidimensional independent product:

```text
M[k_0, ..., k_D-1] = number_density * product_d E[x_d^k_d]
M[0, ..., 0] = number_density
```

## Supported distributions

QMOMKit provides analytical raw moments for these one-dimensional distribution
families:

- normal, parameterized by mean and standard deviation;
- log-normal, parameterized by the underlying normal mean and standard
  deviation;
- beta on `[0, 1]`, parameterized by alpha and beta;
- gamma, parameterized by shape and scale;
- chi, parameterized by degrees of freedom;
- chi-squared, parameterized by degrees of freedom;
- Weibull, parameterized by shape and scale.

All parameters must be finite. Normal and log-normal location parameters may
be any finite value, but their standard deviations must be positive. Beta
requires positive `alpha` and `beta` and has support `[0, 1]`; gamma and
Weibull require positive shape and scale; chi and chi-squared require positive
degrees of freedom. In `.qmom` files, use the canonical family names
`lognormal` and `chi-squared`; `log-normal`, `log_normal`, `chi_squared`, and
`chisquared` are also accepted aliases.

The C++ API uses typed distribution specifications, for example
`NormalDistribution`, `LogNormalDistribution`, and `GammaDistribution`. The
configuration-file guide uses the corresponding field names accepted by the
parser.

Moment generation reports a status. Check `ok()` before passing generated
moments into QMOM or CQMOM so invalid parameters or unsupported requests are
caught before reconstruction.

## QMOM initialization

For one-dimensional QMOM, generate moments through `compute_raw_moments`, then
pass the returned `MomentSet` to `compute_qmom`. This mirrors the physical
workflow: choose an initial distribution, convert it to raw moments, then
reconstruct the quadrature representation your application will inspect or
evolve.

```cpp
#include <qmomkit/qmomkit.hpp>

#include <cstddef>
#include <iostream>

int main()
{
    const std::size_t node_count = 2;
    const std::size_t maximum_order = 2 * node_count - 1;
    const double number_density = 2.5;

    const qmomkit::DistributionSpecification distribution =
        qmomkit::LogNormalDistribution{-0.2, 0.4};

    const qmomkit::RawMomentResult moments =
        qmomkit::compute_raw_moments(distribution, maximum_order, number_density);
    if (!moments.ok())
    {
        std::cerr << "Moment generation failed: " << moments.message << '\n';
        return 1;
    }

    const qmomkit::QmomResult qmom =
        qmomkit::compute_qmom(moments.moments, node_count);
    if (!qmom.ok())
    {
        std::cerr << "QMOM failed: " << qmom.message << '\n';
        return 1;
    }

    std::cout << "m0 = " << moments.moments[0] << '\n';
    for (std::size_t node = 0; node < qmom.weights.size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << qmom.weights[node]
                  << ", abscissa=" << qmom.abscissas[node] << '\n';
    }
}
```

`maximum_order` is inclusive. A `maximum_order` of `3` generates `m0` through
`m3`, which is enough for a two-node QMOM inversion.

## Mixtures

Mixtures are useful when one analytical distribution is not enough to describe
a population. For example, a size distribution might have two peaks, or a
sample population might combine multiple regimes.

QMOMKit mixtures use absolute component densities:

```text
m_k = sum_r number_density_r * E_r[xi^k]
m_0 = sum_r number_density_r
```

QMOMKit does not normalize component densities.
Component densities must be finite and non-negative, and a mixture must contain
at least one component with positive density.

```cpp
#include <qmomkit/qmomkit.hpp>

#include <cstddef>
#include <iostream>

int main()
{
    const std::size_t node_count = 3;
    const std::size_t maximum_order = 2 * node_count - 1;

    const qmomkit::DistributionMixture mixture{{
        {qmomkit::GammaDistribution{2.0, 0.4}, 1.5},
        {qmomkit::LogNormalDistribution{-0.2, 0.5}, 0.8},
    }};

    const qmomkit::RawMomentResult moments =
        qmomkit::compute_raw_moments(mixture, maximum_order);
    if (!moments.ok())
    {
        std::cerr << "Mixture moment generation failed: " << moments.message << '\n';
        return 1;
    }

    std::cout << "mixture m0 = " << moments.moments[0] << '\n';
}
```

The generated `MomentSet` can be passed directly to `compute_qmom`.

## CQMOM initialization

For multidimensional CQMOM, independent-product initialization covers cases
where the initial coordinates can be described as independent marginal
distributions. Each phase-space dimension gets one simple distribution or
mixture marginal, normalized through its zeroth raw moment, and the product has
one global `number_density`. Mixture component densities retain their usual
absolute-density meaning.

This is not a general joint-distribution interface: correlated distributions,
copulas, multidimensional mixture components, and mixtures of independent
products are not represented. See the configuration guide for the nested
`.qmom` marginal-mixture syntax.

For the normal configured workflow, describe those marginals in a
multidimensional `.qmom` configuration file and load it with
`load_cqmom_prepared(...)`. QMOMKit then generates the required CQMOM moments
for the requested hierarchy. Users do not need to build or manage that required
moment set themselves.

The direct C++ example below takes the explicit full-tensor path, for
applications that already manage rectangular tensor moments.

```cpp
#include <qmomkit/qmomkit.hpp>

#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
    const qmomkit::IndependentProductDistribution distribution{
        {
            qmomkit::GammaDistribution{3.0, 0.2},
            qmomkit::NormalDistribution{-0.5, 1.25},
        },
        2.5};

    const std::vector<std::size_t> maximum_orders{3, 3};
    const qmomkit::TensorMomentResult tensor =
        qmomkit::compute_raw_moments(distribution, maximum_orders);
    if (!tensor.ok())
    {
        std::cerr << "Tensor moment generation failed: " << tensor.message << '\n';
        return 1;
    }

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {0, 1});
    const qmomkit::CqmomCase cqmom_case{
        qmomkit::extract_cqmom_moments(tensor.moments, qmomkit::create_cqmom_moment_layout(hierarchy))};
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
    if (!cqmom.ok())
    {
        std::cerr << "CQMOM failed: " << cqmom.message << '\n';
        return 1;
    }

    std::cout << "M[0,0] = " << tensor.moments.values.front() << '\n';
}
```

`maximum_orders` are inclusive. In the example above, the generated tensor has
extents `{4, 4}`.
Choose tensor extents that satisfy the requested CQMOM node counts.
`extract_cqmom_moments(...)` then selects the moments the hierarchy
reconstructs from, and throws if the tensor is too small to contain one.

Keep these settings separate:

- Marginal order follows phase-space dimension order.
- `maximum_orders` controls how many raw tensor moments are generated for a
  full tensor.
- `node_counts` controls the requested CQMOM quadrature size.
- `conditioning_order` controls the CQMOM conditioning order. It does not
  change tensor storage order or select different marginals.

## Configuration files (`.qmom`)

The same distribution workflows are available through QMOMKit configuration
files. For one-dimensional QMOM, use:

```ini
[initial]
type = distribution
family = lognormal
number_density = 2.5
maximum_order = 5
mu_log = -0.2
sigma_log = 0.4
```

For multidimensional CQMOM independent products, use one section per
phase-space dimension. The normal CQMOM configuration path derives the required
raw moments from the configured CQMOM hierarchy:

```ini
[initial]
type = independent_product
number_density = 2.5

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = 0.0
standard_deviation = 1.0
```

See [Configuration](configuration.md) for complete configuration-file
examples and field semantics.

## Practical considerations

Keep the raw-moment convention consistent. Passing moments with their physical
$m_0$ preserves the total number density; choosing $m_0 = 1$ represents a
unit-density setup because the caller chose that normalization.

Choose `number_density` with the same physical meaning as the total weight in
your model. For mixtures, component densities are absolute contributions to
the total $m_0$.

Keep parameter units consistent with the phase-space coordinate. If a gamma
distribution represents a size coordinate, its scale has the same unit as that
coordinate, and generated moments carry the corresponding powers of that unit.

Generate enough moment orders for the requested quadrature. For $N$ QMOM nodes,
generate at least $2N - 1$ as the inclusive maximum order. For normal CQMOM
configuration, `node_counts` and `conditioning_order` define the selected
moment request. For explicit full-tensor CQMOM workflows, choose
`maximum_orders[d]` large enough for the node count in each phase-space
dimension.

Check both moment generation and quadrature reconstruction statuses. A valid
analytical distribution can still lead to a quadrature request that fails if
too few moments were generated or if the configured reconstruction is not
accepted.

Realizability checks are optional and separate. They are useful when moments
come from external data, approximate calculations, or runtime updates. Moments
generated directly from supported analytical distributions should satisfy the
distribution's own moment convention, but downstream transformations can still
require explicit status handling.
