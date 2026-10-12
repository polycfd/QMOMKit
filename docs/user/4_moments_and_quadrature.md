# Moments and quadrature

This page covers the first half of every workflow: setting up a case,
reconstructing its quadrature, and checking whether moments are realizable.
The examples continue those of [Getting started](2_getting_started.md).

## Set up a case

A case holds the initial moments, the node counts, and the reconstruction
options. Its moments can come from your own data, from distributions, or from a
configuration file.

Choose units in which the coordinates are of order one, for example droplet
volumes in cubic micrometres rather than cubic metres. High powers of very large
or very small values do not fit in a double.

### From your own moments

A QMOM case takes the moments in ascending order, $m_0, m_1, \ldots$:

```cpp
qmomkit::QmomCase qmom_case;
qmom_case.moments = qmomkit::MomentSet({1.0, 0.8, 2.8, 4.4});
qmom_case.node_count = 2;
```

A CQMOM case needs only the moments its node counts and conditioning order
require. `create_cqmom_moment_layout(hierarchy)` lists them in
`layout.required_indices`; supply one value per entry, in that order:

```cpp
const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
const qmomkit::CqmomCase cqmom_case{qmomkit::CqmomMomentData{layout, values}};
const qmomkit::CqmomPreparedCase prepared_case{cqmom_case};
```

If you hold a full tensor of moments instead, `extract_cqmom_moments(tensor,
layout)` selects the required ones and throws if one is missing. Prepare a
CQMOM case once, as `CqmomPreparedCase`, and reuse it.

### From distributions

`compute_raw_moments(...)` generates moments from a distribution, scaled by a
number density. For two nodes, generate $m_0$ to $m_3$, so the maximum order is
$2N - 1 = 3$:

```cpp
const qmomkit::RawMomentResult lognormal_moments = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{-0.2, 0.4}, 3, 2.5);
if (!lognormal_moments.ok())
{
    std::cerr << lognormal_moments.message << '\n';
    return 1;
}
qmom_case.moments = lognormal_moments.moments;
```

A mixture adds the moments of its components. Each component carries its own
number density, so this mixture has $m_0 = 1.5 + 0.8 = 2.3$:

```cpp
const qmomkit::DistributionMixture mixture{{
    {qmomkit::GammaDistribution{2.0, 0.4}, 1.5},
    {qmomkit::LogNormalDistribution{-0.2, 0.5}, 0.8},
}};
const qmomkit::RawMomentResult generated = qmomkit::compute_raw_moments(mixture, 5);
```

For CQMOM, an independent product combines one distribution per coordinate
with one total number density. Ask for exactly the moments the layout requires:

```cpp
const qmomkit::IndependentProductDistribution distribution{
    {qmomkit::GammaDistribution{3.0, 0.2}, qmomkit::NormalDistribution{-0.5, 1.25}}, 2.5};

const qmomkit::SelectedRawMomentResult generated =
    qmomkit::compute_raw_moments(distribution, qmomkit::RawMomentRequest{layout.required_indices.to_multi_indices()});
if (!generated.ok())
{
    std::cerr << generated.message << '\n';
    return 1;
}
const qmomkit::CqmomCase cqmom_case{qmomkit::create_cqmom_moment_data(layout, generated)};
```

A coordinate may also follow a mixture. Correlated distributions cannot be
expressed this way; supply their moments yourself. The
[Reference](7_reference.md#distribution-families) lists the families and
their parameters.

### From a configuration file

`load_qmom(...)` and `load_cqmom_prepared(...)` read a case from a `.qmom` file,
as in [Getting started](2_getting_started.md). The file gives either explicit
moments in `[moments]` or distributions in `[initial]`; the
[Reference](7_reference.md#configuration-file) lists every key.
`load_cqmom(...)` returns an editable `CqmomCase` instead, to change in code
before preparing it.

For example, this file, `qmom.qmom`, holds the one-dimensional case from
[Getting started](2_getting_started.md) and, for its evolution, the splitting
scheme:

```ini
[phase_space]
node_counts = 2

[moments]
values = 1.0, 0.8, 2.8, 4.4

[time_integration]
scheme = strang
```

It is loaded with:

```cpp
const auto configured = qmomkit::config::load_qmom("qmom.qmom");
if (!configured)
{
    std::cerr << configured.error() << '\n';
    return 1;
}
const qmomkit::QmomCase qmom_case = configured.value();
```

The following pages use this file.

For CQMOM, to keep a loaded case's settings but replace its moments, for
example with moments from another code:

```cpp
const auto replaced = qmomkit::config::with_initial_moments(prepared_case, moments);
if (!replaced)
{
    std::cerr << replaced.error() << '\n';
    return 1;
}
const qmomkit::CqmomPreparedCase new_case = replaced.value();
```

`moments` may be a full tensor or the required moment data. For QMOM, assign
the new moments to `qmom_case.moments`.

## Reconstruct

`compute_qmom(...)` returns the weights and abscissas, in ascending order of
abscissa:

```cpp
const qmomkit::QmomResult result = qmomkit::compute_qmom(qmom_case);
```

`compute_cqmom(...)` returns a `CqmomQuadrature` whose nodes you iterate:

```cpp
const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(prepared_case);
if (cqmom.ok())
{
    const qmomkit::CqmomQuadrature& quadrature = cqmom.quadrature;
    for (std::size_t i = 0; i < quadrature.total_nodes(); ++i)
    {
        std::cout << quadrature.weights()[i] << " at (" << quadrature.abscissas()[0][i] << ", " << quadrature.abscissas()[1][i] << ")\n";
    }
}
```

### Adaptive or fixed node counts

Reconstruction is adaptive by default, as
[Concepts](3_concepts.md#adaptive-and-fixed-node-counts) describes. Request a
fixed one with:

```cpp
qmom_case.qmom_options.adaptive = false;
```

For CQMOM the setting applies to every branch. Only a fixed CQMOM quadrature
lets you address a node by its multi-index, with `weight(index)` and
`abscissa(dimension, index)`.

Choose node counts by the resolution your model needs. More nodes need more,
and higher, moments; fewer nodes are cheaper and more robust when only the low
moments are reliable. In CQMOM the number of nodes is the product of the node
counts of the coordinates, and the number of required moments grows with it, so
each added coordinate or node makes every step more expensive.

### Bounds on the nodes

When a coordinate cannot leave an interval, require the nodes to lie in it:

```cpp
qmom_case.abscissa_bounds = qmomkit::AbscissaBounds{0.0, 10.0};
```

Nodes outside the interval are never moved into it. An adaptive reconstruction
then retries with fewer nodes; a fixed one fails. The bounds hold for every
later reconstruction of the states made from the case; node advection between
reconstructions does not enforce them. Bounds are never inferred from a
distribution or a coordinate name. For CQMOM, set one optional interval per
coordinate on the `CqmomCase` before preparing it, for example
`abscissa_bounds = {std::nullopt, qmomkit::AbscissaBounds{0.0, 1.0}}`. In a
file, `[phase_space] support_lower` and `support_upper` set the bounds; they
are unrelated to the `support` of a realizability check.

### Conditioning order

The second argument of `PhaseSpaceHierarchy` sets the CQMOM conditioning order;
here the second coordinate is reconstructed first:

```cpp
const qmomkit::PhaseSpaceHierarchy reordered({2, 3}, {1, 0});
```

The node counts still belong to the coordinates in their declared order.
Different orders give different quadratures of the same moments, and can
differ in numerical robustness.

### Solver and diagnostics

CQMOM solves a small linear system for each conditional reconstruction. The
default solver suits most cases; `[vandermonde] solver` selects another, and
`[cqmom] diagnostic_level` records more detail about each branch. See the
[Reference](7_reference.md#cqmom-and-vandermonde), and the
[developer documentation](../developer/linear_algebra.md) for the solvers.

### Derived quantities

Moments of any real order can be evaluated from a quadrature, here $m_{1/3}$
of the log-normal moments generated above and, in two dimensions,
$M_{1/3,\,2}$ of the CQMOM result:

```cpp
const qmomkit::QmomResult lognormal_quadrature = qmomkit::compute_qmom(lognormal_moments.moments, 2);
const double m_third = qmomkit::evaluate_generalized_moment(lognormal_quadrature, 1.0 / 3.0);
const double cqmom_m = qmomkit::evaluate_generalized_moment(cqmom.quadrature, {1.0 / 3.0, 2.0});
```

Fractional orders need non-negative abscissas, as here, and throw
`std::domain_error` otherwise. These values are computed from the quadrature
only; they are not added to the moments a case or state holds.

## Check realizability

Not every list of numbers is the moments of a non-negative number density.
Whether it is depends on the *support*, the range of the coordinate.
`check_realizability(...)` decides this for one coordinate:

```cpp
const qmomkit::RealizabilityResult check =
    qmomkit::check_realizability(qmomkit::MomentSet({1.0, 0.8, 2.8, 4.4}), qmomkit::MomentSupport::nonnegative);
if (!check.realizable())
{
    std::cerr << check.message << '\n';
}
```

```text
Moment sequence is not realizable: S_1 has a negative eigenvalue beyond tolerance. Support: nonnegative.
```

These moments have a node at $-1$, so they are realizable on the real line but
not on the non-negative numbers. The supports are `MomentSupport::real_line`,
`MomentSupport::nonnegative`, and an interval, `BoundedInterval{lower, upper}`.
`realizable()` is true for moments strictly inside the realizable set and on
its boundary, which `on_boundary()` reports. The message names the failing
moment matrix.

A successful reconstruction is no realizability check: QMOM works on the real
line, and an adaptive reconstruction drops nodes the moments do not support.
To keep nodes inside an interval, use bounds as above.

For a distribution far from the origin compared with its width, at high
orders, the moments carry too few digits to decide, and the check reports
`boundary` or `numerically_indeterminate`. Measuring the coordinate from a
nearby reference point avoids this.

### During evolution

A realizability policy checks the moments after each source step of an
evolution:

```cpp
qmomkit::RealizabilityPolicy policy;
policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
policy.support = qmomkit::MomentSupport::nonnegative;

qmomkit::QmomPhaseSpaceEvolutionOptions qmom_options;
qmom_options.realizability_policy = policy;

qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_options;
cqmom_options.source_options.coordinate_realizability_policies = {policy, std::nullopt};
```

`FinalOnly` checks the result of the source step, `EveryStage` also its
intermediate stages. Boundary moments are accepted unless `accept_boundary` is
false, indeterminate ones only if `accept_indeterminate` is set; accepted, they
may lie just outside the realizable set, and an adaptive reconstruction then
drops nodes. A rejected source step fails, or is retried with a smaller step if
you allow retries, as described in [Evolving in time](5_evolving_in_time.md).
For CQMOM, each coordinate's policy checks only the moments of that coordinate,
$M_{0,\ldots,k,\ldots,0}$. In a file, these settings form the
`[time_integration.realizability]` section.
