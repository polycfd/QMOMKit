# CQMOM

CQMOM, the Conditional Quadrature Method of Moments [[Yuan & Fox,
2011](https://www.doi.org/10.1016/j.jcp.2011.07.020)], extends quadrature-based
moment reconstruction to multidimensional phase spaces. A multidimensional
population may depend on variables such as size, velocity, temperature, or any
other application-defined properties. QMOMKit treats these generically as
phase-space dimensions.

Use CQMOM when the population is described by a joint distribution over several
phase-space variables and the dependence between those variables matters. The
reconstructed quadrature consists of terminal nodes, each with one global
weight and one abscissa in every configured phase-space dimension. This
provides a compact representation that can be used to evaluate multidimensional
closures or evolve coupled phase-space dynamics.

In one dimension, QMOM reconstructs a quadrature from a raw moment sequence. In
multiple dimensions, the CQMOM implementation of QMOMKit reconstructs the
quadrature recursively [[Bouguezzoul,
2026](https://publications.polymtl.ca/76505/)] through marginal and conditional
one-dimensional reconstructions.

A CQMOM case defines the phase-space dimensions, maximum node counts,
conditioning order, initial moments, and reconstruction options used
consistently throughout the calculation.

QMOMKit uses raw, unnormalized tensor moments:
$$M[k_0, ..., k_{D-1}] = \sum_j W_j \prod_d xi_{d,j}^{k_d}$$
The zeroth tensor moment `M[0, ..., 0]` is the total quadrature weight or
number density.

## CQMOM workflow

The typical CQMOM workflow starts from a modelling description rather than from
internal storage. That description can come from a `.qmom` configuration file
or from direct C++ setup. In either case, users name or order the coordinates,
set the requested node counts, choose the conditioning order, and describe the
initial moments. File-based setup builds a `CqmomPreparedCase`; direct C++
setup can build a `CqmomLeanCase` when the application assembles the required
moments itself, then promote it to the normal `CqmomPreparedCase` workflow.

Both case types hold the same moments and produce the same quadrature. A
prepared case also carries the recursive branch structure implied by the
phase-space hierarchy, so reconstructing repeatedly from the same hierarchy
does not rebuild it. Prefer the prepared case unless you have a reason to hold
the moments yourself.

```text
configuration file:
    .qmom file
    -> load_cqmom_prepared(...)
    -> CqmomPreparedCase

direct API:
    C++ setup
    -> CqmomLeanCase
    -> CqmomPreparedCase

reconstruction:
    CqmomPreparedCase or CqmomLeanCase
    -> compute_cqmom(...)
    -> CqmomResult / CqmomQuadrature
    -> inspect / evaluate / output

normal optional phase-space evolution:
    CqmomPreparedCase::initial_state()
    -> CqmomState
    -> advance_cqmom_phase_space(...)
```

At the user level, this means:

1. Define the CQMOM setup in a `.qmom` file or directly in C++.
2. Load it with `load_cqmom_prepared(...)`, or construct a `CqmomLeanCase`
   and pass it to `CqmomPreparedCase` in C++.
3. Keep the case object as the configured problem definition.
4. Call `compute_cqmom(...)` to compute the initial quadrature.
5. Inspect the quadrature, evaluate derived quantities, or write output.
6. If the case evolves in time, call `CqmomPreparedCase::initial_state()` to
   create the runtime state for phase-space evolution and pass that
   `CqmomState` to `advance_cqmom_phase_space(...)`.

The result is a set of weighted multidimensional nodes. Each terminal node has
one global weight and one abscissa per phase-space dimension.

## Automatic reconstruction scaling

Set `QmomOptions::automatic_scaling = true` on the shared CQMOM/QMOM options,
or set `automatic_scaling = true` in the configuration `[phase_space]`
section. The default is `false`. CQMOM selects one scale for each physical
coordinate from the complete current required moment set, transforms the
recursive problem, then returns nodes and bounded-support diagnostics in
physical coordinates. Raw moments, weights, callbacks, states, and configured
coordinate supports keep their physical meaning.

`CqmomPreparedCase` is the immutable setup produced by the normal configuration
loader or by promoting a `CqmomLeanCase`. Normal CQMOM phase-space evolution
creates `CqmomState` from `CqmomPreparedCase` and advances it with
`advance_cqmom_phase_space(...)`. `CqmomLeanCase` also remains available for
static inversion with `compute_cqmom(...)`. In all cases, `compute_cqmom(...)`
returns quadrature; it does not return a phase-space evolution state.

For the usual CQMOM path, users specify the phase space, initial condition, and
reconstruction settings. QMOMKit determines which moments are needed for the
configured hierarchy and avoids exposing the internal moment storage.

### Replacing initial moments programmatically

When an application obtains a complete `RawMomentTensor` after loading its
ordinary CQMOM setup, it can make a separate phase-space evolution-capable case
without rewriting configuration text:

```cpp
const auto replaced = qmomkit::config::with_initial_moments(base_case, std::move(raw_moments));
if (!replaced)
{
    std::cerr << replaced.error() << '\n';
    return 1;
}
const qmomkit::CqmomPreparedCase cqmom_case = replaced.value();
```

The new case retains the base hierarchy, conditioning order, support, and every
reconstruction policy, and selects the moments the hierarchy needs from the
supplied tensor. The transformation
checks rank, extents, value count, finiteness, and required coverage only. It
does not test realizability: `compute_cqmom(...)` reports numerical
reconstruction status, and `initial_state()` keeps its existing
exception-on-failure behavior. The producer of the raw moments can be a
sampler, measurement pipeline, or any other application code.

When the producer can supply the required CQMOM moment set directly, pass
`CqmomMomentData` instead of a tensor:

```cpp
const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(base_case.hierarchy());
const qmomkit::RawMomentRequest request{layout.required_indices.to_multi_indices()};
const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(product_distribution, request);
const auto replaced = qmomkit::config::with_initial_moments(base_case, qmomkit::create_cqmom_moment_data(layout, selected));
```

This asks for exactly the moments the hierarchy needs rather than a full
rectangular tensor, which matters when the two differ substantially: a
five-dimensional hierarchy with three nodes per coordinate requires 606 moments,
while the full tensor holds 7776. The resulting case is the same either way.

## Using CQMOM in QMOMKit

Include the umbrella header when following the user guide examples:

```cpp
#include <qmomkit/qmomkit.hpp>
```

The example below loads a two-dimensional CQMOM case from configuration and
computes the initial quadrature. It stops at reconstruction; phase-space
evolution is shown afterward as an optional step.

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    const auto configured = qmomkit::config::load_cqmom_prepared("configuration.qmom");
    if (!configured)
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    const qmomkit::CqmomPreparedCase cqmom_case = configured.value();

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);
    if (!result.ok())
    {
        std::cerr << "CQMOM failed: " << result.message << '\n';
        return 1;
    }

    const qmomkit::CqmomQuadrature& quadrature = result.quadrature;

    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << quadrature.weights()[node];
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            std::cout << ", x" << dimension << "=" << quadrature.abscissas()[dimension][node];
        }
        std::cout << '\n';
    }
}
```

If this reconstructed distribution evolves in time, create the phase-space
evolution state from the same case and pass it to
`advance_cqmom_phase_space(...)`:

```cpp
qmomkit::CqmomState state = cqmom_case.initial_state();
const qmomkit::CqmomQuadrature& quadrature = state.quadrature();
```

Write the moments represented by the current state with a state-aware moment
writer. A state writes the moments its hierarchy reconstructs from:

```cpp
qmomkit::OutputOptions output_options;
qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom("moments.txt", state, output_options);
writer.write(time, state, quadrature);
```

For a simple independent-product setup, the configuration file can look like
this. The example describes a two-coordinate population with one global number
density and one analytical marginal distribution per phase-space dimension:

```ini
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 0, 1

[initial]
type = independent_product
number_density = 2.5

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
```

### Choosing the Vandermonde solver

CQMOM reconstructs conditional quantities by solving branch-local
conditional Vandermonde systems. QMOMKit provides several solver policies for
those systems through the optional `[vandermonde]` section:

```ini
[vandermonde]
solver = automatic
```

This line may be omitted when using the default. Supported values are:

| Value | Practical meaning |
|---|---|
| `automatic` | Default policy; uses Björck-Pereyra. |
| `dense` | Uses dense Gaussian elimination. |
| `pivoted_qr` | Uses pivoted QR with numerical-rank detection. |
| `svd` | Uses SVD with rank and convergence checks; rank-deficient systems are rejected. |
| `moore_penrose` | Uses a truncated SVD Moore-Penrose pseudoinverse. |
| `bjorck_pereyra` | Uses Björck-Pereyra. |

QR and SVD report rank-deficient systems; SVD also reports nonconvergence.
`moore_penrose` uses a truncated SVD with cutoff controlled by `rcond`; see the
[Configuration Guide](configuration.md) for its detailed `rcond` contract and
exact policy behavior.
Developer-level numerical details are in the
[Linear Algebra Guide](../developer/linear_algebra.md).

### Reconstruction diagnostics

CQMOM always runs the numerical checks required to accept or reject a branch,
including solver acceptance, support validation, and
adaptive branch selection. The optional `[cqmom] diagnostic_level` controls
only retained observability data and defaults to `none`:

```ini
[cqmom]
diagnostic_level = basic
```

Use `basic` for structural branch records and `detailed` when residual and
condition-estimate diagnostics are needed; detailed reporting may perform extra
work. A failed reconstruction still returns its status and available failure
context at every level. See the [Configuration Guide](configuration.md) for the
exact syntax.

For this normal configuration path, the CQMOM hierarchy determines the
required raw moments automatically. `initial.maximum_orders` is therefore not a
reconstruction requirement and does not affect a configured case. It applies to
applications that build a full raw moment tensor themselves, where larger values
retain a larger rectangular initial raw-moment tensor. They do not change the
hierarchy, reconstruction order, or branch counts. Sample-based initialization
can use the larger tensor extent when evaluating raw moments from samples.

For fixed rectangular quadratures, you can also access a node by its
phase-space-order multi-index:

```cpp
const qmomkit::MultiIndex node_index{0, 1};
const double weight = result.quadrature.weight(node_index);
const double x0 = result.quadrature.abscissa(0, node_index);
const double x1 = result.quadrature.abscissa(1, node_index);
```

## Conditioning and ordering

Phase-space dimensions are the stable coordinate order used by moment tensors,
node counts, output coordinates, and callbacks. Choose this order to match the
rest of your application and keep it unchanged across configuration, physics
callbacks, and validation data. For a two-dimensional problem:

```text
dimension 0: first coordinate
dimension 1: second coordinate
```

`PhaseSpaceHierarchy` stores node counts by phase-space dimension:

```cpp
const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
```

This requests two nodes in dimension `0` and three nodes in dimension `1`.

The conditioning order controls the numerical order in which dimensions are
processed. It is a numerical reconstruction choice, not a renaming of
coordinates. The default order is identity order:

```text
0, 1, 2, ...
```

To specify a different order, pass a second vector:

```cpp
const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3}, {1, 0});
```

This still means dimension `0` has two nodes and dimension `1` has three
nodes. It only changes the conditioning order. It does not change tensor
storage order, coordinate order, or the meaning of the dimensions.

Ordering matters because different conditioning orders may expose different
numerical conditioning in a multidimensional reconstruction. Keep the order
explicit in configuration files, examples, and validation output so results
are reproducible and callback coordinate order stays easy to audit.

## Fixed and adaptive quadrature

The default CQMOM setup uses a fixed branch quadrature. This is the easiest
mode to reason about: each branch must produce the requested number of nodes
for its phase-space dimension, and fixed results support rectangular node
access by multi-index.

Adaptive branch quadrature is available when a workflow can tolerate fewer
active nodes on branches where the supplied moments support only a lower-order
rule. This can be useful near degenerate or boundary-like states where a fixed
node count may be too ambitious. Enable it explicitly in configuration with
adaptive branch counts:

```ini
[phase_space]
adaptive_branch_counts = true
```

Adaptive CQMOM branch counts imply adaptive Wheeler internally for all
marginal and conditional inversions. The standalone QMOM adaptive option
remains independent.

Adaptive results expose terminal weights, terminal abscissas, and terminal
conditioning paths. When adaptive branch counts are enabled, prefer iterating
over `weights()` and `abscissas()` instead of assuming rectangular multi-index
storage.

Adaptive branch counts are supported in both reconstruction and the normal
configured `CqmomState` phase-space evolution workflow. `CqmomState` keeps the
required moment representation for the configured hierarchy, while each
reconstructed quadrature may contain adaptive active node counts.

As with QMOM, always check `result.ok()` before using the quadrature. If CQMOM
fails, inspect `result.status`, `result.message`, and diagnostics to determine
which phase-space dimension or conditioning step was not accepted.

## Moment reconstruction from quadrature

Native CQMOM moments remain integer-order polynomial raw moments with
phase-space-order multi-indices. QMOMKit reconstructs those moments from
terminal quadrature when output, validation, or phase-space evolution handoff
needs raw moment values.

After checking that a CQMOM quadrature result is valid, callers may also
evaluate generalized moments with finite real exponent vectors as derived
post-processing quantities:

```cpp
const std::vector<double> exponents{1.0 / 3.0, 2.0};
const double derived =
    qmomkit::evaluate_generalized_moment(result.quadrature, exponents);
```

Generalized moment evaluation applies the terminal quadrature formula directly.
These values are not CQMOM reconstruction inputs, phase-space-evolved state, or native
`CqmomMomentData`/`RawMomentTensor` entries. Non-integer exponents require
non-negative abscissas in the corresponding dimensions.

## Advanced representation workflows

The workflows in this section are for static inversion, specialized
initialization, validation, and expert use. They are for applications that
assemble the moments themselves, and they are not required for
configuration-driven CQMOM onboarding.

### Starting from a full moment tensor

Some applications already own full tensor moments, from external data exchange,
validation against analytical tensors, debugging, or interoperability with codes
that work in tensors. Select the moments the hierarchy needs from the tensor and
build a case from those:

```cpp
const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});

const qmomkit::RawMomentTensor moments{
    {4, 4},
    {
        1.0, 0.44, 0.424, 0.5552,
        0.55, 0.242, 0.2332, 0.30536,
        1.2125, 0.5335, 0.5141, 0.67318,
        1.121875, 0.493625, 0.475675, 0.622865,
    }};

const qmomkit::CqmomLeanCase case_{
    qmomkit::extract_cqmom_moments(moments, qmomkit::create_cqmom_moment_layout(hierarchy))};
const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);
```

`extract_cqmom_moments(...)` throws if the tensor is the wrong rank or too small
to contain a moment the hierarchy needs. If the tensor comes from outside the
application and you would rather check a status than catch an exception, load or
build a case first and replace its moments with
`config::with_initial_moments(case, tensor)`, which performs the same selection
and returns a configuration error instead.

### Explicit lean representation workflow

Some applications can produce exactly the raw moments required by a CQMOM
hierarchy and want direct control over that compact required-moment
representation. In those cases, use `CqmomMomentLayout` to describe the
required moment multi-indices, store values in `CqmomMomentData`, and construct
`CqmomLeanCase` explicitly.

The explicit lean representation workflow is an expert workflow for callers
that intentionally manage the required CQMOM moment representation. Lean is not
a different physical model; it is the compact representation used for the
moments required by the configured hierarchy and conditioning order. Normal
configuration workflows build and hide this representation automatically.

To evolve a case in time, create its runtime state with
`CqmomPreparedCase::initial_state()` and pass that `CqmomState` to
`advance_cqmom_phase_space(...)`, which performs a complete step. If you would
rather drive the stages yourself, the same state goes through
`advance_cqmom_advection(...)` and `advance_cqmom_moment_source(...)`.

The state carries your configuration forward between calls — reconstruction
settings and, in adaptive mode, the node counts each branch has settled on — so
you do not re-supply them each step. `CqmomState::quadrature()` gives the
currently accepted quadrature at any point. See
[Phase-space evolution](phase_space_evolution.md) for the step itself.
