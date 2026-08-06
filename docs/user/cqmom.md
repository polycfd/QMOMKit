# Multidimensional CQMOM

CQMOM, the Conditional Quadrature Method of Moments, extends moment-based
quadrature reconstruction to multidimensional phase spaces. A multidimensional
distribution may depend on coordinates such as size, velocity, or any other
application-defined variables. QMOMKit treats those coordinates generically as
phase-space dimensions.

Use CQMOM when the modelled population cannot be represented adequately by one
coordinate. The result is still a compact quadrature representation, but each
terminal node has one global weight and one abscissa in every configured
phase-space dimension. This lets another part of the simulation evaluate
multidimensional closures or evolve a coupled set of coordinates while
continuing to store raw moments between reconstructions.

In one dimension, QMOM reconstructs nodes and weights from a raw moment
sequence. In multiple dimensions, CQMOM reconstructs a multidimensional
quadrature from a CQMOM case. The case records the phase-space dimensions,
requested node counts, conditioning order, initial moments, and reconstruction
options so these choices remain consistent across initialization, phase-space evolution,
output, and validation.

QMOMKit uses raw, unnormalized tensor moments:

```text
M[k_0, ..., k_D-1] =
    sum_j W_j product_d xi_(d,j)^k_d
```

The zeroth tensor moment `M[0, ..., 0]` is the total quadrature weight or
number density.

## CQMOM Workflow

The typical CQMOM workflow starts from a modelling description rather than from
internal storage choices. That description can come from a `.qmom` configuration
file or from direct C++ setup. In either case, users name or order the
coordinates, set the requested node counts, choose the conditioning order, and
describe the initial moments. File-based setup builds a `CqmomCase`; direct C++
setup can provide advanced explicit-representation case types such as
`CqmomDenseCase` or `CqmomLeanCase` when the application intentionally manages
dense or lean moment storage.

```text
configuration file:
    .qmom file
    -> load_cqmom(...)
    -> CqmomCase

direct API:
    C++ setup
    -> CqmomDenseCase or CqmomLeanCase

reconstruction:
    CqmomCase, CqmomDenseCase, or CqmomLeanCase
    -> compute_cqmom(...)
    -> CqmomResult / CqmomQuadrature
    -> inspect / evaluate / output

normal optional phase-space evolution:
    CqmomCase::initial_state()
    -> CqmomState
    -> advance_cqmom_phase_space(...)
```

At the user level, this means:

1. Define the CQMOM setup in a `.qmom` file or directly in C++.
2. Load it with `load_cqmom(...)` or construct an explicit public CQMOM case in
   C++.
3. Keep the case object as the configured problem definition.
4. Call `compute_cqmom(...)` to compute the initial quadrature.
5. Inspect the quadrature, evaluate derived quantities, or write output.
6. If the case evolves in time, call `CqmomCase::initial_state()` to create the
   runtime state for phase-space evolution and pass that `CqmomState` to
   `advance_cqmom_phase_space(...)`.

The result is a set of weighted multidimensional nodes. Each terminal node has
one global weight and one abscissa per phase-space dimension.

`CqmomCase` is the representation-neutral immutable setup produced by the normal
configuration loader. Normal CQMOM phase-space evolution creates `CqmomState` from
`CqmomCase` and advances it with `advance_cqmom_phase_space(...)`, independent of
whether the case stores moments densely or leanly. Advanced direct C++ workflows
can instead pass `CqmomDenseCase` or `CqmomLeanCase` to `compute_cqmom(...)` for
explicit static inversion. In all cases, `compute_cqmom(...)` returns
quadrature; it does not return a phase-space evolution state.

For the usual CQMOM path, users specify the phase space, initial condition, and
reconstruction settings. QMOMKit determines which moments are needed for the
configured hierarchy and avoids exposing the internal moment storage.

### Replacing Initial Moments Programmatically

When an application obtains a complete `RawMomentTensor` after loading its
ordinary CQMOM setup, it can make a separate phase-space evolution-capable case without
rewriting configuration text:

```cpp
const auto replaced = qmomkit::config::with_initial_moments(base_case, std::move(raw_moments));
if (!replaced)
{
    std::cerr << replaced.error() << '\n';
    return 1;
}
const qmomkit::CqmomCase cqmom_case = replaced.value();
```

The new case retains the base hierarchy, conditioning order, support, storage
selection, and every reconstruction policy. A dense base retains the complete
tensor; a lean base extracts its fixed required layout. The transformation checks
rank, extents, value count, finiteness, and required coverage only. It does not
test realizability: `compute_cqmom(...)` reports numerical reconstruction
status, and `initial_state()` keeps its existing exception-on-failure behavior.
The producer of the raw moments can be a sampler, measurement pipeline, or any
other application code.

## Using CQMOM In QMOMKit

Include the umbrella header when following the user guide examples:

```cpp
#include <qmomkit/qmomkit.hpp>
```

The example below loads a two-dimensional CQMOM case from configuration and
computes the initial quadrature. It stops at reconstruction; phase-space evolution is shown
afterward as an optional step.

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    const auto configured = qmomkit::config::load_cqmom("configuration.qmom");
    if (!configured)
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    const qmomkit::CqmomCase cqmom_case = configured.value();

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

If this reconstructed distribution evolves in time, create the phase-space evolution state
from the same case and pass it to `advance_cqmom_phase_space(...)`:

```cpp
qmomkit::CqmomState state = cqmom_case.initial_state();
```

For a simple independent-product setup, the configuration file can look like
this. The example describes a two-coordinate population with one global number
density and one analytical marginal distribution per phase-space dimension:

```text
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

### Choosing The Vandermonde Solver

Recursive CQMOM reconstructs conditional quantities by solving branch-local
conditional Vandermonde systems. QMOMKit provides several solver policies for
those systems through the optional `[vandermonde]` section:

```ini
[vandermonde]
solver = automatic
```

This line may be omitted when using the default. Supported values are:

| Value | Practical meaning |
|---|---|
| `automatic` | Recommended general default; currently uses Björck-Pereyra with eligible dense fallback. |
| `dense` | Explicit dense Gaussian-elimination solve. |
| `pivoted_qr` | Explicit pivoted-QR solve with numerical-rank detection and no fallback. |
| `svd` | Explicit SVD solve with rank and convergence checks and no fallback. |
| `bjorck_pereyra` | Explicit structured Vandermonde solve with no fallback. |
| `bjorck_pereyra_with_dense_fallback` | Structured solve followed by eligible dense fallback. |

QR and SVD rank-deficient systems fail, and SVD nonconvergence also fails.
QMOMKit does not silently return a pseudoinverse or regularized solution in
these paths. See [configuration.md](configuration.md) for exact syntax and
policy behavior. Developer-level numerical details are in
[Internal Linear Algebra](../developer/linear_algebra.md).

For this normal configuration path, the CQMOM hierarchy determines the
required raw moments automatically. Advanced full-tensor workflows may also
omit `maximum_orders`; QMOMKit then derives the minimum dense extents from the
configured quadrature, while explicit larger extents remain available.

For fixed rectangular quadratures, you can also access a node by its
phase-space-order multi-index:

```cpp
const qmomkit::MultiIndex node_index{0, 1};
const double weight = result.quadrature.weight(node_index);
const double x0 = result.quadrature.abscissa(0, node_index);
const double x1 = result.quadrature.abscissa(1, node_index);
```

## Conditioning And Ordering

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

## Fixed And Adaptive Quadrature

The default CQMOM setup uses a fixed branch quadrature. This is the easiest
mode to reason about: each branch must produce the requested number of nodes
for its phase-space dimension, and fixed results support rectangular node
access by multi-index.

Adaptive branch quadrature is available when a workflow can tolerate fewer
active nodes on branches where the supplied moments support only a lower-order
rule. This can be useful near degenerate or boundary-like states where a fixed
node count may be too ambitious. Enable it explicitly in configuration with
adaptive branch counts:

```text
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
configured `CqmomState` phase-space evolution workflow. `CqmomState` keeps the required
moment representation for the configured hierarchy, while each reconstructed
quadrature may contain adaptive active node counts.

As with QMOM, always check `result.ok()` before using the quadrature. If CQMOM
fails, inspect `result.status`, `result.message`, and diagnostics to determine
which phase-space dimension or conditioning step was not accepted.

## Moment Reconstruction From Quadrature

Native CQMOM moments remain integer-order polynomial raw moments with
phase-space-order multi-indices. QMOMKit reconstructs those moments from
terminal quadrature when output, validation, or phase-space evolution handoff needs raw
moment values.

After checking that a CQMOM quadrature result is valid, callers may also
evaluate generalized moments with finite real exponent vectors as derived
post-processing quantities:

```cpp
const std::vector<double> exponents{1.0 / 3.0, 2.0};
const double derived =
    qmomkit::evaluate_generalized_moment(result.quadrature, exponents);
```

Generalized moment evaluation applies the terminal quadrature formula directly.
These values are not CQMOM reconstruction inputs, phase_space_evolved state, or native
`CqmomMomentData`/`RawMomentTensor` entries. Non-integer exponents require
non-negative abscissas in the corresponding dimensions.

## Advanced Representation Workflows

The workflows in this section are advanced explicit dense and lean
representations for static inversion, specialized initialization, validation,
and expert workflows. They are for applications that deliberately manage moment
storage themselves, and they are not required for configuration-driven CQMOM
onboarding.

### Explicit Dense Representation Workflow

Some applications already own full tensor moments. This is useful for external
data exchange, validation against analytical tensors, debugging, and
interoperability with dense tensor codes. In those cases, construct a
`RawMomentTensor` and use `CqmomDenseCase` explicitly:

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

const qmomkit::CqmomDenseCase dense_case{hierarchy, moments};
const qmomkit::CqmomResult result = qmomkit::compute_cqmom(dense_case);
```

This path is not the recommended beginner workflow. It is the explicit dense
representation workflow for callers that deliberately need full rectangular
tensor input. Dense remains a first-class supported expert workflow when the
application intentionally manages dense moment storage.

### Explicit Lean Representation Workflow

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

For phase-space evolution, use `CqmomCase::initial_state()` and pass the resulting
`CqmomState` to `advance_cqmom_phase_space(...)` when QMOMKit should own the
configured macro-step sequence. Applications that own stage orchestration can
pass the same `CqmomState` through `advance_cqmom_advection(...)` and
`advance_cqmom_moment_source(...)`. `CqmomState` preserves the configured
reconstruction policy, selected moment storage, and adaptive branch-count
topology between successful phase-space evolution calls. Dense and lean case types remain
available for callers that intentionally manage static inversion inputs, but
public CQMOM phase-space evolution is storage-independent.
