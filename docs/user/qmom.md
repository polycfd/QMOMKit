# One-Dimensional QMOM

QMOM, the Quadrature Method of Moments, represents a one-dimensional
distribution by a small set of weighted quadrature nodes. Instead of carrying a
full distribution function, a caller provides raw moments such as `m0`, `m1`,
and `m2`. QMOMKit reconstructs weights and abscissas whose raw moments match
the supplied sequence.

This is useful when an application naturally stores or transports moments, but
later needs representative nodes for evaluation, evolution, output, or
coupling to another numerical step.

In a modelling workflow, the split of responsibilities is simple. The
application decides what the coordinate represents, how moments are obtained,
and how reconstructed nodes are used. QMOMKit checks the numerical request,
performs the inversion, keeps weights aligned with sorted abscissas, and
reports whether the requested quadrature was accepted.

QMOMKit uses raw, unnormalized moments:

```text
m_k = sum_i w_i x_i^k
```

`m0` is the total quadrature weight or number density. It is not normalized to
one unless you supply moments with `m0 = 1`.

## QMOM Workflow

The usual one-dimensional workflow begins with the quantity your application
already knows how to produce: raw moments of the number density function. QMOM
then replaces that moment sequence with a small quadrature rule that preserves
the supplied moments up to the requested order.

```text
raw moments
    ↓
moment set
    ↓
quadrature inversion
    ↓
nodes and weights
```

The inputs and outputs are:

- `moments`: raw moment values in ascending order, `m0, m1, m2, ...`.
- `node_count`: requested number of quadrature nodes.
- `weights`: quadrature weights, whose sum reproduces `m0`.
- `abscissas`: node locations associated with the weights.

After reconstruction, a QMOM workflow may inspect the quadrature, evaluate
derived quantities, write output, or optionally evolve the distribution. QMOM
transport advances an accepted `QmomState` from `QmomCase::initial_state()`
with `advance_transport_qmom(...)`. `QmomState` exposes read-only `moments()`
and `quadrature()` accessors; its stored moments are canonicalized from its
accepted quadrature. Users may copy and move accepted states, but direct
construction from independent moments and quadrature is intentionally
unavailable.

For `N` nodes, provide at least `2 * N` moments:

```text
m0, m1, ..., m(2N-1)
```

Returned abscissas are sorted in ascending order. Weights are stored in the
same order, so `weights[i]` belongs to `abscissas[i]`.

Choose `node_count` according to the resolution your closure or downstream
model needs. More nodes require more moments and can represent richer
one-dimensional structure; fewer nodes are cheaper and may be more robust when
only low-order moments are trustworthy.

## Using QMOMKit

Include the umbrella header when following the user guide examples:

```cpp
#include <qmomkit/qmomkit.hpp>
```

Then create a `MomentSet`, request a node count, and check the result status:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    // Raw moments m0..m3 for weights {0.4, 0.6}
    // at abscissas {-1.0, 2.0}.
    const qmomkit::MomentSet moments({1.0, 0.8, 2.8, 4.4});
    const std::size_t node_count = 2;

    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, node_count);
    if (!result.ok())
    {
        std::cerr << "QMOM failed: " << result.message << '\n';
        return 1;
    }

    for (std::size_t node = 0; node < result.weights.size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << result.weights[node]
                  << ", abscissa=" << result.abscissas[node] << '\n';
    }
}
```

Use explicit moments when another code already owns the moment calculation, or
when measured or tabulated moment data is the starting point. When the initial
condition is easier to describe as an analytical distribution, QMOMKit can
generate the raw moments first; that workflow is covered separately in the
distribution guide.

## Moment reconstruction from quadrature

Native QMOM and CQMOM moments remain integer-order polynomial moments. QMOMKit
reconstructs integer-order moments from quadrature when it needs raw moments
for reconstruction handoff, transport, output, or validation workflows. These
moments use orders such as `m0`, `m1`, and `m2`, or integer-valued
multidimensional raw-moment indices.

After checking that a quadrature result is valid, callers may also evaluate
generalized moments with arbitrary finite real exponents as derived
quantities:

```cpp
if (!result.ok())
{
    return 1;
}

const qmomkit::Quadrature1D quadrature{
    result.abscissas,
    result.weights,
    qmomkit::QuadratureStatus::success,
    ""};

const double m_one_third =
    qmomkit::evaluate_generalized_moment(quadrature, 1.0 / 3.0);
```

Generalized moment evaluation applies the quadrature formula directly. The same
idea is available for CQMOM quadrature with a phase-space-order exponent vector.
These values are post-processing quantities, not QMOM or CQMOM inversion inputs,
transported state, native `MomentSet` values, `CqmomMomentData`, or
`RawMomentTensor` entries. Non-integer exponents require non-negative abscissas.

## Options And Diagnostics

The simplest call uses default numerical options. This is the right starting
point when the requested number of nodes is fixed by the model:

```cpp
const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, node_count);
```

For applications that can accept a lower-order rule near difficult moment
states, adaptive node reduction can be enabled with `QmomOptions`:

```cpp
qmomkit::QmomOptions options;
options.adaptive = true;

const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, node_count, options);
```

`result.ok()` identifies whether the quadrature is a successful
reconstruction. A failed result reports a status and message so the caller can
decide whether to reject a state, reduce a timestep, change a requested node
count, or inspect the supplied moments:

- `success`: the quadrature was reconstructed.
- `invalid_input`: the input was structurally invalid, such as too few
  moments or a zero node count.
- `not_realizable`: the supplied moment sequence could not represent a valid
  quadrature under the configured reconstruction.
- `no_adaptive_quadrature`: adaptive reconstruction could not find an accepted
  multi-node rule.

Moment realizability is a mathematical property of the supplied sequence.
QMOMKit does not silently change moments to make them realizable. If your
workflow can approach the boundary of moment space, use explicit status
handling and consider the optional realizability checks described in the
realizability user guide.
