# QMOM

QMOM, the Quadrature Method of Moments [[McGraw,
1997](https://www.doi.org/10.1080/02786829708965471)], represents a
one-dimensional distribution by a small set of weighted quadrature nodes.
Instead of carrying a full distribution function, a caller provides raw moments
such as $m_0$, $m_1$, and $m_2$. QMOMKit reconstructs weights and abscissas
whose raw moments match the supplied sequence using the Wheeler algorithm
[[Wheeler, 1974](https://www.doi.org/10.1216/RMJ-1974-4-2-287)].

QMOMKit uses raw, unnormalized moments:
$$m_k = \sum_i w_i x_i^k$$

$m_0$ is the total quadrature weight or number density.

## QMOM workflow

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

- `moments`: raw moment values in ascending order, $m_0, m_1, m_2, ...$.
- `node_count`: requested number of quadrature nodes.
- `weights`: quadrature weights, whose sum reproduces $m_0$.
- `abscissas`: node locations associated with the weights.

After reconstruction, a QMOM workflow may inspect the quadrature, evaluate
derived quantities, write output, or optionally evolve the distribution. QMOM
phase-space evolution advances an accepted `QmomState` from `QmomCase::initial_state()`
with `advance_qmom_phase_space(...)`. `QmomState` exposes read-only `moment_values()`,
`quadrature()`, `maximum_node_count()`, and `qmom_options()` accessors. Its
moments and quadrature always agree with each other, and it keeps the
reconstruction settings of the case it came from. You can copy and move an
accepted state, but you cannot build one from a moment vector and a quadrature
of your own choosing: that would allow the two to disagree. When your
application changes the moments between steps, for example by transport in
physical space, `QmomCase::state_candidate(...)` creates the state of the new
moments; see
[Moments changed by the host](phase_space_evolution.md#moments-changed-by-the-host).

For $N$ nodes, provide at least $2N$ moments: $m_0, m_1, ..., m_{2N-1}$


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
for reconstruction handoff, phase-space evolution, output, or validation
workflows. These moments use orders such as $m_0$, $m_1$, and $m_2$, or
integer-valued multidimensional raw-moment indices.

Each contribution $w_i\,x_i^k$ is evaluated with the power $x_i^k$ formed
first, so that power must itself lie within the range of a double, about
$10^{\pm 308}$. A power that overflows makes the reconstruction fail with
`std::overflow_error`, and a power that underflows contributes zero, even
where $w_i\,x_i^k$ would be representable. This arises only for coordinates
far from unity at high orders, for example droplet volumes of about
$10^{-18}\,\mathrm{m}^3$ from order 18 on. Moment sets like that span hundreds
of decades and are poorly suited to inversion in any case, so express such a
coordinate in units that keep it near unity, for example cubic micrometres
instead of cubic metres.

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
idea is available for CQMOM quadrature with a phase-space-order exponent
vector. These values are derived post-processing quantities. Non-integer
exponents require non-negative abscissas.

## Options and diagnostics

The simplest call uses the default numerical options, which reduce the node
count adaptively: near degenerate moment states, or where a node's weight falls
below the minimum weight ratio `rmin` relative to the largest, the result has
fewer nodes than requested.

```cpp
const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, node_count);
```

When the model needs exactly the requested number of nodes, request a fixed
count; the reconstruction then returns every node or fails. It fails with
`not_realizable` when the moments, within their rounding, support fewer nodes
than requested:

```cpp
qmomkit::QmomOptions options;
options.adaptive = false;

const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, node_count, options);
```

Adaptive reduction compares the recurrence with a bound on its own rounding,
so the number of nodes it keeps depends neither on the unit of the coordinate
nor on the number density.

### Optional abscissa bounds

`AbscissaBounds` optionally requires returned QMOM nodes to lie in a closed
interval. QMOMKit still performs ordinary Wheeler inversion; it accepts the
candidate only when every returned node is inside the interval within the
numerical reconstruction tolerance. Nodes are never clipped or projected.

```cpp
const qmomkit::QmomResult result =
    qmomkit::compute_qmom(moments, node_count, {}, qmomkit::AbscissaBounds{0.0, 1.0});
```

In fixed-count mode, an out-of-interval candidate fails. In adaptive mode,
QMOMKit tries lower explicit node counts and accepts the highest admissible
count. These bounds constrain reconstructed abscissas only; they are separate
from the `MomentSupport` choices used by realizability checks.

Store bounds in `QmomCase::abscissa_bounds` when creating an evolving state.
The accepted `QmomState` retains them across advection and source evolution.
Advection itself does not clip or reject nodes for crossing the interval;
bounds are enforced at subsequent reconstruction boundaries.

`result.ok()` identifies whether the quadrature is a successful
reconstruction. In adaptive mode, Wheeler may retry candidates at lower node
counts. A valid reduced reconstruction, including a one-node fallback, is a
successful result. A failed result reports a status and message so the caller
can decide whether to reject a state, reduce a timestep, change a requested
node count, or inspect the supplied moments:

- `success`: the quadrature was reconstructed, possibly at a lower node count
  in adaptive mode.
- `invalid_input`: the input was structurally invalid, such as too few
  moments or a zero node count.
- `not_realizable`: the supplied moment sequence could not represent a valid
  quadrature under the configured reconstruction.

Moment realizability is a mathematical property of the supplied sequence.
QMOMKit does not silently change moments to make them realizable. If your
workflow can approach the boundary of moment space, use explicit status
handling and consider the optional realizability checks described in the
realizability user guide.

## Next steps

For a population described by more than one coordinate, CQMOM extends this
reconstruction to multiple dimensions by conditioning them in turn — see
[Multidimensional CQMOM](cqmom.md). To evolve a one-dimensional population in
time, see [Phase-space evolution](phase_space_evolution.md), and when a
reconstruction is rejected, [Troubleshooting](troubleshooting.md).
