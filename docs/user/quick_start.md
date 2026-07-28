# QMOMKit Quick Start

QMOMKit turns raw moments into quadrature representations that can be used by a
larger application or validation workflow. This is useful when a research code
stores a reduced set of moments, but a modelling step needs representative
weighted nodes: for example to evaluate a closure, move a population through
phase space, write interpretable output, or compare with sampled data.

QMOMKit does not define a physical model, mesh, PDE discretization, solver loop,
or persistent probability-density object. It converts moments to quadrature,
stores the case setup, and, when transport is used, tracks the state advanced by
the transport helper. Your application defines the meaning of the coordinates
and the equations that generate or evolve the moments.

The basic workflow is:

```text
Moment representation
        ↓
Quadrature inversion
        ↓
Quadrature representation
        ↓
Evolution/output
```

QMOMKit uses raw, unnormalized moments. In one dimension, `m0` is the total
quadrature weight or number density. In multiple dimensions, configured
phase-space dimensions keep a stable order for moment indices, nodes,
callbacks, and output columns.

In practice, a first QMOMKit integration usually asks four questions:

- What moment representation does my model already have?
- Should this case use one-dimensional QMOM or multidimensional CQMOM?
- Are initial moments supplied explicitly, generated from distributions, or
  read from configuration?
- After reconstruction, will the quadrature only be inspected, or will it be
  evolved through transport?

## QMOM

Use QMOM when the represented population has one phase-space coordinate and
your model provides a raw moment sequence for that coordinate. QMOMKit
reconstructs a one-dimensional quadrature rule: weights whose sum is `m0`, and
abscissas ordered along the coordinate.

```cpp
#include <qmomkit/qmom.hpp>

#include <iostream>

int main()
{
    // Raw moments m0..m3 for two nodes:
    // weights {0.4, 0.6}, abscissas {-1.0, 2.0}.
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

For an `N`-node QMOM inversion, provide at least `2 * N` raw moments:
`m0, m1, ..., m(2N-1)`. Returned abscissas are sorted in ascending order, and
weights use the same order.

## CQMOM

Use CQMOM when the represented population depends on several phase-space
coordinates, such as a size coordinate and a velocity coordinate. QMOMKit
keeps the phase-space order, CQMOM conditioning order, required moments, and
quadrature layout consistent. This quick start uses a configuration file because
it is concise and reproducible, then loads a CQMOM case and computes its
quadrature. The same setup can also be created directly through the C++ API when
an application already manages setup data in code.

Save this minimal configuration as `configuration.qmom` next to the program
below:

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

`compute_cqmom(...)` returns the quadrature for the configured case. Transport is
optional. When the CQMOM distribution evolves in time, create the runtime state
separately:

```cpp
qmomkit::CqmomState state = cqmom_case.initial_state();
```

That state is passed to `advance_transport_cqmom(...)` in a transport loop.

The configuration file defines the parts of the numerical setup that should be
repeatable: phase-space names, node counts, conditioning order, and either an
initial distribution or explicit moments. It does not define the physics
callbacks used later for transport. See `configuration.md` for complete
configuration syntax and the relationship between `.qmom` files and direct C++
setup.

The quadrature result exposes weights and abscissas by terminal node. For fixed
rectangular quadratures, you can also access a node by its phase-space-order
multi-index:

```cpp
const qmomkit::MultiIndex node_index{0, 1};
const double weight = quadrature.weight(node_index);
const double x0 = quadrature.abscissa(0, node_index);
const double x1 = quadrature.abscissa(1, node_index);
```

## After The First Inversion

Once you have a quadrature representation, a larger workflow can:

- write moments or quadrature histories with QMOMKit output helpers;
- advance quadrature nodes with a caller-supplied continuous velocity;
- advance raw moments with a caller-supplied discontinuous source term;
- run optional realizability checks;
- compare moment histories against Monte Carlo companion outputs.

`MomentWriter` writes native integer-order raw moments. `QuadratureWriter`
writes reconstructed nodes and weights. Generalized moments are derived
quantities evaluated separately from a valid quadrature result.

See `configuration.md` for configuration-file examples, `transport.md` for
time evolution, `realizability.md` for support-aware checks, and
`monte_carlo.md` for Monte Carlo companion workflows.
