# QMOMKit Quick Start

QMOMKit turns raw moments into quadrature representations for use in a larger
application or validation workflow. Use the resulting weighted nodes to
evaluate closures, evolve a population through phase space, write output, or
compare with sampled data. Your application supplies the coordinate meanings
and physical rate laws.

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

QMOMKit uses raw, unnormalized moments. In one dimension, $m_0$ is the total
quadrature weight or number density. In multiple dimensions, configured
phase-space dimensions keep a stable order for moment indices, nodes,
callbacks, and output columns.

In practice, a first QMOMKit integration usually asks four questions:

- What moment representation does my model already have?
- Should this case use one-dimensional QMOM or multidimensional CQMOM?
- Are initial moments supplied explicitly, generated from distributions, or
  read from configuration?
- After reconstruction, will the quadrature only be inspected, or will it be
  evolved through phase-space evolution?

## QMOM

Use QMOM [[McGraw, 1997](https://www.doi.org/10.1080/02786829708965471)] when the represented population has one phase-space coordinate and
your model provides a raw moment sequence for that coordinate. QMOMKit
reconstructs a one-dimensional quadrature rule: weights whose sum is $m_0$, and
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

For an $N$-node QMOM inversion, provide at least $2N$ raw moments:
$m_0, m_1, ..., m_{2N-1}$. The returned abscissas are sorted in ascending order, with
the corresponding weights using the same order.

## CQMOM

Use CQMOM [[Yuan & Fox, 2011](https://www.doi.org/10.1016/j.jcp.2011.07.020)] 
when the represented population depends on several phase-space
coordinates, such as a size coordinate and a velocity coordinate. QMOMKit
keeps the phase-space order, CQMOM conditioning order, required moments, and
quadrature layout consistent. This quick start uses a configuration file because
it is concise and reproducible, then loads a CQMOM case and computes its
quadrature. The same setup can also be created directly through the C++ API when
an application already manages setup data in code.

Save this minimal configuration as `configuration.qmom` next to the program
below:

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

`compute_cqmom(...)` returns the quadrature for the configured case. Phase-space evolution is
optional. When the CQMOM distribution evolves in time, create the runtime state
separately:

```cpp
qmomkit::CqmomState state = cqmom_case.initial_state();
```

That state is passed to `advance_cqmom_phase_space(...)` in a phase-space evolution loop.

The configuration file defines the parts of the numerical setup that should be
repeatable: phase-space names, node counts, conditioning order, and either an
initial distribution or explicit moments. It does not define the physics
callbacks used later for phase-space evolution. See [configuration.md](./configuration.md)
 for complete
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
- advance quadrature nodes with a caller-supplied advection velocity;
- advance raw moments with a caller-supplied moment-source term;
- run optional realizability checks;
- compare moment histories against Monte Carlo companion outputs.

`MomentWriter` writes native integer-order raw moments. A state-aware writer
configured with `OutputOptions` can append generalized quadrature-derived
columns. `QuadratureWriter` writes reconstructed nodes and weights. Generalized
moments remain derived quantities evaluated from a valid quadrature result.

See [configuration.md](./configuration.md) for configuration-file examples, 
[phase_space_evolution.md](./phase_space_evolution.md) for time evolution, 
[realizability.md](./realizability.md) for support-aware checks, and
[monte_carlo.md](./monte_carlo.md) for Monte Carlo companion workflows.
