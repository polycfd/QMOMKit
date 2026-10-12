# Getting started

This page walks through three short programs, each building on the one before:

1. reconstruct a quadrature from four moments;
2. evolve that population in time;
3. do the same for a two-dimensional population defined in a configuration
   file.

To build one, install QMOMKit as described in [Installation](1_installation.md)
and put the program next to a `CMakeLists.txt` such as:

```cmake
cmake_minimum_required(VERSION 3.16)
project(getting_started LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(step1 step1.cpp)
target_link_libraries(step1 PRIVATE qmomkit::qmomkit)
```

Change the file name for the other programs, and configure with
`-DCMAKE_PREFIX_PATH` pointing to the QMOMKit installation.

## 1. Reconstruct a quadrature

Consider a population of two groups: 0.4 at $\xi = -1$ and 0.6 at $\xi = 2$.
Its first four raw moments are

```math
m_0 = 1, \qquad m_1 = 0.8, \qquad m_2 = 2.8, \qquad m_3 = 4.4 .
```

Given only these moments, QMOM recovers the two groups:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    // Raw moments m0..m3 of weights 0.4 and 0.6 at abscissas -1 and 2.
    const qmomkit::MomentSet moments({1.0, 0.8, 2.8, 4.4});

    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, 2);
    if (!result.ok())
    {
        std::cerr << result.message << '\n';
        return 1;
    }
    for (std::size_t i = 0; i < result.weights.size(); ++i)
    {
        std::cout << "weight " << result.weights[i] << " at " << result.abscissas[i] << '\n';
    }
}
```

```text
weight 0.4 at -1
weight 0.6 at 2
```

`compute_qmom(...)` takes the moments and the number of nodes. $N$ nodes need
the $2N$ moments $m_0, \ldots, m_{2N-1}$, so four moments allow two nodes. The
result holds the weights and the abscissas, the node locations in phase space,
in ascending order. If the reconstruction fails, `ok()` is false and `message`
says why, so always check it before using a result.

## 2. Evolve it in time

Now let the population change. Its nodes drift towards zero at a rate
proportional to their abscissa, and the population as a whole decays:

```math
\frac{\mathrm{d}\xi}{\mathrm{d}t} = -0.5\,\xi ,
\qquad
\frac{\mathrm{d}m_k}{\mathrm{d}t} = -0.1\,m_k .
```

The drift moves the nodes, so it is a *velocity*. The decay changes the
moments, so it is a *source*. You supply both as callbacks:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    qmomkit::QmomCase qmom_case;
    qmom_case.moments = qmomkit::MomentSet({1.0, 0.8, 2.8, 4.4});
    qmom_case.node_count = 2;
    qmomkit::QmomState state = qmom_case.initial_state();

    // Nodes drift towards zero: dxi/dt = -0.5 xi.
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi) { return std::vector<double>{-0.5 * xi[0]}; });

    // The population decays: dm_k/dt = -0.1 m_k.
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size())};
        for (std::size_t k = 0; k < context.size(); ++k)
        {
            rate.values[k] = -0.1 * context.value(k);
        }
        return rate;
    };

    double time = 0.0;
    for (int step = 0; step < 10; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, time, 0.1, velocity, source);
        if (!result.ok())
        {
            std::cerr << result.message << '\n';
            return 1;
        }
        state = result.state;
        time = result.final_time;
    }
    const std::vector<double>& m = state.moment_values();
    std::cout << "t = " << time << ": m0 = " << m[0] << ", mean = " << m[1] / m[0] << '\n';
}
```

```text
t = 1: m0 = 0.904837, mean = 0.485224
```

The program has three parts:

- a *case*, the moments and the node count, from which `initial_state()`
  creates the *state* that is evolved: moments together with their quadrature;
- the two callbacks. The velocity receives the coordinates of one node and
  returns how fast it moves. The source returns the rate of change of each
  moment; `context.value(k)` is the current $m_k$;
- the time loop. Each call to `advance_qmom_phase_space(...)` advances the
  state by one step and returns the new state and the time it has reached.

The results match the exact solution: the decay leaves $m_0 = e^{-0.1}$ after
one time unit, and the drift scales every abscissa, and so the mean, by
$e^{-0.5}$, giving $0.8\,e^{-0.5}$.

The loop and the step size are yours.

## 3. Two dimensions from a configuration file

A population can depend on several coordinates, here size and velocity. With
more than one coordinate, QMOMKit uses CQMOM. Such cases are usually described
in a `.qmom` configuration file:

```ini
[phase_space]
names = size, velocity
node_counts = 2, 2

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

`[phase_space]` names the two coordinates and asks for two nodes in each.
`[initial]` describes the initial population: a total number density of 2.5,
with sizes following a gamma distribution and velocities a normal distribution,
independently of each other. QMOMKit generates the moments CQMOM needs from
these distributions.

Save it as `configuration.qmom` in the directory you run the program from:

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
    qmomkit::CqmomState state = configured.value().initial_state();

    // Sizes shrink and velocities relax: d(size)/dt = -0.1 size, d(velocity)/dt = -velocity.
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi) { return std::vector<double>{-0.1 * xi[0], -xi[1]}; });
    // No source in this example.
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    double time = 0.0;
    for (int step = 0; step < 10; ++step)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, time, 0.1, velocity, source);
        if (!result.ok())
        {
            std::cerr << result.message << '\n';
            return 1;
        }
        state = result.state;
        time = result.final_time;
    }

    const qmomkit::CqmomQuadrature& quadrature = state.quadrature();
    for (std::size_t i = 0; i < quadrature.total_nodes(); ++i)
    {
        std::cout << "weight " << quadrature.weights()[i] << " at size " << quadrature.abscissas()[0][i] << ", velocity "
                  << quadrature.abscissas()[1][i] << '\n';
    }
}
```

```text
weight 0.9375 at size 0.361935, velocity -0.643788
weight 0.9375 at size 0.361935, velocity 0.275909
weight 0.3125 at size 1.0858, velocity -0.643788
weight 0.3125 at size 1.0858, velocity 0.275909
```

Compared with the one-dimensional program:

- `load_cqmom_prepared(...)` reads the file and reports a problem with it
  through `configured.error()`, with the line number where there is one;
- the phase-space velocity returns two values, one per coordinate, in the order
  the file declares them;
- this population has no source, which is expressed as a source that returns
  zero rates;
- the result is a set of four nodes, two sizes times two velocities, each with
  one weight and one value per coordinate.

Initially the nodes were at sizes 0.4 and 1.2 and velocities $-1.75$ and
$0.75$. They moved as the velocity prescribes, sizes by $e^{-0.1}$ and
velocities by $e^{-1}$, while the weights stayed the same.

## Next

[Concepts](3_concepts.md) explains the ideas behind these programs.
[Moments and quadrature](4_moments_and_quadrature.md) and
[Evolving in time](5_evolving_in_time.md) cover each step in more depth.
