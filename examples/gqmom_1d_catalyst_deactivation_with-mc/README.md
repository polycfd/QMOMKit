# 1D GQMOM catalyst deactivation with a Monte Carlo reference

This example evolves a population of catalyst particles whose surface cokes up
over time, and which are withdrawn for regeneration once they are heavily
coked. It compares GQMOM, with its `beta` continuation, and classical QMOM.
`main.cpp` runs both from the raw moments of a finite sample of particles, and
`main_mc.cpp`, the Monte Carlo reference, evolves the sampled particles
themselves with the Monte Carlo companion. Particles are withdrawn only above a
critical coked fraction, and the many nodes of GQMOM follow the particles
above it far better than the few nodes of classical QMOM.

## Model

The phase-space coordinate is the coked fraction $\theta \in [0, 1]$ of a
particle's surface and $n(\theta)$ is the number density; time $t$ and all
parameters are dimensionless. The raw moments are
$m_k = \int \theta^k n(\theta)\,\mathrm{d}\theta$. Coke builds up on every
particle as

```math
\frac{\mathrm{d}\theta}{\mathrm{d}t} = \kappa\,(1 - \theta),
```

with $\kappa = 0.1$, which keeps every particle in $[0, 1]$. A particle is
withdrawn at the frequency

```math
g(\theta) = \lambda\,\max(0,\ \theta - \theta_c),
```

with the critical coked fraction $\theta_c = 0.8$ and $\lambda = 20$, so
particles below $\theta_c$ stay. This is a model problem, a simple picture of a
catalyst whose particles are sent to regeneration once their activity is too
low. QMOMKit applies coking by moving the nodes, and removal changes the
moments as

```math
\frac{\mathrm{d}m_k}{\mathrm{d}t} = -\int g(\theta)\,\theta^k n(\theta)\,\mathrm{d}\theta.
```

$g$ is zero below $\theta_c$, so the right-hand side is not a combination of the
moments, and the moment source sums it over the nodes of the quadrature QMOMKit
reconstructs at each Runge–Kutta stage. Each step of $\Delta t$ is split as
coking over $\Delta t/2$, removal over $\Delta t$ and coking over $\Delta t/2$.
In the Monte Carlo run, each particle is removed within the step with
probability $1 - \exp(-g(\theta)\,\Delta t)$. The parameters, the velocity, the
moment source and the stochastic removal are in [physics.hpp](./physics.hpp)
and [physics.cpp](./physics.cpp). Both drivers run from $t = 0$ to $10$ with
$\Delta t = 0.1$.

[configuration.qmom](./configuration.qmom) defines the initial number density
as a beta distribution with number density $1$, $\alpha = 2$ and $\beta = 3$,
$n(\theta) \propto \theta(1 - \theta)^2$, with mean coked fraction $0.4$. It
sets the node count $N = 2$, so that the GQMOM state holds $m_0, \ldots, m_4$,
and the `generalized` family with $N_G = 20$ nodes and the `beta` continuation,
which is for a number density on $[0, 1]$. `main.cpp` evolves this case and
the same case with the `classical` family, which has $N$ nodes and evolves
$m_0, \ldots, m_3$. Both executables draw the same $50000$ samples
$\theta_p$ with a fixed seed ([sampling.hpp](./sampling.hpp)) and give each the
weight $1/50000$, and `main.cpp` replaces the moments of the configured case
with those of the sample.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./gqmom_1d_catalyst_deactivation_with_mc_example
./gqmom_1d_catalyst_deactivation_with_mc_mc_example
```

## Output

The QMOMKit run prints the initial and final raw moments of GQMOM,
$m_0, \ldots, m_4$, and of classical QMOM, $m_0, \ldots, m_3$, labeled `M0` to
`M4`. The Monte Carlo run prints those of the particles, $m_0, \ldots, m_4$,
and the initial and final numbers of particles. Every $1$ in time, the QMOMKit
run writes the moments and quadratures of both states to `build/output/`, in
`gqmom_moment_history.txt`, `qmom_moment_history.txt`,
`gqmom_quadrature_history.txt` and `qmom_quadrature_history.txt`, and the Monte
Carlo run writes `build/output/moment_history_mc.txt`. The `M[4]` of the
classical QMOM history is the moment of its two nodes, which classical QMOM does
not evolve.

## What to expect

All three runs start from the same sample, so their initial moments agree. At
$t = 10$:

```text
         GQMOM      QMOM       Monte Carlo
  M0     0.760858   0.602164   0.76288
  M1     0.569816   0.435785   0.571022
  M2     0.428846   0.316073   0.429441
  M3     0.32432    0.22984    0.324469
  M4     0.246441   —          0.246272
```

GQMOM is within 0.4% of the Monte Carlo moments. Classical QMOM ends with 21%
fewer particles and moments 24 to 29% lower than the Monte Carlo run. Its two
nodes start at $\theta = 0.23$ and $0.63$, and its upper node reaches
$\theta_c$ only between $t = 6$ and $7$. Until then, classical QMOM removes no
particles, while the particles of the tail above $\theta_c$ are already being
withdrawn. Afterwards, its upper node, which carries 43% of the particles, is
removed almost entirely by $t = 10$. The 20 nodes of GQMOM start between
$0.007$ and $0.986$ and follow the tail. The Monte Carlo run is a stochastic
reference: its random removal changes its moments by about 0.1 to 0.2% from
run to run.
