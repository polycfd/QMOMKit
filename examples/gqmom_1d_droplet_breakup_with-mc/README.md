# 1D GQMOM droplet breakup with a Monte Carlo reference

This example evolves a population of droplets that break up when they are
larger than a critical volume. It compares GQMOM with its `gamma`
continuation, with classical QMOM. `main.cpp` runs both from the raw
moments of a finite sample of droplets, and `main_mc.cpp`, the Monte Carlo
reference, breaks the sampled droplets themselves with the Monte Carlo
companion. The breakup frequency is zero below the critical volume, and the
many nodes of GQMOM follow the droplets above it far better than the few nodes
of classical QMOM.

## Model

The phase-space coordinate is the droplet volume $V \ge 0$ and $n(V)$ is the
number density; $V$, time $t$ and all parameters are dimensionless. The raw
moments are $m_k = \int V^k n(V)\,\mathrm{d}V$. A droplet of volume $V$ breaks
into two droplets of volume $V/2$ at the frequency

```math
g(V) = \lambda\,\max(0,\ V - V_c),
```

with the critical volume $V_c = 1$ and $\lambda = 1$, so droplets smaller than
$V_c$ are stable. This is a model problem, a simple picture of droplets in a
turbulent flow, where only droplets above a critical size break up. Breakup
conserves the total volume $m_1$, and the moments change as

```math
\frac{\mathrm{d}m_k}{\mathrm{d}t} = \left(2^{1-k}-1\right)\int g(V)\,V^k n(V)\,\mathrm{d}V.
```

$g$ is zero below $V_c$, so the right-hand side is not a combination of the
moments, and QMOMKit sums it over the nodes of the quadrature it reconstructs
at each Runge–Kutta stage. In the Monte Carlo run, each droplet breaks within a
step $\Delta t$ with probability $1 - \exp(-g(V)\,\Delta t)$. $g$, $V_c$,
$\lambda$, the moment source and the stochastic breakup are in
[physics.hpp](./physics.hpp) and [physics.cpp](./physics.cpp). Both drivers
run from $t = 0$ to $8$ with $\Delta t = 0.01$.

[configuration.qmom](./configuration.qmom) defines the initial number density
as a gamma distribution with number density $1$, shape $4$ and scale $0.75$,
$n(V) \propto V^3 e^{-V/0.75}$, with mean volume $3$. It sets the node count
$N = 2$, so that the GQMOM state holds $m_0, \ldots, m_4$, and the `generalized`
family with $N_G = 20$ nodes and the `gamma` continuation, which is for a
number density on $[0, \infty)$. `main.cpp` evolves this case and the same case
with the `classical` family, which has $N$ nodes and evolves
$m_0, \ldots, m_3$, with source steps only: breakup changes the number of
droplets, not the volume of a droplet between breakups. Both executables draw
the same $50000$ samples $V_p$ with a fixed seed ([sampling.hpp](./sampling.hpp))
and give each the weight $1/50000$, and `main.cpp` replaces the moments of the
configured case with those of the sample.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./gqmom_1d_droplet_breakup_with_mc_example
./gqmom_1d_droplet_breakup_with_mc_mc_example
```

## Output

The QMOMKit run prints the initial and final raw moments of GQMOM,
$m_0, \ldots, m_4$, and of classical QMOM, $m_0, \ldots, m_3$, labeled `M0` to
`M4`. The Monte Carlo run prints those of the droplets, $m_0, \ldots, m_4$, and
the initial and final numbers of droplets. Every $1$ in time, the QMOMKit run
writes the moments and quadratures of both states to `build/output/`, in
`gqmom_moment_history.txt`, `qmom_moment_history.txt`,
`gqmom_quadrature_history.txt` and `qmom_quadrature_history.txt`, and the Monte
Carlo run writes `build/output/moment_history_mc.txt`. The `M[4]` of the
classical QMOM history is the moment of its two nodes, which classical QMOM does
not evolve.

## What to expect

All three runs start from the same sample, so their initial moments agree, and
all three keep $m_1 = 2.99134$. At $t = 8$:

```text
         GQMOM      QMOM      Monte Carlo
  M0     3.82592    2.99162   3.83364
  M2     2.4954     2.99181   2.45957
  M3     2.19463    2.99375   2.13301
  M4     2.02349    —         1.95041
```

GQMOM is 0.2% below the Monte Carlo $m_0$ and within 4% of its $m_2$ to
$m_4$. Classical QMOM ends with 22% fewer droplets and a 40% larger $m_3$. Its
larger node keeps breaking and losing weight, while its other node settles just
below $V_c = 1$, where the breakup frequency vanishes, so classical QMOM ends
with nearly all its droplets near volume $1$ and its $m_0$ stops near $m_1$.
The droplets themselves spread over volumes below $V_c$: each one above $V_c$
keeps breaking into halves until they are stable, which the 20 nodes of GQMOM
resolve. Raising `generalized_node_counts` in `build/configuration.qmom`
changes the GQMOM results little, so the remaining difference is mostly that of
the closure itself. The Monte Carlo run is a stochastic reference with sampling
noise, and it lets each droplet break at most once per step, an error of first
order in $\Delta t$, so it is not an exact solution: halving $\Delta t$ changes
its $m_0$ to $m_4$ by at most 0.6%.
