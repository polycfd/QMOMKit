# 1D QMOM breakage and advection with a Monte Carlo reference

This example runs the breakage and advection model of
[qmom_1d_breakage_advection](../qmom_1d_breakage_advection/README.md), with a
constant fragmentation rate, twice from the same finite sample of particles.
`main.cpp` starts QMOMKit from the raw moments of the sample and advances them
with `qmomkit::advance_qmom_phase_space(...)`; `main_mc.cpp`, the Monte Carlo
reference, advances the sampled particles themselves with the Monte Carlo
companion.

## Model

The phase-space coordinate is the particle volume $V$ and $n(V)$ is the number
density; $V$, time $t$ and all parameters are dimensionless. The raw moments
are $m_k = \int V^k n(V)\,\mathrm{d}V$. Volumes relax toward $V = 1$,

```math
\frac{\mathrm{d}V}{\mathrm{d}t} = 1 - V,
```

and particles break at the constant rate $\lambda = 0.5$, one particle of
volume $V$ into two of volume $V/2$, which conserves the total volume $m_1$.
The moments change as

```math
\frac{\mathrm{d}m_k}{\mathrm{d}t} = \lambda\left(2^{1-k}-1\right)m_k.
```

In the Monte Carlo run each particle breaks within a step $\Delta t$ with
probability $1-\exp(-\lambda\,\Delta t)$. Both runs use Strang splitting: half
a step of advection, the breakage over the whole step, and a second half step
of advection. The velocity and the moment source are in
[physics.cpp](./physics.cpp), $\lambda$ in [physics.hpp](./physics.hpp), the
stochastic breakage in [particles_mc.cpp](./particles_mc.cpp), and the time
span $0 \le t \le 1$ with $\Delta t = 0.01$ in both drivers.

[configuration.qmom](./configuration.qmom) defines the initial number density
as two normal modes with number densities $0.5$ each, means $1$ and $3$ and
standard deviations $0.1$ and $0.2$, so the total number density is
$n_0 = 1$. It sets three quadrature nodes, so QMOMKit carries
$m_0, \ldots, m_5$. Both executables draw the same $N = 10000$ samples $V_p$
with a fixed seed ([sampling.hpp](./sampling.hpp)) and give each the weight
$n_0/N$, and the Monte Carlo companion estimates

```math
m_k = \frac{n_0}{N}\sum_{p=1}^{N} V_p^k.
```

`main.cpp` replaces the moments of the configured case with these before it
creates the initial `qmomkit::QmomState`.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_1d_breakage_advection_with_mc_example
./qmom_1d_breakage_advection_with_mc_mc_example
```

## Output

Both executables print the initial and final raw moments $m_0, \ldots, m_5$,
labeled `M0` to `M5`; the Monte Carlo run also prints the initial and final
particle counts. Every $0.1$ in time, the QMOMKit run writes
`build/output/moment_history.txt` and `build/output/quadrature_history.txt`,
and the Monte Carlo run `build/output/moment_history_mc.txt` and
`build/output/particle_history.txt`, which lists the time, identifier and
volume of every particle.

## What to expect

Both runs start from the same sample, so their initial moments agree. At
$t = 1$, the QMOMKit run gives $m_0 = 1.64872$, which is $e^{0.5}$ because
$\mathrm{d}m_0/\mathrm{d}t = 0.5\,m_0$, and the Monte Carlo run
$m_0 = 1.6578$ from 16578 particles; $m_1$ to $m_5$ agree to within 0.2
percent. The Monte Carlo run is a stochastic reference with sampling noise, not
an exact solution, so the two runs are not expected to agree to more digits.
