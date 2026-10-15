# 2D CQMOM breakage and advection with a Monte Carlo reference

This example is the two-dimensional counterpart of
[qmom_1d_breakage_advection_with-mc](../qmom_1d_breakage_advection_with-mc/README.md):
particles that break and relax in volume and velocity, run twice from the same
finite sample. `main.cpp` loads a `qmomkit::CqmomPreparedCase`, replaces its
initial moments with those of the sample through
`qmomkit::config::with_initial_moments(...)` and advances them with
`qmomkit::advance_cqmom_phase_space(...)`; `main_mc.cpp`, the Monte Carlo
reference, advances the sampled particles themselves with the Monte Carlo
companion.

## Model

The phase-space coordinates are, in this order, the particle volume $V$ and the
particle velocity $u$, and $n(V,u)$ is the number density; $V$, $u$, time $t$
and all parameters are dimensionless. The raw moments are
$M_{i,j} = \iint V^i u^j\, n(V,u)\,\mathrm{d}V\,\mathrm{d}u$. Volumes relax
toward $1$ and velocities toward $0$,

```math
\frac{\mathrm{d}V}{\mathrm{d}t} = 1 - V,
\qquad
\frac{\mathrm{d}u}{\mathrm{d}t} = -u,
```

and particles break at the constant rate $\lambda = 0.5$, one particle of
volume $V$ into two of volume $V/2$ with the same velocity, which conserves
the total volume $M_{1,0}$. The moments change as

```math
\frac{\mathrm{d}M_{i,j}}{\mathrm{d}t} = \lambda\left(2^{1-i}-1\right)M_{i,j}.
```

In the Monte Carlo run each particle breaks within a step $\Delta t$ with
probability $1-\exp(-\lambda\,\Delta t)$. Both runs use Strang splitting: half
a step of advection, the breakage over the whole step, and a second half step
of advection. The velocity and the moment source are in
[physics.cpp](./physics.cpp), $\lambda$ in [physics.hpp](./physics.hpp), the
stochastic breakage in [particles_mc.cpp](./particles_mc.cpp), and the time
span $0 \le t \le 1$ with $\Delta t = 0.01$ in both drivers.

[configuration.qmom](./configuration.qmom) defines the initial number density
as the total number density $n_0 = 1$ times the product of the normal
distributions $V \sim \mathcal{N}(2, 0.5^2)$ and $u \sim \mathcal{N}(0, 0.5^2)$.
It sets `node_counts = 2, 2`, at most two nodes in $V$ and, conditioned on
each, two in $u$, with the default adaptive reconstruction, and
`maximum_orders = 3, 3`, so the sample provides $M_{i,j}$ for $i, j \le 3$.
Both executables draw the same $N = 10000$ samples $(V_p, u_p)$ with a fixed
seed ([sampling.hpp](./sampling.hpp)) and give each the weight $n_0/N$, and
the Monte Carlo companion estimates

```math
M_{i,j} = \frac{n_0}{N}\sum_{p=1}^{N} V_p^i u_p^j.
```

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_2d_breakage_advection_with_mc_example
./cqmom_2d_breakage_advection_with_mc_mc_example
```

## Output

Both executables print the initial and final raw moments $M_{0,0}$,
$M_{1,0}$, $M_{0,1}$, $M_{2,0}$, $M_{1,1}$ and $M_{0,2}$; the Monte Carlo run
also prints the initial and final particle counts. Every $0.1$ in time, the
QMOMKit run writes `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`, and the Monte Carlo run
`build/output/moment_history_mc.txt`.

## What to expect

Both runs start from the same sample, so their initial moments agree, and the
QMOMKit run keeps all four nodes. At $t = 1$, the QMOMKit run gives
$M_{0,0} = 1.64872$, which is $e^{0.5}$ because
$\mathrm{d}M_{0,0}/\mathrm{d}t = 0.5\,M_{0,0}$, and the Monte Carlo run
$M_{0,0} = 1.6578$ from 16578 particles. $M_{1,0}$ and $M_{2,0}$ agree to
within 0.2 percent and $M_{0,2}$ to within 1.3 percent; $M_{0,1}$ and $M_{1,1}$
stay near zero in both runs, where sampling noise dominates their difference.
The Monte Carlo run is a stochastic reference with sampling noise, not an exact
solution.
