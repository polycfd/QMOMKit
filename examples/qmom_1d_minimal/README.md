# Minimal 1D QMOM evolution

This example is the smallest QMOM time loop: it creates a `QmomState` with
`QmomCase::initial_state()` and advances it with
`advance_qmom_phase_space(...)`, given a `PhaseSpaceVelocity` callback that
moves the quadrature nodes and a `QmomMomentSource` callback that returns
moment rates. Everything is set in [main.cpp](./main.cpp).

## Model

The phase-space coordinate is the particle volume $V$ and $n(V)$ is the number
density, with raw moments $m_k = \int V^k\,n(V)\,\mathrm{d}V$. Volume and time
$t$ are dimensionless. Volumes relax toward $V = 1$,

```math
\frac{\mathrm{d}V}{\mathrm{d}t} = 1 - V,
```

and each particle breaks at the rate $0.5$ into two daughters of half its
volume, $V \to V/2$, which conserves the total volume $m_1$:

```math
\frac{\mathrm{d}m_k}{\mathrm{d}t} = 0.5\left(2^{1-k}-1\right)m_k.
```

The initial moments $m_0, \dots, m_5$ are those of the bimodal mixture
$0.5\,\mathcal{N}(1, 0.1^2) + 0.5\,\mathcal{N}(3, 0.2^2)$ of
[qmom_1d_breakage_advection](../qmom_1d_breakage_advection/README.md), with
three nodes. The run goes from $t = 0$ to $1$ with steps of $0.01$ and the
default options, which use Strang splitting, as described in
[Evolving in time](../../docs/user/5_evolving_in_time.md).

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_1d_minimal_example
```

## Output

The example prints three lines: that the run is complete, and the names of the
two files it wrote. `MomentWriter` and `QuadratureWriter` write the moments and
the quadrature every $0.1$ to `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`; see
[Output files](../../docs/user/7_reference.md#output-files).

## What to expect

Breakage raises $m_0$ as $e^{0.5t}$, from 1 to about 1.65 at $t = 1$. The mean
volume $m_1/m_0$ falls from 2 to about 0.96, and the three abscissas at
$t = 1$ lie at about 0.46, 0.95 and 1.74.
