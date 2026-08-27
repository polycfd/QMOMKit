# QMOM: 1D local moment-rate integration

## Purpose

This example demonstrates the local moment-rate interface for CG and DG
finite-element coupling. In that setting, the host code advances spatial
degrees of freedom or element coefficients. At each spatial quadrature point
and Runge-Kutta stage, it may reconstruct a local raw-moment vector, ask
QMOMKit for the corresponding phase-space moment rate, project that rate back
onto the finite-element basis, and advance the coefficient-level ODE with the
host solver's time integrator.

This example contains only the QMOMKit-facing local operation. It is spatially
homogeneous and uses the same one-dimensional advection-breakage case as
`examples/qmom_1d_breakage_advection_with-mc/`: the same configured
bimodal distribution, fixed-seed Monte Carlo initialization, advection
velocity, binary-fragmentation source, timestep, output cadence, and output
files. It also includes a separate Monte Carlo reference executable that
evolves the sampled particles directly with the same physical model, including
a simple nucleation boundary contribution.

## Difference From The Native Workflow

The native-workflow example starts from sampled raw moments, reconstructs a
QMOM quadrature, advances quadrature nodes for the advection process, advances
raw moments for the moment-source, and reconstructs updated moments
from the evolved quadrature.

This example keeps the raw moments as application-owned state. Its RHS
callbacks call:

```cpp
qmomkit::evaluate_qmom_advection_moment_rates(...)
qmomkit::evaluate_qmom_moment_source_rates(..., source, boundary)
```

inside application-managed SSP-RK3 steps. These functions return
instantaneous raw-moment rates; they do not advance time. The example owns the
SSP-RK3 stages and the advection/source sequence.

Because the advection updates are discrete in different variables, the two
examples are not expected to be bitwise-identical. The native workflow moves
quadrature nodes directly. This example advances raw moments using rates
reconstructed from the quadrature at each SSP-RK3 RHS evaluation. Both
represent the same underlying phase-space process.

## Numerical Workflow

After Monte Carlo initialization, each macro step uses the same
advection/source ordering as the native example:

```text
advection half-step
    -> moment-source stage
    -> advection half-step
```

The executable loads the configured phase-space evolution options and supports the default
`strang` sequence used here. Other configured phase-space evolution schemes are rejected
explicitly by this local-rate driver.

For the advection half-steps,
`evaluate_qmom_advection_moment_rates(...)` inverts the supplied raw moments,
evaluates the phase-space velocity at QMOM nodes, and reconstructs
instantaneous raw-moment rates. For the moment-source stage,
`evaluate_qmom_moment_source_rates(...)` calls the source model, which
already returns raw-moment rates; that rate evaluation does not require a
quadrature inversion. The moment-source rate is combined with a lower-boundary
nucleation contribution by the four-argument QMOM moment-source rate evaluator,
which returns the complete instantaneous rate to the application.

The moment convention uses raw, unnormalized number-density moments:

$$
m_k = \int_{\mathbb{R}} d^k n(d)\,\mathrm{d}d.
$$

So constant-rate nucleation at diameter $d_{\mathrm{nuc}}$ contributes

$$
B_k = J_{\mathrm{nuc}}d_{\mathrm{nuc}}^k
$$

to raw moment order $k$. This example uses
$J_{\mathrm{nuc}} = 0.05$ and $d_{\mathrm{nuc}} = 0.2$. The rate is visible in the histories over
the unit simulation interval, but it remains small relative to the existing
growth and breakage dynamics.

Moments remain the evolved state. The driver calls `qmomkit::compute_qmom(...)`
only when it needs weights and abscissas for output.

## Monte Carlo Reference

`main_mc.cpp` starts from the same configured distribution, sample count, and
initialization seed as `main.cpp`. It interprets those samples as particles,
uses the shared advection velocity, applies the shared binary breakage event
rate and daughter-state rule, and injects nucleated particles at the shared
nucleation diameter.

The Monte Carlo reference follows the same macro-step sequence:

```text
advection particle half-step
    -> stochastic binary breakage and Poisson nucleation
    -> advection particle half-step
```

Each computational particle has the fixed statistical weight
$w_p = n_0/N$. Raw moments are weighted sums over particles; they are not
divided by the current particle count, so increases in the raw zeroth moment
$m_0$ from breakage and nucleation remain physical number-density changes.
Over a source interval $\Delta t$, the Monte Carlo nucleation count is sampled
from a Poisson distribution with mean
$J_{\mathrm{nuc}}\Delta t / w_p$, where $w_p$ is the particle weight. The
expected weighted birth contribution is the same
$J_{\mathrm{nuc}}\Delta t$ used by the QMOM boundary rate. The Monte Carlo result is stochastic after initialization because
breakage and nucleation events are sampled with fixed deterministic seeds. The
fixed seeds make the example reproducible, but the final moments are not
expected to match the QMOM local-rate result exactly.

## Stage Policy

The public rate evaluators do not define SSP-RK3 stage acceptability. This
example explicitly chooses successful `qmomkit::compute_qmom(...)` inversion
as its stage and final acceptance policy so later advection-rate calls receive
invertible moment states. That is one application policy, not a universal
definition of acceptability. A finite-element application may add bounds,
limiters, or other coefficient-level criteria.

## Source Organization

`physics.*` contains the shared physical model: advection velocity, binary
breakage rate, daughter state, nucleation rate, nucleation diameter, and the
QMOM nucleation moment contribution. `local_rate_workflow.cpp` contains the
application-managed rate evaluation, SSP-RK3 integration, stage policy, and
QMOM simulation loop; `main.cpp` only dispatches to that workflow.
`main_mc.cpp` contains the Monte Carlo particle workflow.
`particles_mc.*` adapts the shared physical model to the sampled-particle
representation. `nucleation_balance_check.cpp` verifies the analytical
correspondence between the QMOM nucleation moment rate and the MC birth
intensity. `utilities.*` contains only output, filesystem, and reporting
support.

## Output

The example writes the same QMOM moment and quadrature history files as the
native-workflow example:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

The files are written with QMOMKit's public `qmomkit::MomentWriter` and
`qmomkit::QuadratureWriter`. They are whitespace-delimited text tables with
commented headers.

The Monte Carlo executable writes:

- `output/moment_history_mc.txt`
- `output/mc_diagnostics.txt`
- `output/particle_history.txt`

`moment_history_mc.txt` uses the same time column and raw-moment ordering as
`moment_history.txt`. `mc_diagnostics.txt` records particle count, total number
density, first raw moment, cumulative breakage events, cumulative nucleation
particles, expected cumulative nucleation number, realized cumulative
nucleation number, and the nucleation contribution to $m_0$.
`particle_history.txt` is a custom sampled-particle snapshot table.

## Build And Run

Requires:

- QMOMKit
- QMOMKitMonteCarlo

From the repository root, install both packages into the shared
repository-local `install/` prefix:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build build
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

Then configure this example from its directory:

```bash
cd examples/qmom_1d_breakage_advection_local-rates_with-mc

cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`. Run from the build
directory so generated output stays there. Run the QMOM local-rate executable
first when you want the Monte Carlo executable to print direct comparison
metrics against `output/moment_history.txt`:

```bash
cd build
./qmom_1d_breakage_advection_local_rates_with_mc_example
./qmom_1d_breakage_advection_local_rates_with_mc_mc_example
```

The build also produces a standalone nucleation-balance executable. Run it
directly when you want to check the analytical correspondence between the QMOM
boundary rate and the Monte Carlo birth intensity:

```bash
./qmom_1d_breakage_advection_local_rates_with_mc_nucleation_balance_check
```

## What to expect

The QMOM history is produced by application-owned SSP-RK3 integration of the
local rates; QMOMKit supplies those rates but does not own the moment state or
time integrator. The source history includes binary-breakage growth and the
small lower-boundary nucleation contribution. The balance-check executable
verifies the analytical correspondence between that QMOM boundary rate and the
Monte Carlo birth intensity. The direct Monte Carlo companion is a stochastic
particle-based companion/reference realization, so its later moments need not
match the local-rate history point-for-point.

## Limitations

This is not a CG or DG solver. It does not include a mesh, basis functions,
spatial quadrature, residual assembly, a mass matrix, physical-space fluxes,
or spatial realizability limiting. Its purpose is to show the local QMOMKit
rate call that a finite-element code would place inside its residual
evaluation.
