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
homogeneous and uses the same one-dimensional continuous-breakage case as
`examples/qmom_1d_continuous_breakage_mc_initialization/`: the same configured
bimodal distribution, fixed-seed Monte Carlo initialization, continuous
velocity, binary-fragmentation source, timestep, output cadence, and output
files.

## Difference From The Native Workflow

The native-workflow example starts from sampled raw moments, reconstructs a
QMOM quadrature, advances quadrature nodes for the continuous process, advances
raw moments for the discontinuous source, and reconstructs updated moments
from the evolved quadrature.

This example keeps the raw moments as application-owned state. Its RHS
callbacks call:

```cpp
qmomkit::evaluate_qmom_continuous_moment_rates(...)
qmomkit::evaluate_qmom_discontinuous_moment_rates(...)
```

inside application-managed SSP-RK3 steps. These functions return
instantaneous raw-moment rates; they do not advance time. The example owns the
SSP-RK3 stages and the continuous/source sequence.

Because the continuous updates are discrete in different variables, the two
examples are not expected to be bitwise-identical. The native workflow moves
quadrature nodes directly. This example advances raw moments using rates
reconstructed from the quadrature at each SSP-RK3 RHS evaluation. Both
represent the same underlying phase-space process.

## Numerical Workflow

After Monte Carlo initialization, each macro step uses the same
continuous/source ordering as the native example:

```text
continuous half-step
    -> discontinuous source step
    -> continuous half-step
```

For the continuous half-steps,
`evaluate_qmom_continuous_moment_rates(...)` inverts the supplied raw moments,
evaluates the phase-space velocity at QMOM nodes, and reconstructs
instantaneous raw-moment rates. For the discontinuous source step,
`evaluate_qmom_discontinuous_moment_rates(...)` calls the source model, which
already returns raw-moment rates; that rate evaluation does not require a
quadrature inversion.

Moments remain the evolved state. The driver calls `qmomkit::compute_qmom(...)`
only when it needs weights and abscissas for output.

## Stage Policy

The public rate evaluators do not define SSP-RK3 stage acceptability. This
example explicitly chooses successful `qmomkit::compute_qmom(...)` inversion
as its stage and final acceptance policy so later continuous-rate calls receive
invertible moment states. That is one application policy, not a universal
definition of acceptability. A finite-element application may add bounds,
limiters, or other coefficient-level criteria.

## Source Organization

`main.cpp` contains the application-managed rate evaluation, SSP-RK3
integration, stage policy, and simulation loop. `utilities.*` contains only
output, filesystem, and reporting support.

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
cd examples/qmom_1d_continuous_breakage_mc_local_rates

cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`. Run from the build
directory so generated output stays there:

```bash
cd build
./qmom_1d_continuous_breakage_mc_local_rates_example
```

The helper script follows the same configure, build, and run pattern:

```bash
examples/qmom_1d_continuous_breakage_mc_local_rates/build_and_run.sh "$PWD/install"
```

## Limitations

This is not a CG or DG solver. It does not include a mesh, basis functions,
spatial quadrature, residual assembly, a mass matrix, physical-space fluxes,
or spatial realizability limiting. Its purpose is to show the local QMOMKit
rate call that a finite-element code would place inside its residual
evaluation.
