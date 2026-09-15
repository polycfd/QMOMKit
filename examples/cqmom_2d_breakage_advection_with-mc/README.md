# CQMOM: 2D MC Initialization with Phase-Space Advection and Moment-Source Evolution

## Purpose

This example is the Monte Carlo sampled-initialization companion to
`examples/cqmom_2d_minimal/`.
It contains two workflows that start from the same fixed-seed sampled ensemble.
The first converts the ensemble into sampled raw moments, initializes a CQMOM
case through the configuration API, and runs CQMOM. The second evolves the
sampled ensemble directly. This demonstrates the Monte Carlo companion tool in
both of its main roles: drawing a finite sample realization and producing an
independent sampled-ensemble reference path.

Sample realization workflow:

```text
configuration
        |
        v
fixed-seed MC sampling
        |
        v
finite sample ensemble
```

The Monte Carlo tool represents a finite sample ensemble. In this example, each
sample is interpreted as a physical particle for the breakage model.

CQMOM workflow:

```text
configuration
        |
        v
fixed-seed MC sampling
        |
        v
ensemble moments
        |
        v
CQMOM case
        |
        v
compute_cqmom(...) -> initial quadrature
        |
        v
CqmomCase::initial_state() -> CqmomState
        |
        v
CQMOM phase-space evolution
```

The CQMOM driver loads this numerical setup as a normal `CqmomCase`, estimates
the complete sampled tensor, and creates an immutable replacement case with
`with_initial_moments(...)`. The hierarchy, storage selection, and
reconstruction policies therefore remain those declared by
`configuration.qmom`; sampling is only one possible producer of replacement raw
moments.

Monte Carlo sampled-distribution companion:

```text
configuration
        |
        v
fixed-seed MC sampling
        |
        v
sampled ensemble
        |
        v
phase-space advection + stochastic moment-source evolution
```

Both workflows use the same `configuration.qmom`, distribution, sample count,
initialization seed, phase-space ordering, advection velocity, binary
fragmentation model, and Strang splitting.

The phase-space ordering is $(d,u)$,

where `d` is particle size and `u` is velocity. The same ordering is used in
`configuration.qmom`, the sampled ensemble, the CQMOM case, and the output.

## Initialization

Edit `configuration.qmom` to change the independent-product distribution and
CQMOM hierarchy. The default configuration uses:

$$
\begin{aligned}
d &\sim \mathcal{N}(2, 0.5^2), \\
u &\sim \mathcal{N}(0, 0.5^2).
\end{aligned}
$$

This example draws 10,000 samples with a fixed initialization seed. QMOMKit
computes raw moments from that finite ensemble using the absolute
number-density convention:

$$
M_{i,j} = \frac{n_0}{N}\sum_{p=1}^{N}d_p^i u_p^j.
$$

Equivalently, each initial sample has statistical weight $n_0/N$ in the
moment estimate. The finite ensemble has ordinary Monte Carlo sampling
uncertainty. Increasing the sample count generally reduces that uncertainty,
but high-order tensor moments can require substantially more samples.

The CQMOM driver treats the sampled raw moment tensor as an initialization
object and asks the configuration API for a `qmomkit::CqmomCase` consumed by
`qmomkit::compute_cqmom(...)` and `qmomkit::CqmomCase::initial_state()`. The
direct MC driver uses the sampled ensemble itself. Because the two drivers are
separate executables, they reproduce the same initial ensemble by using the
same configuration, sample count, seed, and sampling utility.

The default CQMOM hierarchy is `node_counts = 2, 2`, meaning:

- two quadrature nodes in `d`;
- two conditional quadrature nodes in `u`;
- four terminal CQMOM nodes.

The CQMOM driver also supports adaptive branch counts through the same
`CqmomState` phase-space evolution workflow. With adaptive branch counts
enabled, these node counts are maximums and output may contain fewer active
terminal nodes.

The sampled moment tensor contains $M_{i,j}$ for $i=0,\ldots,3$ and
$j=0,\ldots,3$.

## Advection velocity

The shared velocity is in `physics.cpp`:

$$
\begin{aligned}
\frac{\mathrm{d}d}{\mathrm{d}t} &= 1-d, \\
\frac{\mathrm{d}u}{\mathrm{d}t} &= -u.
\end{aligned}
$$

Diameter relaxes toward `1`, and velocity relaxes toward `0`. The CQMOM driver
uses `qmomkit::advance_cqmom_phase_space(...)` with `qmomkit::CqmomState`. The
direct MC driver advances sampled ensemble coordinates with
`qmomkit::monte_carlo::evolve_ensemble(...)`.

## Moment source

The shared fragmentation model is:

$$
d \to \frac{d}{2}.
$$

Velocity is unchanged by fragmentation. In the CQMOM driver, the source acts
diagonally on each CQMOM moment:

$$
\frac{\mathrm{d}M_{i,j}}{\mathrm{d}t}
= \lambda\left(2^{1-i}-1\right)M_{i,j}.
$$

In the direct MC driver, each sample breaks stochastically over the source
step with probability

$$
1-\exp(-\lambda\,\Delta t).
$$

When a sample breaks, the parent is replaced by two daughters with half the
parent diameter and the same velocity.

Both workflows use the same Strang sequence:

```text
advection half-step
        |
        v
moment-source or stochastic breakage
        |
        v
advection half-step
```

The time controls are local constants in each driver:

- `initial_time = 0.0`
- `final_time = 1.0`
- `time_step = 0.01`

The output directory, enabled files, and output cadence are read from
`[output]` in `configuration.qmom` through `qmomkit::OutputOptions`.

## Build and run

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
cd examples/cqmom_2d_breakage_advection_with-mc

cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`.

## Run the CQMOM workflow

Run from the build directory so generated output stays there:

```bash
cd build
./cqmom_2d_breakage_advection_with_mc_example
```

Running the CQMOM executable writes:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

It also prints the sampled initialization metadata, selected initial and final
raw moments, and the final terminal quadrature nodes and weights.

`initial_moments.txt` records raw moments reconstructed from the initial CQMOM
quadrature built from the sampled initialization. `moment_history.txt` records
the deterministic CQMOM raw-moment evolution. `initial_quadrature.txt` and
`quadrature_history.txt` record terminal node weights and abscissas for
inspecting the reconstructed representation. All four files are written through
QMOMKit-native writers as whitespace-delimited text tables with commented
headers.

## Run the direct Monte Carlo workflow

From the same build directory, run:

```bash
./cqmom_2d_breakage_advection_with_mc_mc_example
```

Running the direct MC executable writes:

- `output/moment_history_mc.txt`

It also prints the sampled initialization metadata, selected initial and final
raw moments, and the final sample count.

`moment_history_mc.txt` stores Monte Carlo estimated moments with
`qmomkit::MomentWriter`. It is a whitespace-delimited text table with a
commented header line. It uses a constant per-sample statistical weight from
the initial sampled ensemble, so $M_{0,0}$ reflects the sample population
increase caused by fragmentation.
