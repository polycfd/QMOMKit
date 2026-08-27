# QMOM: 1D MC initialization and sample-based evolution

## Purpose

This tutorial-style example is the Monte Carlo sampled-initialization companion
to `examples/qmom_1d_breakage_advection/`.

This is an optional reference workflow. It is not required for ordinary QMOM
phase-space evolution onboarding.

Both examples use the same configured bimodal distribution,
$0.5\,\mathcal{N}(1, 0.1^2) + 0.5\,\mathcal{N}(3, 0.2^2)$. This directory contains two workflows
that start from the same finite sampled ensemble. The Monte Carlo companion
tool is used in both of its main roles here: as a sample realization tool and
as a direct sampled-ensemble reference path.

Sample realization workflow:

```text
distribution
    -> fixed-seed Monte Carlo sample
    -> finite sample ensemble
```

The Monte Carlo tool represents a finite sample ensemble. In this example, each
sample is interpreted as a physical particle for the breakage model.

Analytical initialization example:

```text
distribution
    -> analytical moments
    -> QMOM
```

MC initialization example:

```text
distribution
    -> sampled ensemble
    -> statistical moments
    -> QMOM
    -> phase-space advection + moment-source QMOM evolution
```

Monte Carlo sampled-distribution companion:

```text
distribution
    -> same sampled ensemble
    -> sample-based stochastic evolution
    -> statistical moments
```

The distribution remains defined in `configuration.qmom`. This example draws
10,000 samples with a fixed seed. With `d` as the particle-diameter coordinate,
the QMOM workflow estimates raw, unnormalized moments from the ensemble as

$$
M_k = \frac{n_0}{N}\sum_{p=1}^{N} d_p^k,
$$

where $n_0$ is the configured total number density and $N$ is the sample count.
It uses those moments in `qmomkit::compute_qmom(...)`. The Monte Carlo
workflow evolves the sampled ensemble directly with the same advection and
binary-fragmentation models.

The sample is finite, so its raw moments include ordinary Monte Carlo sampling
uncertainty. Increasing the sample count generally reduces this uncertainty,
although high-order moments and broad distributions may converge more slowly.
The initial total number density is distributed across the samples as a
constant statistical weight.

The Monte Carlo workflow is a direct stochastic particle-based
companion/reference realization, not an automated formal validation test. It
does not estimate uncertainty or declare one history to be exact truth.

## QMOM evolution

After sampled initialization, the QMOM evolution matches the analytical
initialization example:

```text
Phase-space advection half-step
        ↓
Moment-source update
        ↓
QMOM reconstruction
        ↓
Phase-space advection half-step
```

The phase-space advection physics and moment-source binary fragmentation source are in
`physics.cpp`. The time controls are local constants in `main.cpp`:

- `initial_time`
- `final_time`
- `time_step`

The output directory, enabled files, and output cadence are read from
`[output]` in `configuration.qmom`.

## Monte Carlo Moment Evolution

The Monte Carlo workflow uses the same configured distribution, sample count, and
initialization seed as the QMOM workflow. It then applies the same Strang
sequence directly to the sampled ensemble:

```text
Advection sampled-ensemble half-step
        ↓
Stochastic binary fragmentation
        ↓
Advection sampled-ensemble half-step
```

The advection sampled-ensemble update uses
`qmomkit::monte_carlo::evolve_ensemble`.
The stochastic breakage mechanics are local to this example in
`particles_mc.cpp`. The advection velocity and fragmentation rate are shared
with the QMOM workflow through `physics.cpp`.

## Build

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
cd examples/qmom_1d_breakage_advection_with-mc

cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`.

## Run the QMOM example

Run from the build directory so generated output stays there:

```bash
cd build
./qmom_1d_breakage_advection_with_mc_example
```

Running the executable writes:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

These files are written with `qmomkit::MomentWriter` and
`qmomkit::QuadratureWriter`. They are whitespace-delimited text tables with a
commented header line.

## Run the Direct Monte Carlo Companion

From the same build directory, run:

```bash
./qmom_1d_breakage_advection_with_mc_mc_example
```

Running the executable writes:

- `output/moment_history_mc.txt`
- `output/particle_history.txt`

The Monte Carlo moments use a constant per-sample statistical weight from the
initial sampled ensemble and are written with `qmomkit::MomentWriter`.
`particle_history.txt` is a custom particle snapshot table, not a QMOMKit
moment or quadrature output.

## What to expect

Both executables reproduce the same fixed-seed sampled initialization, so their
initial sampled moments describe the same finite ensemble. The direct Monte
Carlo companion evolves particles stochastically, whereas QMOM evolves a finite
quadrature closure of sampled moments. Compare trends such as number-density
growth and moment evolution, but do not expect final histories to agree
point-for-point: finite sampling variation and the different representations
produce differences during evolution.
