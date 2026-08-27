# QMOM: 1D phase-space advection with a moment-source

## Purpose

This tutorial-style example initializes and evolves a one-dimensional QMOM
representation for a particle-size distribution with advection node motion
and a moment-source fragmentation source. The initial state is specified in
`configuration.qmom`. The phase-space coordinate is $\xi=d$, where `d` is
particle diameter and $f(d)$ is the number density. The raw, unnormalized
moments are

$$
M_k = \int_{\mathbb{R}} d^k f(d)\,\mathrm{d}d.
$$

The current workflow is:

```text
initial distribution
    -> raw moments
    -> QMOM inversion
    -> quadrature nodes and weights
    -> phase-space advection and moment-source evolution with Strang splitting
```

The executable uses a small bimodal distribution,
$0.5\,\mathcal{N}(1, 0.1^2) + 0.5\,\mathcal{N}(3, 0.2^2)$, and reconstructs a three-node QMOM
quadrature from moments $M_0$ through $M_5$. QMOMKit reads the mixture
configuration, generates the analytical raw moments, and then the example calls
`qmomkit::compute_qmom(...)`.

Users can modify the number of QMOM nodes and the mode densities, locations,
and standard deviations directly in `configuration.qmom`. For $n$ QMOM nodes,
the configuration must request moments through order `2*n - 1`.

The phase-space evolution defined in `physics.cpp` is:

$$
\frac{\mathrm{d}d}{\mathrm{d}t} = 1-d.
$$

Particles smaller than $d=1$ grow, and particles larger than $d=1$ shrink.
The example passes this callback to `qmomkit::advance_qmom_phase_space(...)`,
which manages the phase-space evolution as part of the configured
phase-space evolution step. The underlying advection stage advances quadrature nodes, and
weights are unchanged during that stage.

The moment-source is also defined in `physics.cpp`. It is a simple
binary fragmentation model,

$$
d \to \frac{d}{2},
$$

with constant rate $\lambda=0.5$. For raw moments, the source callback
returns

$$
\frac{\mathrm{d}M_k}{\mathrm{d}t}
= \lambda\left(2^{1-k}-1\right)M_k.
$$

The example passes this callback to `qmomkit::advance_qmom_phase_space(...)`,
which manages the moment-source update as part of the configured
phase-space evolution step.

The time loop in `main.cpp` calls `qmomkit::advance_qmom_phase_space(...)` for
each configured phase-space evolution step. With the shipped `strang` setting, the native
phase-space evolution operation uses:

```text
Advection half-step
        ↓
Moment-source update
        ↓
QMOM reconstruction
        ↓
Advection half-step
```

This keeps the distinction visible: phase-space advection acts on quadrature
nodes, while the fragmentation source acts on raw moments.

These two callbacks are choices made by this breakage model. A different
QMOMKit workflow may use only phase-space advection, only moment-source
evolution, or add a boundary contribution callback when the physics includes
phase-space boundary exchange.

The time controls are local constants in `main.cpp`:

- `initial_time`
- `final_time`
- `time_step`

These define the physical simulation duration and timestep size. The output
cadence is read from `[output].interval` in `configuration.qmom`.

The quadrature weights represent the amount of distribution mass carried by
each node. The abscissas are the particle diameters at which that mass is
concentrated in the QMOM approximation.

## Run the example

This example covers initialization, phase-space evolution, a
moment-source fragmentation source, and Strang splitting. It is one
illustrative callback combination, not a requirement that every phase-space evolution
workflow include both advection and source physics.

Requires: QMOMKit.

Required: `configuration.qmom`, `main.cpp`, and `physics.cpp`,
and an installed QMOMKit package. Optional: output-file analysis after the run.

Build from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`. Run from the build
directory so generated output stays there:

```bash
cd build
./qmom_1d_breakage_advection_example
```

Running the executable writes:

- `output/moments.txt`
- `output/quadrature.txt`

It also prints a concise console summary with the initial and final raw
moments, plus initial and final quadrature nodes and weights.

The output directory, enabled files, and write cadence are configured through
the standard QMOMKit output section in `configuration.qmom`:

```text
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
```

`moments.txt` records the raw-moment evolution of the represented
distribution. `quadrature.txt` records the reconstructed node weights and
abscissas at the same output times so users can inspect the numerical
representation. Both files are written through QMOMKit-native writers as
whitespace-delimited text tables with commented headers.

## What to expect

The fragmentation source increases $M_0$, so the represented number density
grows. At the same time, advection moves node diameters toward `d = 1`.
Consequently, the higher raw moments and the final quadrature show the combined
effect of population multiplication and diameter relaxation rather than either
mechanism in isolation.
