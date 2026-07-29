# QMOM: 1D continuous evolution with a discontinuous source

## Purpose

This tutorial-style example initializes and evolves a one-dimensional QMOM
representation for a particle-size distribution with continuous node motion
and a discontinuous fragmentation source. The initial state is specified in
`configuration.qmom`. The phase-space coordinate is `xi = d`, where `d` is
particle diameter, and the raw moments are

```text
M_k = integral d^k f(d) dd
```

The current workflow is:

```text
initial distribution
    -> raw moments
    -> QMOM inversion
    -> quadrature nodes and weights
    -> continuous and discontinuous evolution with Strang splitting
```

The executable uses a small bimodal distribution,
`0.5 N(1, 0.1^2) + 0.5 N(3, 0.2^2)`, and reconstructs a three-node QMOM
quadrature from moments `M0` through `M5`. QMOMKit reads the mixture
configuration, generates the analytical raw moments, and then the example calls
`qmomkit::compute_qmom(...)`.

Users can modify the number of QMOM nodes and the mode densities, locations,
and standard deviations directly in `configuration.qmom`. For `n` QMOM nodes,
the configuration must request moments through order `2*n - 1`.

The continuous physics is defined in `physics.cpp` as:

```text
dd/dt = 1 - d
```

Particles smaller than `d = 1` grow, and particles larger than `d = 1` shrink.
QMOMKit advances the quadrature nodes with `qmomkit::advance_nodes(...)`.
Weights are unchanged during this continuous node evolution step.

The discontinuous source is also defined in `physics.cpp`. It is a simple
binary fragmentation model,

```text
d -> d/2
```

with constant rate `0.5`. For raw moments, the source callback returns

```text
dM_k/dt = lambda * (2^(1-k) - 1) * M_k
```

QMOMKit advances those moments with
`qmomkit::advance_moments_with_source(...)`.

The time loop in `main.cpp` uses explicit Strang splitting:

```text
Continuous half-step
        ↓
Moment source update
        ↓
QMOM reconstruction
        ↓
Continuous half-step
```

This keeps the distinction visible: continuous evolution acts on quadrature
nodes, while the fragmentation source acts on raw moments.

These two callbacks are choices made by this breakage model. A different
QMOMKit workflow may use only continuous node advection, only moment-source
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

This example covers initialization, continuous phase-space evolution, a
discontinuous fragmentation source, and explicit Strang splitting. It is one
illustrative callback combination, not a requirement that every transport
workflow include both continuous and source physics.

Requires: QMOMKit.

Required: `configuration.qmom`, `main.cpp`, `physics.cpp`, `utilities.cpp`,
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
./qmom_1d_continuous_breakage_example
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
