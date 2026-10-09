# QMOM: 1D phase-space advection with a moment-source

## Purpose

This tutorial-style example initializes and evolves a one-dimensional QMOM
representation for a particle-size distribution with advection node motion
and a size-dependent fragmentation source written over the quadrature nodes. The initial state is specified in
[configuration.qmom](./configuration.qmom). The phase-space coordinate is $\xi=d$, where $d$ is
particle diameter and $n(d)$ is the number density. The raw, unnormalized
moments are
```math
M_k = \int_{\mathbb{R}} d^k n(d)\,\mathrm{d}d.
```

The workflow is:

```text
initial distribution
    -> raw moments
    -> QMOM inversion
    -> quadrature nodes and weights
    -> phase-space advection and moment-source evolution with Strang splitting
```

The executable uses a small bimodal distribution,
$`0.5\,\mathcal{N}(1, 0.1^2) + 0.5\,\mathcal{N}(3, 0.2^2)`$, and reconstructs a three-node QMOM
quadrature from moments $M_0$ through $M_5$. QMOMKit reads the mixture
configuration, generates the analytical raw moments, and then the example calls
`qmomkit::QmomCase::initial_state()` to construct the initial `QmomState`.

Users can modify the number of QMOM nodes and the mode densities, locations,
and standard deviations directly in [configuration.qmom](./configuration.qmom).
For $n$ QMOM nodes, the configuration must request moments through order $2N -
1$.

The phase-space evolution defined in [physics.cpp](./physics.cpp) is:

$$\frac{\mathrm{d}d}{\mathrm{d}t} = 1-d. $$

Particles smaller than $d=1$ grow, and particles larger than $d=1$ shrink. The
example passes this callback to `qmomkit::advance_qmom_phase_space(...)`, which
manages the phase-space evolution as part of the configured phase-space
evolution step. The underlying advection stage advances quadrature nodes, and
weights are unchanged during that stage.

The moment-source is also defined in [physics.cpp](./physics.cpp). It is a
binary fragmentation model,

$$d \to \frac{d}{2},$$

in which one particle of diameter $d$ becomes two of diameter $d/2$, at a rate
that grows with size,

$$a(d) = a_\mathrm{ref}\left(\frac{d}{d_\mathrm{ref}}\right)^{p},
\qquad a_\mathrm{ref} = 0.5,\quad d_\mathrm{ref} = 1,\quad p = \tfrac{3}{2}.$$

The raw moments then change as

```math
\frac{\mathrm{d}M_k}{\mathrm{d}t}
= \left(2^{1-k}-1\right)\int a(d)\,d^k\,n(d)\,\mathrm{d}d
= \left(2^{1-k}-1\right)\frac{a_\mathrm{ref}}{d_\mathrm{ref}^{\,p}}\,M_{k+3/2}.
```

The right-hand side needs moments of order $k + 3/2$, which are not among the
transported moments $M_0$ through $M_5$. The source callback closes it with the
quadrature: it evaluates the integral over the nodes of the quadrature that
QMOMKit builds from the current moments, through `context.for_each_node(...)`:

```math
\frac{\mathrm{d}M_k}{\mathrm{d}t} \approx \left(2^{1-k}-1\right)
\sum_i w_i\,a(d_i)\,d_i^k.
```

QMOMKit reconstructs that quadrature at every Runge-Kutta stage of the
moment-source update, so the node values follow the moments as they change. The
example passes this callback to `qmomkit::advance_qmom_phase_space(...)`, which
manages the moment-source update as part of the configured phase-space
evolution step. With $p = 0$ the rate is constant and the source reduces to
$\mathrm{d}M_k/\mathrm{d}t = a_\mathrm{ref}\left(2^{1-k}-1\right)M_k$, the
model used by
[`qmom_1d_breakage_advection_with-mc`](../qmom_1d_breakage_advection_with-mc/README.md).

The time loop in [main.cpp](./main.cpp) calls
`qmomkit::advance_qmom_phase_space(...)` for each configured phase-space
evolution step. With the shipped `strang` setting, the native phase-space
evolution operation uses:

```text
Advection half-step
        ↓
Moment-source update
        ↓
QMOM reconstruction
        ↓
Advection half-step
```

This keeps the distinction visible: phase-space advection moves the quadrature
nodes, while the fragmentation source changes the raw moments, using the nodes
of each stage quadrature only to evaluate its rates.

These two callbacks are choices made by this breakage model. A different
QMOMKit workflow may use only phase-space advection, only moment-source
evolution, or add a boundary contribution callback when the physics includes
phase-space boundary exchange.

The time controls are local constants in [main.cpp](./main.cpp):

- `initial_time`
- `final_time`
- `time_step`

These define the physical simulation duration and timestep size. The output
cadence is read from `[output].interval` in
[configuration.qmom](./configuration.qmom).

The quadrature weights represent the amount of distribution mass carried by
each node. The abscissas are the particle diameters at which that mass is
concentrated in the QMOM approximation.

## Build and run

Build from this example directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies [configuration.qmom](./configuration.qmom) into `build/`. Run
from the build directory so generated output stays there:

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
the standard QMOMKit output section in
[configuration.qmom](./configuration.qmom):

```ini
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
grows. Because large particles break faster, the large mode near $d = 3$
initially fragments about five times as fast as the small mode near $d = 1$.
Fragmentation and advection both shrink the large particles, so the growth of
$M_0$ slows over time: it gains about 0.14 in the first 0.1 time units and
about 0.09 in the last. By $t = 1$, $M_0$ has grown from 1 to about 2.04,
compared with about 1.65 at the constant rate $a_\mathrm{ref}$, while $M_5$
has fallen from about 127 to about 2.8, and the third node carries the least
weight. The final quadrature shows the combined effect of population
multiplication and diameter relaxation rather than either mechanism in
isolation.
