# 1D QMOM aerosol evaporation

This example, inspired by [Feng et al.
(2023)](https://www.doi.org/10.1016/j.ijmultiphaseflow.2023.104500), evaporates
an aerosol with a six-node QMOM quadrature and works without a state: it calls
`compute_qmom(...)` once and then moves the nodes of a `Quadrature1D` with
`advance_nodes(...)`, as described in [Without a
state](../../docs/user/6_coupling_qmomkit.md#without-a-state). Its phase-space
velocity is written per node and wrapped with `make_phase_space_velocity(...)`.

## Model

The phase-space coordinate is the droplet diameter $d$ in micrometers, time
$t$ is in seconds. The diameters shrink by the $d^2$-law,

```math
\frac{\mathrm{d}(d^2)}{\mathrm{d}t} = -K
\quad\Longleftrightarrow\quad
\frac{\mathrm{d}d}{\mathrm{d}t} = -\frac{K}{2d},
```

with $K = 500\,\mu\mathrm{m}^2/\mathrm{s}$. Evaporation only moves droplets to
smaller diameters, so the weights, the number densities of the nodes, stay
unchanged. Droplets are held at a residual diameter of $0.8\,\mu\mathrm{m}$:
the velocity is zero there, and after each step the driver raises any node
below it back to it. $K$ and the residual diameter are set in
[physics.cpp](./physics.cpp), which defines the velocity and the residual
rule.

[configuration.qmom](./configuration.qmom) defines the initial number density
$n(d)$ as a mixture of three log-normal modes with a total number density
$m_0 = 2.23 \times 10^5$, from which QMOMKit computes the raw moments
$m_0, \ldots, m_{11}$. It sets `adaptive = false` so that all six nodes are
kept. [main.cpp](./main.cpp) advances from $t = 0$ to $20\,\mathrm{s}$ with
steps of $0.2\,\mathrm{s}$.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_1d_aerosol_evaporation_example
```

## Output

The example prints the initial and final total number density $m_0$ and mean
diameter $m_1/m_0$, and the final quadrature. It writes the moments and the
quadrature every $2\,\mathrm{s}$ to `build/output/moments.txt` and
`build/output/quadrature.txt`, as set in the `[output]` section of
`configuration.qmom`.

## What to expect

The number density stays at $2.23 \times 10^5$ and the weights are unchanged.
The mean diameter falls from $12.30$ to $9.08\,\mu\mathrm{m}$. The smallest
node, which carries most of the droplets, starts at $2.61\,\mu\mathrm{m}$ and
sits at the residual diameter from the first output at $2\,\mathrm{s}$ on. The
next node shrinks from $175.6$ to $144.3\,\mu\mathrm{m}$, as the $d^2$-law
predicts, and the larger nodes move less.
