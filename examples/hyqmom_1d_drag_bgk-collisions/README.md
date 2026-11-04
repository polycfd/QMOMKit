# 1D HyQMOM drag and BGK collisions

This example evolves the velocity distribution of particles with a `QmomState`
of the `hyperbolic` family, HyQMOM [[Fox & Laurent,
2022](https://www.doi.org/10.1137/21M1406143)]. With node count 2, its
quadrature has three nodes, reconstructed from $m_0, \ldots, m_4$.
`advance_qmom_phase_space(...)` combines two callbacks: a velocity, the drag,
which moves the abscissas, and a source, the collisions, which change the
moments. At the start and the end, the example prints $m_5$ of the quadrature,
which HyQMOM provides beyond the moments of the state, and the characteristic
speeds, which `compute_hyqmom_characteristic_speeds(...)` returns. A code that
transports $m_0, \ldots, m_4$ in physical space needs both: the transport of
$m_4$ involves $m_5$, and the slowest and fastest characteristic speeds are
those at which changes in the moments travel, as
[hyqmom_characteristic_speeds](../hyqmom_characteristic_speeds/README.md)
explains.

## Model

The phase-space coordinate is the particle velocity $u$ in m/s, named `u` in
the configuration, and time $t$ is in seconds. The particles are glass beads of
radius $r_p = 1\,\mathrm{mm}$ and density $\rho_p = 2500\,\mathrm{kg/m^3}$ in
still air of density $\rho_g = 1.2\,\mathrm{kg/m^3}$ and viscosity
$1.8 \times 10^{-5}\,\mathrm{Pa\,s}$. Drag with a constant drag coefficient
$C_D = 0.44$ slows each particle down,

```math
\frac{\mathrm{d}u}{\mathrm{d}t} = -c\,|u|\,u,
\qquad
c = \frac{3\,C_D\,\rho_g}{8\,\rho_p\,r_p} = 0.0792\,\mathrm{m^{-1}}.
```

A constant $C_D$ holds for particle Reynolds numbers above about $10^3$, which
these beads reach above about $7.5\,\mathrm{m/s}$; for the slower particles,
present from the start, it is an assumption of the model. The drag changes the
velocity of each particle and not their number, so it is the velocity
callback: it moves the abscissas and leaves the weights unchanged.

The collisions are modeled after Bhatnagar, Gross and Krook (BGK): they relax
the velocity distribution toward the normal distribution with the same number
density, mean and variance, at the constant rate $1/\theta$ with
$\theta = 0.1\,\mathrm{s}$,

```math
\frac{\mathrm{d}m_k}{\mathrm{d}t} = \frac{m_k^{\mathrm{eq}} - m_k}{\theta},
```

where $m_k^{\mathrm{eq}}$ are the moments of that normal distribution. They
leave $m_0$, $m_1$ and $m_2$, the number, momentum and energy of the particles,
unchanged. The constants are set in [physics.cpp](./physics.cpp), which
defines both callbacks.

[configuration.qmom](./configuration.qmom) selects the family with
`qmom_family = hyperbolic` and defines the initial number density as two
beams, normal distributions of $2 \times 10^4\,\mathrm{m^{-3}}$ particles at
$20 \pm 2\,\mathrm{m/s}$ and $10^4\,\mathrm{m^{-3}}$ at
$-10 \pm 4\,\mathrm{m/s}$. QMOMKit computes $m_0, \ldots, m_4$ from them, as
the family requires. The file has no `[time_integration]` section, so each step
uses the default Strang splitting, as described in [Evolving in
time](../../docs/user/5_evolving_in_time.md). [main.cpp](./main.cpp) advances
from $t = 0$ to $2\,\mathrm{s}$ with steps of $0.01\,\mathrm{s}$.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./hyqmom_1d_drag_bgk_collisions_example
```

## Output

At $t = 0$ and $t = 2\,\mathrm{s}$, the example prints the node count, the
number density $m_0$, the mean and standard deviation of the velocity, $m_5$
of the quadrature, and the five characteristic speeds. It writes the moments
and the quadrature every $0.05\,\mathrm{s}$ to `build/output/moments.txt` and
`build/output/quadrature.txt`, as set in the `[output]` section of
`configuration.qmom`.

## What to expect

The number density stays at $3 \times 10^4\,\mathrm{m^{-3}}$ and the
quadrature keeps three nodes. The drag lowers the mean velocity from $10$ to
$2.41\,\mathrm{m/s}$ and the standard deviation from $14.4$ to
$2.53\,\mathrm{m/s}$; the collisions change neither. The initial speeds, from
$-14.09$ to $21.86\,\mathrm{m/s}$, are those of the two beams. The collisions
merge the beams within a few $\theta$, and the five speeds then approach
$\mu - \sqrt{6}\,\sigma$, $\mu - \sigma$, $\mu$, $\mu + \sigma$ and
$\mu + \sqrt{6}\,\sigma$, those of a normal distribution of mean $\mu$ and
standard deviation $\sigma$ with node count 2, as $\mu$ and $\sigma$ keep
falling. At $t = 2\,\mathrm{s}$ they are $-3.97$, $-0.24$, $2.39$, $4.82$ and
$8.44\,\mathrm{m/s}$, against $-3.78$, $-0.12$, $2.41$, $4.94$ and
$8.60\,\mathrm{m/s}$ for a normal distribution. They are slightly lower
because the drag slows the fastest particles most, which leaves more particles
in a tail toward low velocities than a normal distribution has; the
collisions remove that tail only at the rate $1/\theta$. $m_5$ falls from
$6.86 \times 10^{10}$ to $6.88 \times 10^{7}\,\mathrm{m^2/s^5}$.
