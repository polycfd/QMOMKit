# CQMOM spherical breakup: moment-based source term (5D)

## Purpose and workflows

This standalone example evolves a spatially homogeneous spherical-particle
population in CQMOM and compares it with a direct Monte Carlo (MC) reference.
It combines Schiller–Naumann drag, Ranz–Marshall heat transfer, and equal-volume
binary breakup. Its CQMOM breakup source is written directly in moment form;
it does not traverse the terminal quadrature.

The directory contains two independent executables:

- `cqmom_spherical_breakup_moments_5d_example` samples the initial population,
  converts its raw moments into the CQMOM initial state, and advances CQMOM.
- `cqmom_spherical_breakup_moments_5d_mc_example` regenerates the same samples and
  advances them with direct stochastic breakup events.

The MC reference is a finite-sample reference. It represents event sampling
directly, while the CQMOM source evolves the analytical conditional moment
rates. Close agreement in compatible tracked raw moments is expected, but not exact
agreement at finite MC sample count.

## Phase space and raw moments

The configured coordinate order is

$$
\boldsymbol{\xi}=(V,T,u,v,w),
$$

where $V$ is particle volume, $T$ is particle temperature, and $(u,v,w)$ is
particle velocity. The raw, unnormalized number-density moments are

$$
M_{p,q,r,s,t} =
\int V^p T^q u^r v^s w^t f(V,T,u,v,w)
\,\mathrm{d}V\,\mathrm{d}T\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w.
$$

In particular, $M_{0,0,0,0,0}$ is number density and $M_{1,0,0,0,0}$ is total
particle volume. The physical diameter used by drag and heat transfer is

$$
d(V)=\left(\frac{6V}{\pi}\right)^{1/3}.
$$

The fractional power of $V$ is evaluated directly at quadrature nodes by the
continuous drag and heat-transfer callback. The constant-rate breakup source
uses only tracked integer raw moments, avoiding artificial generalized-moment
queries.

## Deterministic sampled initialization

Both executables first regenerate 16,000 samples from the local
independent-product initialization using seed `2026082311`. The CQMOM path
uses the complete raw tensor estimated from this population with
`with_initial_moments(...)`; the MC path evolves the population itself. They
therefore start from the same finite-sample population rather than separately
drawn realizations of equal distributions. Each executable verifies that the
moments recomputed from its ensemble exactly match the CQMOM input tensor
before advancing.

The configured total number density is $n_0=1$, with

$$
\begin{aligned}
V &\sim \operatorname{LogNormal}(-28.278651581805676,0.36^2),\\
T &\sim \mathcal N(350,4^2),\\
u &\sim \mathcal N(0.2,0.08^2),\\
v &\sim \mathcal N(-0.1,0.08^2),\\
w &\sim \mathcal N(0.05,0.08^2).
\end{aligned}
$$

Each initial sample carries weight $n_0/N$. The direct MC breakup stream uses
the independent deterministic seed `2026082312`.

## Drag and heat transfer

The uniform carrier gas has velocity

$$
\boldsymbol U_g=(2.0,1.5,0.8)\ \mathrm{m\,s^{-1}},
$$

density $\rho_g=1.2\ \mathrm{kg\,m^{-3}}$, viscosity
$\mu_g=1.8\times10^{-5}\ \mathrm{Pa\,s}$, temperature $T_g=300\ \mathrm{K}$,
thermal conductivity $k_g=0.026\ \mathrm{W\,m^{-1}\,K^{-1}}$, and Prandtl
number $\mathrm{Pr}=0.71$. The particle material density is
$\rho_p=2000\ \mathrm{kg\,m^{-3}}$ and heat capacity is
$c_p=1000\ \mathrm{J\,kg^{-1}\,K^{-1}}$.

With $\boldsymbol s=\boldsymbol U_g-\boldsymbol u$,

$$
\mathrm{Re}=\frac{\rho_g d|\boldsymbol s|}{\mu_g},\qquad
C_D=\frac{24}{\mathrm{Re}}\left(1+0.15\,\mathrm{Re}^{0.687}\right),
$$

and the velocity rate is

$$
\frac{\mathrm{d}\boldsymbol u}{\mathrm{d}t}
=\frac{3\rho_g}{4\rho_p d}C_D|\boldsymbol s|\boldsymbol s.
$$

The zero-slip limit uses the equivalent Stokes response. Ranz–Marshall heat
transfer supplies

$$
\mathrm{Nu}=2+0.6\,\mathrm{Re}^{1/2}\mathrm{Pr}^{1/3},
$$

and the lumped particle-temperature equation is

$$
\frac{\mathrm{d}T}{\mathrm{d}t}
=\frac{6k_g\mathrm{Nu}}{\rho_p c_p d^2}(T_g-T).
$$

Volume has no continuous advection rate. These rates are provided through the
application-local phase-space velocity callback; QMOMKit provides the generic
phase-space integration and splitting machinery.

## Constant-rate equal-volume breakup

Particles undergo constant-rate binary breakup,

$$
\beta=\beta_0, \qquad \beta_0=20\ \mathrm{s^{-1}}.
$$

Thus $\beta_0$ directly means breakup events per particle per unit time, with
no size, velocity, slip, or empirical-correlation scaling. Each event replaces
a parent by two daughters with $V'=V/2$; temperature and velocity are inherited
exactly, and there is no recoil.

The source of $M_{p,q,r,s,t}$ is

$$
\dot M_{p,q,r,s,t}^{\mathrm{br}}
=\beta_0\left(2^{1-p}-1\right)M_{p,q,r,s,t}.
$$

The factor $2^{1-p}-1$ is the exact daughter-minus-parent factor. This is a
moment-based source: it reads each tracked raw moment directly and does not
call `for_each_terminal_node(...)`. For $p=1$, the factor is zero, so total
volume is exactly conserved by breakup. The inherited daughter velocity and
temperature likewise conserve volume-weighted momentum and volume-weighted
temperature under breakup alone. Starting from $N_0$,

$$
N(t)=N_0\exp(\beta_0t),
$$

so the $0.1\ \mathrm{s}$ run has expected number-density increase
$\exp(2)\approx7.39$. The MC executable samples this same constant-rate binary
event process, allowing further daughter events in the same source interval.

## CQMOM configuration

`configuration.qmom` controls the standard QMOMKit setup:

- phase-space names are `V, T, u, v, w`;
- maximum node counts are `(3,3,3,3,3)` with conditioning order `(0,1,2,3,4)`;
- `adaptive_branch_counts` is `true`;
- the automatic Vandermonde solver is configured;
- phase-space evolution uses Strang splitting;
- adaptive RK23 node advection uses the stated relative/absolute tolerances
  and step bounds;
- `[output]` controls history files and their $0.01\ \mathrm{s}$ cadence.

The example-local run extends from $0$ to $0.1\ \mathrm{s}$ with macro step
$0.001\ \mathrm{s}$. Physical properties and these application controls are
compiled into the example; CQMOM hierarchy and standard numerical behavior are
configuration controlled.

## Build and run

Install QMOMKit and QMOMKitMonteCarlo into a common prefix, then configure this
standalone package:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build

cmake -S examples/cqmom_spherical_breakup_moments_5d \
    -B examples/cqmom_spherical_breakup_moments_5d/build \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build examples/cqmom_spherical_breakup_moments_5d/build

cd examples/cqmom_spherical_breakup_moments_5d/build
./cqmom_spherical_breakup_moments_5d_example
./cqmom_spherical_breakup_moments_5d_mc_example
```

## Output and interpretation

Both executables print the initial and final tracked raw moments in the CQMOM
transported-moment order. The console follows the normal QMOMKit example style:
it prints the first twelve tracked moments as `M[p,q,r,s,t]` and directs the
reader to the complete history when more are tracked. The CQMOM and MC labels
and ordering are identical, so their values can be compared directly.

When enabled, the CQMOM executable writes `output/moment_history.txt` and
`output/quadrature_history.txt`. The MC executable writes the corresponding
selected raw moments, in the same order, to `output/moment_history_mc.txt`.
Compare matching rows at equal output times. Number density increases under
breakup, total volume is $M_{1,0,0,0,0}$ and remains constant, and temperature
relaxes toward the carrier temperature. The CQMOM and MC histories should be
close, with finite-sample deviations caused by the direct MC event realization.
