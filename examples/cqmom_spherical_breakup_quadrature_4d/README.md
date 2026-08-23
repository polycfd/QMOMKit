# CQMOM spherical breakup: terminal quadrature source (4D)

## Purpose and workflows

This standalone example is a QMOMKit capability demonstration, not a detailed
physical breakup model. It evolves a spatially homogeneous population in phase
space $(d,u,v,w)$ with CQMOM and compares it with a direct Monte Carlo (MC)
reference solution. It demonstrates terminal-quadrature evaluation of a
nonlinear breakup/recoil source: the source deliberately uses
`CqmomSourceContext::for_each_terminal_node(...)`, rather than a moment
closure.

The breakup physics is intentionally simple: constant-rate binary breakup into
equal-volume daughters. However, the diameter-dependent recoil introduces higher-order 
moments that are not part of the transported finite moment set. The source is 
therefore evaluated by applying the nonlinear daughter mapping directly to the terminal 
quadrature, which provides a natural closure of the population transformation.

The directory contains two independent executables:

- `cqmom_spherical_breakup_quadrature_4d_example` samples an initial population,
  initializes a `CqmomCase` from its raw moments, reconstructs CQMOM, and
  advances the CQMOM state.
- `cqmom_spherical_breakup_quadrature_4d_mc_example` regenerates those same samples
  and advances the particles directly with stochastic binary breakup.

The MC reference is a finite-sample reference, not an expectation of bitwise
agreement with CQMOM. It samples individual breakup events and recoil vectors,
while the CQMOM source uses their conditional expected contribution over the
terminal quadrature.

## Phase space and raw moments

The phase-space coordinate order is

$$
\boldsymbol{\xi}=(d,u,v,w),
$$

where $d$ is diameter and $(u,v,w)$ is particle velocity. The raw,
unnormalized number-density moments are

$$
M_{p,q,r,s} =
\int d^p u^q v^r w^s f(d,u,v,w)
\,\mathrm{d}d\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w.
$$

Thus $M_{0,0,0,0}$ is total number density; QMOMKit does not normalize it. For
spherical particles, total particle volume is

$$
\mathcal V = \frac{\pi}{6}M_{3,0,0,0}.
$$

The same coordinate order is used by the configuration, sampled ensemble,
CQMOM hierarchy, source callback, and output files.

## Deterministic sampled initialization

Both executables first construct the same finite population with 10,000 samples
from the local independent-product initialization.
The CQMOM executable passes the raw tensor estimated from those samples to
`with_initial_moments(...)`; the MC executable uses the samples themselves.
Consequently, the paths begin from the identical finite-sample population, not
merely from equal analytical distributions. Each executable verifies this
identity internally before advancing.

The configured number density is $n_0=1$. The independent coordinates are

$$
\begin{aligned}
d &\sim \operatorname{LogNormal}(-9.210340371976184,0.12^2),\\
u &\sim \mathcal N(0.2,0.08^2),\\
v &\sim \mathcal N(-0.1,0.08^2),\\
w &\sim \mathcal N(0.05,0.08^2).
\end{aligned}
$$

Every initial sample represents the constant statistical weight $n_0/N$. The
independent MC event stream uses its own deterministic seed.

## Continuous phase-space evolution: Schiller–Naumann drag

The carrier gas is uniform and has velocity

$$
\boldsymbol U_g=(2.0,1.5,0.8)\ \mathrm{m\,s^{-1}},
$$

density $\rho_g=1.2\ \mathrm{kg\,m^{-3}}$, and dynamic viscosity
$\mu_g=1.8\times10^{-5}\ \mathrm{Pa\,s}$. The particles have material density
$\rho_p=2000\ \mathrm{kg\,m^{-3}}$. With $\boldsymbol s=\boldsymbol U_g-\boldsymbol u$,

$$
\mathrm{Re}=\frac{\rho_g d|\boldsymbol s|}{\mu_g}, \qquad
C_D=\frac{24}{\mathrm{Re}}\left(1+0.15\,\mathrm{Re}^{0.687}\right).
$$

The phase-space velocity callback returns zero diameter rate and the
componentwise drag acceleration

$$
\frac{\mathrm{d}\boldsymbol u}{\mathrm{d}t}
= \frac{3\rho_g}{4\rho_p d}C_D|\boldsymbol s|\boldsymbol s.
$$

Its low-Reynolds-number limit is evaluated directly as the Stokes response.
This application-specific callback is supplied to QMOMKit's generic
phase-space-advection operation.

## Binary breakup and recoil

Particles undergo constant-rate binary breakup,

$$
\beta=\beta_0, \qquad \beta_0=20\ \mathrm{s^{-1}}.
$$

Here $\beta_0$ directly means breakup events per particle per unit time; it
has no size, velocity, slip, or empirical-correlation scaling. A parent
produces two equal-volume daughters,

$$
d' = 2^{-1/3}d.
$$

For each daughter pair, the velocity is

$$
\boldsymbol u'_{\pm}=\boldsymbol u\pm\frac{d}{d_\mathrm{ref}}\boldsymbol\eta,
\qquad
\eta_i\sim\mathcal N(0,\sigma_r^2),
\qquad \sigma_r=10^{-3}\ \mathrm{m\,s^{-1}},
$$

with independent components and $d_\mathrm{ref}=10^{-4}\ \mathrm{m}$. Thus the
recoil scale is proportional to the parent diameter, while the two daughters
receive opposite recoil vectors. This deliberately nonlinear daughter-state
mapping is the breakup physics used to demonstrate the closure mechanism; it
is not intended as a detailed recoil model. In the CQMOM source, the daughter
monomial is the corresponding Gaussian expectation. For every requested
integer raw moment, the source evaluates the parent loss and two-daughter birth
term at each node yielded by `for_each_terminal_node(...)` and sums with its
terminal weight.

Terminal quadrature is required for this source. A velocity factor of even
Gaussian order $2k$ introduces $(d/d_\mathrm{ref})^{2k}$ into the integrand;
for a moment $M_{p,q,r,s}$, its birth rate therefore needs parent moments with
diameter order $p+2k$ (and analogous contributions from $r$ and $s$). Those
higher diameter moments are outside the finite CQMOM transport layout—for
example, the rate of tracked $M_{2,4,0,0}$ needs $M_{6,0,0,0}$. The terminal
quadrature supplies this reconstruction-dependent closure without a
moment-based breakup formula.

The direct MC counterpart samples the same constant-rate event process and
symmetric pair at each event, allowing later daughter events within the source
interval. Starting from $N_0$, breakup alone gives

$$
N(t)=N_0\exp(\beta_0t).
$$

Over the $0.1\ \mathrm{s}$ run, the expected number-density increase is
$\exp(2)\approx7.39$. Breakup alone conserves total volume and conserves
volume-weighted momentum because the two recoils cancel. The velocity-square
raw moment captures the additional recoil spread. Carrier-gas drag acts
continuously, so it can change the population momentum over the full coupled
run.

## CQMOM configuration

`configuration.qmom` is the source of standard QMOMKit numerical behavior:

- phase-space names are `d, u, v, w`;
- maximum node counts are `(3,3,3,3)` with conditioning order `(0,1,2,3)`;
- `adaptive_branch_counts` is `true`;
- the phase-space/moment-source composition uses Strang splitting;
- node advection uses adaptive RK23 with the listed tolerances and step bounds;
- output cadence and enabled moment/quadrature files come from `[output]`.

The application runs from $t=0$ to $0.1\ \mathrm{s}$ with a macro step of
$0.001\ \mathrm{s}$. These time controls and the physical constants are local
example choices; hierarchy and standard numerical options remain in the
configuration file.

## Build and run

Install QMOMKit and the optional QMOMKitMonteCarlo companion into a common
prefix, then configure this standalone consumer:

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

cmake -S examples/cqmom_spherical_breakup_quadrature_4d \
    -B examples/cqmom_spherical_breakup_quadrature_4d/build \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build examples/cqmom_spherical_breakup_quadrature_4d/build

cd examples/cqmom_spherical_breakup_quadrature_4d/build
./cqmom_spherical_breakup_quadrature_4d_example
./cqmom_spherical_breakup_quadrature_4d_mc_example
```

## Output and interpretation

Both executables print the initial and final tracked raw moments in the CQMOM
transported-moment order. The console follows the normal QMOMKit example style:
it prints the first twelve tracked moments as `M[p,q,r,s]` and directs the
reader to the complete history when more are tracked. The CQMOM and MC labels
and ordering are identical, so their values can be compared directly.

When enabled, the CQMOM executable writes `output/moment_history.txt` and
`output/quadrature_history.txt`. The MC executable writes the corresponding
selected raw moments, in the same order, to `output/moment_history_mc.txt`.
Compare matching rows at equal output times while allowing finite-MC event and
recoil sampling variability. For example, $M_{3,0,0,0}$ is proportional to
total volume and remains constant; the number-density moment rises under
breakup.
