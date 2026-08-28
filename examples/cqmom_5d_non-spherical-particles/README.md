# 5D CQMOM non-spherical-particle evolution

This standalone example evolves a spatially uniform population of
non-spherical particles in a five-dimensional phase space with recursive
CQMOM and an independent stochastic particle-based Monte Carlo (MC)
companion/reference. It combines configuration-controlled bounded support,
adaptive CQMOM branch counts, Strang-split phase-space evolution, and a
moment-source representation of binary breakup. The MC executable uses the
same sampled initial population and physical rules, but realizes breakup as
stochastic particle events.

The example concerns phase-space evolution only. It does not model a mesh,
physical-space transport, particle-particle collisions, or feedback to a
carrier-flow solver.

## Phase space and moments

The coordinate order is fixed everywhere in the example:

```text
(d, chi, u, v, w)
```

The same order is used in `configuration.qmom`, sampled particles, CQMOM
moment indices, quadrature abscissas, the particle-velocity callback, and
output files.

| Coordinate | Meaning | Physical/support information |
|---|---|---|
| `d` | Particle diameter | m |
| `chi` | Aspect-ratio shape parameter used by the oblate-spheroid drag correction | dimensionless; the particle correlation requires `0 <= chi <= 1` |
| `u` | Particle velocity in the first Cartesian direction | m/s |
| `v` | Particle velocity in the second Cartesian direction | m/s |
| `w` | Particle velocity in the third Cartesian direction | m/s |

Diameter and translational velocity alone do not describe the drag response of
this non-spherical model: the drag correction also depends on `chi`. The
implemented correlation treats the near-`chi = 1` limit as spherical and uses
an oblate-spheroid form away from that limit.

For a number-density distribution `f(d, chi, u, v, w)`, the raw moments are

$$
M_{a,b,c,e,g} =
\int d^a\chi^b u^c v^e w^g f(d,\chi,u,v,w)
\,\mathrm{d}d\,\mathrm{d}\chi\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w.
$$

In particular, `M[0,0,0,0,0]` is total number density; QMOMKit does not
implicitly normalize it to one.

## Particle phase-space model

The particle callback returns rates for `(d, chi, u, v, w)`. The prescribed
carrier velocity is time dependent:

$$
\boldsymbol{U}_f(t) =
(2\cos t,\;2\sin t,\;0)\ \mathrm{m/s}.
$$

With particle velocity `\boldsymbol{U}_p = (u,v,w)`, relative velocity
`\boldsymbol{U}_{rel} = \boldsymbol{U}_f-\boldsymbol{U}_p`, and

$$
\mathrm{Re} = \frac{\rho_f |\boldsymbol{U}_{rel}|d}{\mu_f},
$$

the shape-dependent drag factor is evaluated as

$$
K(\chi,\mathrm{Re}) = K_0(\chi)
+ 0.15\,\chi^{95.91}\mathrm{Re}^{0.687}
+ 0.2927\,(1-\chi)^{0.4374}\mathrm{Re}^{0.7512}.
$$

`K_0` is the implemented oblate-spheroid contribution. The velocity rates use
the drag coefficient

$$
\lambda = \frac{18\mu_fK}{\rho_p d^2},
\qquad
\frac{\mathrm{d}\boldsymbol{U}_p}{\mathrm{d}t}
= \lambda\boldsymbol{U}_{rel}-(0,0,g),
$$

with `rho_f = 1.2 kg/m^3`, `mu_f = 1.8e-5 Pa s`, `rho_p = 2000 kg/m^3`, and
`g = 9.81 m/s^2`.

Diameter changes by an attrition rate with a saturated relative-speed response
and a hyperbolic-tangent taper near the configured zero minimum diameter. With
`s = |U_rel|`, its implemented form is

$$
\frac{\mathrm{d}d}{\mathrm{d}t} =
-10^{-6}\,
\frac{s^2}{1+s^2/(5\ \mathrm{m/s})^2}
\tanh\left(\frac{d}{10^{-5}\ \mathrm{m}}\right).
$$

The shape parameter relaxes toward the spherical limit on a two-second time scale:

$$
\frac{\mathrm{d}\chi}{\mathrm{d}t} = \frac{1-\chi}{2\ \mathrm{s}}.
$$

The physical domain of the shape correlation is `0 <= chi <= 1`. The model
does not clamp `chi` in the velocity callback.

## Binary breakup

Particles with `d >= 1e-6 m` break at the diameter-dependent rate

$$
\beta(d) = 10^8 d^2.
$$

Each event removes one parent and produces two equal-volume daughters:

$$
d_d = 2^{-1/3}d_p, \qquad \chi_d = 0.5\chi_p.
$$

The diameter transformation preserves the parent volume because
`2 d_d^3 = d_p^3`. Each daughter receives a centered Gaussian velocity recoil
with standard deviation `0.5 m/s` in every velocity coordinate. In the MC
event, one sampled recoil vector is added to one daughter and subtracted from
the other, so the pair has equal and opposite recoil.

The CQMOM executable represents the same process as a moment source. At each
source-stage terminal quadrature node it adds the breakage contribution

$$
\beta(d)\,[2\,\mathbb{E}(D^{a,b,c,e,g}) - P^{a,b,c,e,g}],
$$

where `P` is the parent monomial and the daughter expectation uses the
equal-volume diameter, halved shape parameter, and centered Gaussian recoil
moments. Thus CQMOM uses a deterministic moment-source representation,
whereas the MC companion samples exponentially distributed event times and
creates daughter particles. Both represent the same underlying breakup model,
but a finite stochastic MC realization is not expected to match the CQMOM
trajectory exactly.

## Shared deterministic initialization

The supplied case is an independent-product distribution with total initial
number density `n0 = 1.0`. The five marginals are sampled independently:

| Coordinate | Configured distribution |
|---|---|
| `d` | Log-normal with `mu_log = -8.9596725473`, `sigma_log = 0.1529473098` |
| `chi` | Beta with `alpha = 6.0`, `beta = 14.0` |
| `u` | Equal-number-density mixture of `N(-1.0, 0.5^2)` and `N(1.0, 0.5^2)` m/s |
| `v` | `N(0.0, 0.5^2)` m/s |
| `w` | `N(0.0, 0.5^2)` m/s |

Both drivers use the shared sampling helper with

```text
sample_count = 10000
initialization_seed = 12345
```

The MC driver retains those sampled particles and assigns each the constant
initial statistical weight `n0 / sample_count`. The CQMOM driver computes its
complete initial raw-moment tensor from that *same* deterministic finite
ensemble, replaces the loaded case's initial moments with that tensor, and
then reconstructs its initial quadrature. Consequently the initial sampled
moments are identical between the two representations; analytical moments of
the configured distributions are not substituted for the production initial
condition.

## CQMOM configuration and evolution

[`configuration.qmom`](configuration.qmom) supplies native QMOMKit numerical
configuration; the C++ drivers load it and do not override these settings.

| Setting | Supplied value and effect |
|---|---|
| Phase-space names/order | `d, chi, u, v, w` |
| Node counts | `3, 3, 3, 3, 3`; up to three conditional nodes per level and at most 243 terminal nodes |
| Conditioning order | `0, 1, 2, 3, 4`, reconstructing in the displayed coordinate order |
| Branch behavior | `adaptive_branch_counts = true`; the configured counts are capacities and the accepted topology may use fewer branches |
| Support | `chi` lower/upper bounds `0.0` and `1.0`; all other coordinate bounds are unset |
| Solver/diagnostics | no solver override in the file (the case uses QMOMKit's automatic policy); `diagnostic_level = none` |
| Sampled-moment orders | no `maximum_orders` entry; the shared sampler uses hierarchy extents, giving orders `0` through `5` in each coordinate |
| Splitting | `strang` |
| Node advection | adaptive `rk23`, relative tolerance `1e-6`, absolute tolerance `1e-9`, initial/max/min substeps `1e-3`, `1e-2`, `1e-9` |
| File output | directory `output`, moments and quadrature enabled, interval `0.1` |

The outer driver controls are local constants, rather than `.qmom` entries:

```text
initial_time = 0.0
final_time = 1.0
time_step = 0.005
```

The CQMOM driver advances a single `CqmomState` with
`advance_cqmom_phase_space(...)`. With Strang splitting this applies a
half-step of phase-space advection, a full breakage-source step, and a second
advection half-step. The MC driver performs the corresponding half advection,
stochastic breakup, half advection sequence for its particle ensemble, using
the same RK23 options for the continuous stages.

## Why the aspect-ratio support matters

The support declarations in the supplied configuration constrain **CQMOM
reconstruction** to `chi in [0,1]`, which is the physical domain required by
the particle correlation. They are not particle clipping, and the particle
callback does not impose a separate `chi` guard. The MC particles are also not
clipped merely because the CQMOM configuration contains support bounds.

This distinction is important when deliberately removing the support. Make a
temporary copy of `configuration.qmom`, remove or comment its two `chi`
support entries, and pass that copy explicitly to the CQMOM executable. Do not
alter the supplied configuration file.

In the currently verified support-free case, an unconstrained reconstruction
produces a node with `chi = 1.0921392410871662` at `t = 0.1175`. The root
cause and the visible terminal error are related as follows:

```text
unbounded reconstruction produces chi > 1
    -> particle state leaves the correlation's physical domain
    -> 1 - chi < 0
    -> (1 - chi)^0.4374 is not real
    -> the particle-velocity callback returns NaN
    -> QMOMKit reports a non-finite phase-space velocity
```

The terminal message therefore names a non-finite `u`, `v`, or `w` velocity
component rather than explicitly reporting `chi > 1`. The non-finite velocity
is the downstream symptom; the original cause is that reconstruction produced
an aspect ratio outside `[0,1]`. This is not a QMOMKit solver failure and
QMOMKit did not create the NaN before the physical callback. The purpose of
the configured support is to keep reconstructed quadrature nodes inside the
model domain, not to hide the issue through clipping.

## Build and run

This is a standalone installed-package consumer. Install both QMOMKit and the
optional QMOMKitMonteCarlo companion, then configure this example from the
example directory (`examples/cqmom_5d_non-spherical-particles`) using

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

This produces the following executables in the build directory:

- `cqmom_5d_non_spherical_particles_example`
- `cqmom_5d_non_spherical_particles_mc_example`

Execute the code from the example root directory with

```sh
./build/cqmom_5d_non_spherical_particles_example configuration.qmom
./build/cqmom_5d_non_spherical_particles_mc_example configuration.qmom
```

Run from a directory where writing `output/` is appropriate: the configured
output directory is relative to the process's working directory.

## Output and expected behavior

For the normal bounded run, the CQMOM executable prints its selected initial
and final moments, final terminal quadrature, configuration/evolution settings,
and output paths. It writes:

```text
output/initial_moments.txt
output/moment_history.txt
output/initial_quadrature.txt
output/quadrature_history.txt
```

The moment files contain time and raw-moment columns. The quadrature files
contain terminal weights, conditioning paths, and abscissas ordered as
`(d, chi, u, v, w)`.

The MC executable prints corresponding initial and final sampled moments and
the final ensemble size. It writes:

```text
output/moment_history_mc.txt
```

Breakup increases number density by replacing parents with daughter pairs,
while equal-volume diameter scaling preserves volume at each breakup event.
Attrition, shape relaxation, drag, gravity, and recoil then change the size,
shape, and velocity distributions. On the currently validated default run,
the final number densities are approximately `3.0924079105` for CQMOM and
`3.0995` for MC; the MC ensemble contains `30995` particles. Small later-time
differences are expected from the finite CQMOM representation, numerical
integration, and stochastic MC breakup realization. The supplied bounded case
should complete successfully; the temporary support-removal experiment is
expected to stop through the non-finite-velocity path described above.

## File guide

| File | Purpose |
|---|---|
| `main.cpp` | Native CQMOM driver: sampled raw moments, initial reconstruction, and phase-space evolution |
| `main_mc.cpp` | Direct stochastic particle-based MC companion |
| `configuration.qmom` | Numerical hierarchy, support, integration, initialization, and output settings |
| `particle_model.*` | Non-spherical particle phase-space velocity callback and example parameters |
| `breakage_model.*` | CQMOM binary-breakage source and matching MC event model |
| `particles_mc.*` | Raw-moment estimation from the evolving MC ensemble |
| `sampling.*` | Shared fixed-seed independent-product initialization |
| `utilities.*` | Output setup and concise console summaries |

The two workflows intentionally remain direct:

```text
configuration.qmom -> main.cpp    -> native CQMOM phase-space evolution
same sampled population -> main_mc.cpp -> stochastic particle evolution
```
