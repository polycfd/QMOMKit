# 5D CQMOM non-spherical particles

This example evolves a spatially uniform population of non-spherical particles
in a five-dimensional phase space with adaptive CQMOM and
`advance_cqmom_phase_space(...)`, following the work of
[Bouguezzoul (2026)](https://publications.polymtl.ca/76505/). It shows a
bounded support for one coordinate, a binary breakage source, and how a host
passes time-varying data, here a carrier flow, to a callback.
`main.cpp` runs QMOMKit, `main_mc.cpp` the Monte Carlo reference, which uses
the same initial sample and models but realizes breakage as stochastic
particle events.

## Model

The coordinates, in this order everywhere in the example, are

| Coordinate | Meaning | Unit |
|---|---|---|
| $d$ | particle diameter | m |
| $\chi$ | aspect ratio of the oblate spheroidal particles; 1 is a sphere | dimensionless |
| $u$, $v$, $w$ | particle velocity components | m/s |

For the number density $n(d,\chi,u,v,w)$, the raw moments are

```math
M_{a,b,c,e,g} = \int d^a\chi^b u^c v^e w^g\, n(d,\chi,u,v,w) \,\mathrm{d}d\,\mathrm{d}\chi\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w.
```

The carrier flow has the prescribed velocity
$\boldsymbol{U}_f(t) = (2\cos t,\;2\sin t,\;0)$ m/s, with $t$ in s. With the
particle velocity $\boldsymbol{U}_p = (u,v,w)$, the relative velocity
$\boldsymbol{U}_{rel} = \boldsymbol{U}_f-\boldsymbol{U}_p$ and
$\mathrm{Re} = \rho_f |\boldsymbol{U}_{rel}|d/\mu_f$, the drag factor of an
oblate spheroid is

```math
K(\chi,\mathrm{Re}) = K_0(\chi) + 0.15\,\chi^{95.91}\mathrm{Re}^{0.687} + 0.2927\,(1-\chi)^{0.4374}\mathrm{Re}^{0.7512},
```

with the creeping-flow shape factor

```math
K_0(\chi) = \frac{\frac{8}{3}\chi^{-1/3}}{\frac{2\chi}{1-\chi^2} + \frac{2(1-2\chi^2)}{(1-\chi^2)^{3/2}}\arctan\!\left(\frac{\sqrt{1-\chi^2}}{\chi}\right)},
```

replaced by the sphere value $K_0 = 1$ for $|1-\chi| < 0.001$. The particles
follow

```math
\begin{aligned}
\frac{\mathrm{d}\boldsymbol{U}_p}{\mathrm{d}t}
&= \frac{18\mu_f K}{\rho_p d^2}\boldsymbol{U}_{rel}-(0,0,g), \\
\frac{\mathrm{d}d}{\mathrm{d}t}
&= -\left(10^{-6}\ \mathrm{s/m}\right)
\frac{|\boldsymbol{U}_{rel}|^2}{1+|\boldsymbol{U}_{rel}|^2/(5\ \mathrm{m/s})^2}
\tanh\left(\frac{d}{10^{-5}\ \mathrm{m}}\right), \\
\frac{\mathrm{d}\chi}{\mathrm{d}t}
&= \frac{1-\chi}{2\ \mathrm{s}},
\end{aligned}
```

that is, drag and gravity, attrition driven by the relative speed, and
relaxation of the shape toward a sphere, with $\rho_f = 1.2$ kg/m³,
$\mu_f = 1.8\times10^{-5}$ Pa s, $\rho_p = 2000$ kg/m³ and $g = 9.81$ m/s².

Particles with $d \geq 10^{-6}$ m break at the rate $\beta(d) = 10^8 d^2$
(with $d$ in m and $\beta$ in 1/s) into two equal-volume daughters with
$d_d = 2^{-1/3}d_p$ and $\chi_d = 0.5\chi_p$. Each daughter receives a velocity
recoil with zero mean and standard deviation 0.5 m/s in every component. The
QMOMKit run adds at every quadrature node of the source stage the contribution
$\beta(d)\,[2\,\mathbb{E}(D^{a,b,c,e,g}) - P^{a,b,c,e,g}]$, where $P$ is the
parent monomial and $\mathbb{E}(D)$ the expected daughter monomial. The Monte
Carlo run samples exponentially distributed event times, and each event draws
one recoil vector and gives it to the two daughters with opposite signs.

The model parameters are set in [particle_model.hpp](particle_model.hpp) and
[breakage_model.hpp](breakage_model.hpp), the carrier flow in
[fluid_field.hpp](fluid_field.hpp). [configuration.qmom](configuration.qmom)
selects up to 3 nodes per coordinate (at most 243 nodes) with
`adaptive = true`, the identity conditioning order, Strang splitting and
adaptive RK23 node advection. Its `support_lower` and `support_upper` keys
bound $\chi$ to $[0,1]$, where the drag factor is defined, so the
reconstruction places no node outside it; without them a node with $\chi > 1$
makes $(1-\chi)^{0.4374}$ undefined, and QMOMKit reports the resulting
non-finite velocity as a failure of the step. The `[initial]` section
describes an independent product with $M_{0,0,0,0,0}=1$: a log-normal diameter,
a beta distribution with $\alpha=6$ and $\beta=14$ for $\chi$, an equal mixture
of two normal distributions with means $\pm1$ m/s and standard deviation
0.5 m/s for $u$, and normal $v$ and $w$ with zero mean and standard deviation
0.5 m/s.

Both executables draw the same 10000 particles from it
([sampling.hpp](sampling.hpp)). The QMOMKit run estimates from them only the
moments its hierarchy requires: `create_cqmom_moment_layout(...)` lists them,
`monte_carlo::sample_initialization(...)` estimates them, and
`config::with_initial_moments(...)` passes them to the case. Both runs advance
from $t=0$ to $t=1$ s with steps of 0.005 s, set in `main.cpp` and
`main_mc.cpp`; the Monte Carlo run applies half a step of advection,
stochastic breakage, and another half step.

### Supplying host data to a callback

A coupled host owns fields that QMOMKit knows nothing about, and a callback
must read their current local values without those fields becoming globals.
This example shows the pattern with the carrier flow, which the drivers own.
The driver hands the callback a reference to a struct that holds the two most
recent host levels:

```cpp
CarrierFlowSamples carrier_samples;
carrier_samples.current = carrier_flow_velocity(initial_time);
carrier_samples.sample_time = initial_time;
const qmomkit::PhaseSpaceVelocity velocity = particle_phase_space_velocity(particle_parameters, carrier_samples);

while (/* time loop */)
{
    if (time > initial_time)
    {
        carrier_samples.previous = carrier_samples.current;      // keep the level being replaced
        carrier_samples.current = carrier_flow_velocity(time);   // the host's new level
        carrier_samples.sample_interval = time - carrier_samples.sample_time;
        carrier_samples.sample_time = time;
    }
    // advance; the callback reads carrier_samples through its reference
}
```

The callback captures the reference once and reads through it on every
call, so the driver publishes new data instead of rebuilding the callback. A
host coupling, say, a concentration and a fluid velocity at cell `i` would
keep `concentration[i]` and `fluid_uvel[i]`, and their previous levels, in the
same struct. The struct must outlive the callback, which is why it is declared
before the callback and in the same scope; passing a temporary does not
compile.

The host cannot evaluate its fields at QMOMKit's stage times, and it does not
know them at `time + step`; it knows the level it just computed and the one
before. A single frozen value is only first-order accurate and would limit the
second-order Strang step. With two levels, the callback extrapolates linearly
to the stage time QMOMKit passes it,

```math
\boldsymbol{U}_f(s) \approx \boldsymbol{U}_f(t_n) + \frac{s-t_n}{\Delta t}\left[\boldsymbol{U}_f(t_n)-\boldsymbol{U}_f(t_{n-1})\right],
```

with an error of $O(\Delta t^2)$. The first step has no previous level and
holds the value constant, which costs first-order accuracy for that step only.

Measured with a copy of this example that evaluates the flow at the true stage
times as the reference, at the same `time_step`, over 5e-3, 2.5e-3 and 1.25e-3:

| Moment | Coordinate | Frozen value | Two-level reconstruction |
|---|---|---|---|
| $M_{0,0,1,0,0}$ | $u$ | 1.32e-02, order 1.02 | 4.30e-05, order 1.95 |
| $M_{0,0,0,1,0}$ | $v$ | 9.04e-03, order 1.02 | 6.24e-05, order 2.00 |
| $M_{0,0,2,0,0}$ | $u^2$ | 2.98e-02, order 1.02 | 9.15e-05, order 1.97 |

Errors are absolute, at $t = 1$ s and `time_step` = 5e-3; the orders are
averaged over the two refinements to 1.25e-3. Refining further, to 6.25e-4,
no longer reduces the two-level errors, which level off at about 2e-6, while
the frozen value stays first order. Mainly $u$ and $v$ are affected, the
components in which this flow varies. Staggering the particle step by half an
outer step, so that the host level falls at the step midpoint, would also give
second order with a single level.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_5d_non_spherical_particles_example ../configuration.qmom
./cqmom_5d_non_spherical_particles_mc_example ../configuration.qmom
```

Both executables take the configuration file as their only argument. The
output directory in it is relative to the working directory.

## Output

Both executables print the initial and final number density, first moments,
second moment of each coordinate, and $M_{1,1,0,0,0}$, which breakage couples;
the Monte Carlo executable also prints the initial and final particle counts.
The QMOMKit run writes `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`, the Monte Carlo run
`build/output/moment_history_mc.txt`, each every 0.1 s from $t=0$.

## What to expect

Breakage increases the number density, while the equal-volume daughters keep
the volume of each parent. Attrition, shape relaxation, drag, gravity and
recoil change the size, shape and velocity distributions; the particles
follow the rotating carrier flow and settle.
At $t = 1$ s the QMOMKit run gives $M_{0,0,0,0,0} = 3.093$ and the Monte Carlo
run 3.0995 from 30995 particles; the other printed moments of the two runs
differ by less than 1.7 percent, the Monte Carlo moments carrying the sampling
noise of a finite sample with stochastic breakage.
