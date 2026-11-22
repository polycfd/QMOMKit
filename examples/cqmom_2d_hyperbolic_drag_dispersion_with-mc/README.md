# 2D CQMOM with a hyperbolic velocity, drag and turbulent dispersion, with a Monte Carlo reference

This example shows a `hyperbolic` coordinate in a CQMOM hierarchy
([QMOM families](../../docs/user/7_qmom_families.md)): particles in a flow,
each with a velocity and a diameter. Stokes drag pulls every particle's velocity
toward the velocity of the flow, small particles faster than large ones, and
turbulent dispersion spreads the velocities. The same case runs in both
conditioning orders, beside a Monte Carlo run that starts from the same finite
sample. `main.cpp` loads a `qmomkit::CqmomPreparedCase` from
[velocity_first.qmom](./velocity_first.qmom) or [size_first.qmom](./size_first.qmom),
replaces its initial moments with those of the sample through
`qmomkit::config::with_initial_moments(...)` and advances them with
`qmomkit::advance_cqmom_phase_space(...)`. `main_mc.cpp`, the Monte Carlo
reference, advances the sampled particles themselves with the Monte Carlo
companion.

## Model

The phase-space coordinates are, in this order, the particle velocity $u$ in
m/s and the particle diameter $d$ in m, and $n(u,d)$ is the number density.
The raw moments are
$M_{k,l} = \iint u^k d^l\, n(u,d)\,\mathrm{d}u\,\mathrm{d}d$. The population is
spatially uniform. Stokes drag relaxes the velocity toward the flow velocity
$U$ with the response time $\tau_p$, and the diameter does not change,

```math
\frac{\mathrm{d}u}{\mathrm{d}t} = \frac{U - u}{\tau_p(d)},
\qquad
\tau_p(d) = \frac{\rho_p d^2}{18\,\mu_f},
\qquad
\frac{\mathrm{d}d}{\mathrm{d}t} = 0 .
```

Turbulent dispersion is modeled as velocity diffusion with the coefficient
$D(d) = \sigma_t^2/\tau_p(d)$, so that the velocity of every size relaxes to the
variance $\sigma_t^2$ about $U$. It changes the moments as

```math
\frac{\mathrm{d}M_{k,l}}{\mathrm{d}t} = k(k-1) \iint D(d)\, u^{k-2} d^l\, n(u,d)\,\mathrm{d}u\,\mathrm{d}d ,
```

which the moment source of the example evaluates over the nodes of the
quadrature that QMOMKit reconstructs at each stage, since $D$ is not a
polynomial in $d$. The particles have the density $\rho_p = 1000$ kg/m³ and the
fluid the viscosity $\mu_f = 1.8 \times 10^{-5}$ Pa s. The flow velocity is
$U = 0.5$ m/s, and $\sigma_t = 0.1$ m/s ([physics.hpp](./physics.hpp)). The
model assumes Stokes drag, which holds for particle Reynolds numbers
$\rho_f |U-u|\, d/\mu_f$ up to about 1, as for most of the population with a
fluid density $\rho_f = 1.2$ kg/m³, and dispersion as white noise. It has no
collisions, breakup, evaporation or gravity.

In the Monte Carlo run, each particle's velocity receives the kick
$\sqrt{2D(d)\,\Delta t}\,\xi$ with $\xi$ standard normal in every step, which is
exact for diffusion with a constant coefficient over the step
([particles_mc.cpp](./particles_mc.cpp)). All three runs use Strang splitting:
half a step of drag, the dispersion over the whole step, and a second half step
of drag, for $0 \le t \le 10$ ms with $\Delta t = 0.05$ ms.

The initial number density is the product of the normal distribution
$u \sim \mathcal{N}(0, (0.02\ \text{m/s})^2)$ and the log-normal distribution of
$d$ with a median of 20 µm and `sigma_log = 0.35`, with a total of one particle
per cubic meter, which only scales the moments, since the particles do not
interact. 95 % of the diameters lie between about 10 and 40 µm, with response
times between about 0.3 and 5 ms. All three executables draw the same
$N = 10000$ samples with a fixed seed ([sampling.hpp](./sampling.hpp)).

## The two conditioning orders

Both configuration files set `qmom_family = hyperbolic, classical` and
`node_counts = 2, 2`, with the default adaptive reconstruction. The velocity
uses HyQMOM, which reads the velocity moments up to order 4 and holds up to
three nodes, and the diameter classical QMOM with up to two nodes. The files
differ only in `conditioning_order`:

- `velocity_first.qmom`, `conditioning_order = 0, 1`: up to three velocity
  nodes, and up to two diameter nodes conditioned on each. The diameter moments
  conditioned on three velocity nodes need the mixed moments with velocity
  orders 0, 1 and 2, such as $M_{2,1}$.
- `size_first.qmom`, `conditioning_order = 1, 0`: up to two diameter nodes, and
  up to three velocity nodes conditioned on each.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_2d_hyperbolic_drag_dispersion_with_mc_velocity_first_example
./cqmom_2d_hyperbolic_drag_dispersion_with_mc_size_first_example
./cqmom_2d_hyperbolic_drag_dispersion_with_mc_mc_example
```

## Output

Each executable prints the initial and final raw moments $M_{0,0}$, $M_{1,0}$,
$M_{0,1}$, $M_{2,0}$, $M_{1,1}$, $M_{0,2}$ and $M_{2,1}$. Every millisecond, the
QMOMKit runs write `build/output/moment_history_velocity_first.txt` and
`build/output/quadrature_history_velocity_first.txt`, or the same files ending
in `size_first`, and the Monte Carlo run writes
`build/output/moment_history_mc.txt`.

## What to expect

All three runs start from the same sample, so their initial moments agree. The
mean velocity $M_{1,0}/M_{0,0}$ rises from zero toward $U$ and reaches about
0.49 m/s at 10 ms, and $M_{1,1}$ ends below $M_{1,0}M_{0,1}/M_{0,0}$, since the
large particles lag behind the small ones. With the diameter first, both
diameter branches keep their three velocity nodes, and the final moments lie
within 1 % of the Monte Carlo run. With the velocity first, the velocity keeps
its three nodes, but during the first 3 ms some diameter branches hold a single
node, where the diameter moments conditioned on a velocity node are not
realizable. The run then keeps the moments of that reduced quadrature, so
$M_{0,2}$, which does not change in this model, ends about 10 % high. The Monte
Carlo run is a stochastic reference with sampling noise, not an exact solution.
