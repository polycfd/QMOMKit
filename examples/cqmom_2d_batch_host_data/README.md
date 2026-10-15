# 2D CQMOM batch evolution with host data

This example advances eight CQMOM states, standing in for the locations one
process of a host application owns, with one
`advance_cqmom_phase_space_batch_with_index(...)` call per step, as described
in [Many states at once](../../docs/user/6_coupling_qmomkit.md#many-states-at-once).
The callbacks receive each state's position in the vector as `batch_index` and
read the host data of that location, including a carrier flow that changes in
time.

## Model

The phase-space coordinates are the particle volume $V$ and the particle
velocity $u$, in that order; $V$, $u$, the position $x$ and the time $t$ are
dimensionless. At location $x_i$ the velocity callback returns

```math
\frac{\mathrm{d}V}{\mathrm{d}t} = 1-V, \qquad
\frac{\mathrm{d}u}{\mathrm{d}t} = u_f(x_i,t)-u,
```

with the traveling-wave carrier flow
$u_f(x,t) = \sin\!\big(2\pi(x - 0.5\,t)\big)$. The moment source is binary
breakage $V \to V/2$ at the local rate $\lambda_i$, leaving $u$ unchanged:

```math
\frac{\mathrm{d}M_{j,k}}{\mathrm{d}t} = \lambda_i\left(2^{1-j}-1\right)M_{j,k}.
```

The eight locations sit at $x_i = i/8$, and $\lambda_i$ rises linearly from
0.25 at $x_0$ to 0.6875 at $x_7$. The carrier flow and the rates are set in
[host_data.hpp](./host_data.hpp) and [host_data.cpp](./host_data.cpp); the
callbacks are defined in [physics.cpp](./physics.cpp).
[configuration.qmom](./configuration.qmom) starts every location from the same
independent product, $V \sim \mathcal{N}(2, 0.5^2)$ and
$u \sim \mathcal{N}(0, 0.5^2)$ with number density 1, on a $2 \times 2$
hierarchy with conditioning order $(V, u)$. Each step uses the default Strang
splitting. [main.cpp](./main.cpp) advances from $t = 0$ to $1$ with steps of
$0.01$.

### One index for host data and states

The host keeps the state of location `i` at position `i` of its state vector
and every array of its `LocalHostData` ([host_data.hpp](./host_data.hpp)) in
the same order, so a callback reads entry `batch_index` directly. The
callbacks capture the host data by reference, so the host publishes a new
carrier-flow level between steps without rebuilding them. Each step it
publishes the level, advances all locations with one call, and checks every
result; this example stops at the first failed location.

### Time-varying host data

QMOMKit evaluates the velocity at stage times inside the step, where the host
has no flow data. The host keeps its last two carrier-flow levels and
`carrier_velocity_at(...)` extrapolates linearly to the requested time, which
keeps the step second order, as described in [Host data that changes in
time](../../docs/user/6_coupling_qmomkit.md#host-data-that-changes-in-time). On
the first step only one level exists, and it is held constant.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_2d_batch_host_data_example
```

## Output

The example prints one row per location at $t = 1$: its position $x$, its
rate $\lambda$, the number density $M_{0,0}$ and the mean velocity
$M_{0,1}/M_{0,0}$, each next to its exact value. It writes the moments and the
quadrature of every location `i` every $0.1$ to
`build/output/location_<i>_moment_history.txt` and
`build/output/location_<i>_quadrature_history.txt`, as set in the `[output]`
section of `configuration.qmom`.

## What to expect

The locations differ only through the host data, so each row shows the
behavior of its own location:

- **Number density.** Advection leaves it unchanged and breakage multiplies
  it by $e^{\lambda_i t}$, from 1.28403 at $x_0$ to 1.98874 at $x_7$. The
  computed and exact columns agree to the printed digits.
- **Mean velocity.** Breakage does not change $u$, so the mean velocity
  $m = M_{0,1}/M_{0,0}$ obeys $\dot m = u_f(x_i,t) - m$. With
  $\phi_i = 2\pi x_i$ and $\omega = \pi$, its exact solution is

  ```math
  m(t) = m_p(t) + \big(m(0) - m_p(0)\big)e^{-t}, \qquad
  m_p(t) = \frac{\sin(\phi_i - \omega t)
  + \omega\cos(\phi_i - \omega t)}{1 + \omega^2}.
  ```

  The particles follow the flow damped and delayed: the periodic part $m_p$
  has amplitude $1/\sqrt{1+\pi^2} \approx 0.30$ and lags the flow by about
  $72^\circ$, and the start-up transient from $m(0) = 0$ decays as $e^{-t}$.
  The computed and exact mean velocities differ by at most about
  $1.1\times10^{-4}$ with the time step of $0.01$. Halving the step in
  `main.cpp` reduces the difference to $2.8\times10^{-5}$, and halving it
  again to $6.9\times10^{-6}$, the second-order convergence of the Strang
  splitting.
- **Opposite locations.** Locations half a wavelength apart see opposite
  flows at all times, so locations $i$ and $i+4$ end with opposite mean
  velocities.
