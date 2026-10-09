# CQMOM 2D batch evolution with host data

This example shows how a host application advances many CQMOM states with the
batch-index-aware batch functions and gives the callbacks time-varying local
data. The host here owns eight locations, standing in for the mesh nodes held
by one processor. It packages its data for those locations, advances all of
them with one `advance_cqmom_phase_space_batch_with_index(...)` call per step,
and lets the callbacks read the data of the location they are evaluated for.

The physics is that of
[`cqmom_2d_breakage_advection_with-mc`](../cqmom_2d_breakage_advection_with-mc/README.md),
without the Monte Carlo companion. The particle velocity now relaxes toward a
local carrier flow supplied by the host, and the fragmentation rate is a local
host property.

## One index for host data and states

The host stores everything its callbacks need in arrays indexed by the local
location index `i` ([`host_data.hpp`](host_data.hpp)):

```cpp
struct LocalHostData
{
    std::vector<double> position;                  // used by the host's flow solver only
    std::vector<double> fragmentation_rate;        // read by the moment source
    std::vector<double> carrier_velocity_current;  // newest carrier-flow level
    std::vector<double> carrier_velocity_previous; // the level before it
    double sample_time = 0.0;
    double sample_interval = 0.0;
};
```

and keeps the CQMOM state of location `i` at position `i` of its state vector:

```cpp
LocalHostData host_data = make_local_host_data(location_count, initial_time);
std::vector<qmomkit::CqmomState> states(location_count, cqmom_prepared_case.initial_state());
```

A `*_batch_with_index(...)` call passes each state's position in `states` to
the callbacks as `batch_index`. Because the host uses the same index for its
own arrays, a callback reads entry `batch_index` directly, with no mapping in
between ([`physics.cpp`](physics.cpp)):

```cpp
return [&host_data](std::size_t location, double time, const std::vector<double>& coordinates,
               std::vector<double>& velocities, std::size_t dimensions)
{
    const double carrier_velocity = carrier_velocity_at(host_data, location, time);
    // ... fill velocities for every node of this location
};
```

The callbacks capture `host_data` by reference and read its current contents on
every invocation. The host therefore publishes new data between steps without
rebuilding the callbacks, and `host_data` must outlive them; passing a temporary
is rejected at compile time.

## The time loop

Each step the host publishes its new carrier-flow level at all locations, then
advances every location with one call, then checks every entry and takes the
new states back ([`main.cpp`](main.cpp)):

```cpp
while (time < final_time - 1.0e-12)
{
    const double step = std::min(time_step, final_time - time);
    if (time > initial_time)
    {
        publish_carrier_flow_level(host_data, time);
    }

    std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> evolved =
        qmomkit::advance_cqmom_phase_space_batch_with_index(states, time, step, velocity, source, evolution_options);

    for (std::size_t location = 0; location < location_count; ++location)
    {
        if (!evolved[location].ok())
        {
            // report the failed location and stop
        }
        states[location] = std::move(evolved[location].state);
    }
    time += step;
}
```

All entries of a batch share the time, the time step and the evolution options.
Each entry keeps its own status, diagnostics, adaptive substeps and retries, so
every result is checked. This example stops at the first failed location.

## Time-varying host data

The host cannot evaluate its flow field at the internal stage times QMOMKit
passes the velocity callback, and it does not know the field at `time + step`.
It keeps the level it computed last and the level before it, and
`carrier_velocity_at(...)` extends the straight line through them to the
requested stage time. This is second-order accurate in the time step, matching
the second-order Strang splitting. On the first step no previous level exists,
and the current level is held constant for that step. The same pattern, for a
single location, is explained in more detail in
[Supplying host data to a callback](../cqmom_5d_non-spherical-particles/README.md#supplying-host-data-to-a-callback)
in `cqmom_5d_non-spherical-particles`.

## Model

The phase-space ordering is $(d,u)$: particle size and particle velocity. At
location $x_i$ the velocity callback returns

$$
\frac{\mathrm{d}d}{\mathrm{d}t} = 1-d, \qquad
\frac{\mathrm{d}u}{\mathrm{d}t} = u_f(x_i,t)-u,
$$

with the travelling-wave carrier flow

```math
u_f(x,t) = \sin\!\big(2\pi(x - 0.5\,t)\big).
```

The moment source is binary fragmentation $d \to d/2$ at the local rate
$\lambda_i$, leaving $u$ unchanged:

$$
\frac{\mathrm{d}M_{j,k}}{\mathrm{d}t} = \lambda_i\left(2^{1-j}-1\right)M_{j,k}.
$$

The eight locations sit at $x_i = i/8$, and $\lambda_i$ rises linearly from
0.25 at $x_0$ to 0.6875 at $x_7$. Every location starts from the configured
independent product, $d \sim \mathcal{N}(2, 0.5^2)$ and
$u \sim \mathcal{N}(0, 0.5^2)$ with number density 1, on a
$2 \times 2$ hierarchy with conditioning order $(d, u)$.

## Expected behavior

The locations differ only through the host data, so the results show directly
that each callback read its own location's entries:

- **Number density.** Advection leaves it unchanged and fragmentation
  multiplies it by $e^{\lambda_i t}$, so its exact value is $n_0 e^{\lambda_i t}$.
  The summary prints the computed $M_{0,0}$ next to it, and the two agree to the
  printed digits.
- **Mean velocity.** Fragmentation does not change $u$, so the mean velocity
  $m = M_{0,1}/M_{0,0}$ obeys $\dot m = u_f(x_i,t) - m$. With
  $\phi_i = 2\pi x_i$ and $\omega = \pi$, its exact solution is

  $$m(t) = m_p(t) + \big(m(0) - m_p(0)\big)e^{-t}, \qquad m_p(t) = \frac{\sin(\phi_i - \omega t) + \omega\cos(\phi_i - \omega t)}{1 + \omega^2}.$$

  The particles follow the flow damped and delayed: the periodic part $m_p$
  has amplitude $1/\sqrt{1+\pi^2} \approx 0.30$ and lags the flow by about
  $72^\circ$, and the start-up transient from $m(0) = 0$ decays as $e^{-t}$.
  The summary prints the computed mean velocity next to $m(1)$. They differ by
  at most about $1.1\times10^{-4}$ with the shipped `time_step = 0.01`.
  Halving the step gives $2.8\times10^{-5}$, and halving it again gives
  $6.9\times10^{-6}$, which is second-order convergence.
- **Opposite locations.** Locations half a wavelength apart see opposite flows
  at all times, so location $i$ and location $i+4$ end with opposite mean
  velocities.

## Build and run

Requires an installed QMOMKit package. From the repository root:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build
```

Then configure and build this example from its directory:

```bash
cd examples/cqmom_2d_batch_host_data
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`. Run from the build
directory so that the output stays there:

```bash
cd build
./cqmom_2d_batch_host_data_example
```

The executable prints one summary row per location: its position and
fragmentation rate, then the computed $M_{0,0}$ and mean velocity, each next to
its exact value. With the shipped `[output]` settings it writes, for every location
`i`,

- `output/location_<i>_moment_history.txt`
- `output/location_<i>_quadrature_history.txt`

every 0.1 time units, through QMOMKit's native writers.

## File guide

| File | Purpose |
|---|---|
| `main.cpp` | Driver: loads the configuration, packages the host data, and runs the batch time loop |
| `configuration.qmom` | Hierarchy, initial condition, and output settings shared by all locations |
| `host_data.*` | Host-side data per location: positions, fragmentation rates, and the two carrier-flow levels |
| `physics.*` | Batch-index-aware velocity and moment-source callbacks reading the host data, and the exact solutions the summary compares with |
| `utilities.*` | Output setup and the per-location console summary |
