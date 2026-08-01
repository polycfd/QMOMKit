# Condensing Droplet Population

## Purpose

This example models a polydisperse population of isolated water droplets. Each
droplet falls under gravity, experiences aerodynamic drag, grows through
condensation, and exchanges sensible and latent heat with the surrounding air.

The example contains continuous phase-space evolution only. It does not model
breakup, coalescence, nucleation, a discontinuous population source, a
zero-radius boundary contribution, or physical-space advection. The signed
mass-transfer equation is retained in the droplet model, but the configured
initial and ambient conditions produce condensation in the default case.

The physical model is based on:

```text
Fabian Denner,
"Water droplet dynamics and evaporation in airtanker firefighting,"
International Journal of Multiphase Flow 199 (2026) 105696.
DOI: 10.1016/j.ijmultiphaseflow.2026.105696
```

## Phase Space

The configured phase-space order is:

```text
(r, w, T)
```

where:

- `r` is droplet radius in metres;
- `w` is vertical droplet velocity in metres per second;
- `T` is droplet temperature in kelvin.

Upward velocity is positive, so falling droplets have negative `w`. The same
positional order is used by the initial distribution, physical callback, CQMOM
moments, Monte Carlo samples, and output.

## Initial Population

The initial condition is an independent product distribution:

```text
r ~ Gamma(shape = 16.0, scale = 1.875e-5)
w ~ Normal(mean = -2.0, standard_deviation = 1.0)
T ~ Normal(mean = 293.15, standard_deviation = 3.0)
```

The number density is `1.0`. The gamma radius distribution has mean:

```text
16.0 * 1.875e-5 = 3.0e-4 m
```

Both executables use the shared sampled-initialization helper with:

```text
sample_count = 20000
initialization_seed = 20260199
```

The driver generates one deterministic finite Monte Carlo realization from the
configured distribution. The empirical raw moments of that same finite
realization initialize the CQMOM case, and the Monte Carlo companion executable
starts from the same sampled realization. Therefore, the initial CQMOM and
Monte Carlo moment rows agree to roundoff; analytical moments are not used for
the production CQMOM initialization in this example.

## CQMOM Configuration

The default `configuration.qmom` uses:

```text
[phase_space]
names = r, w, T
node_counts = 2, 2, 2
conditioning_order = 0, 1, 2
adaptive_branch_counts = true

[time_integration.node_advection]
scheme = adaptive_rk23

[time_integration.node_advection.rk23]
relative_tolerance = 1.0e-6
absolute_tolerance = 1.0e-9
initial_step = 1.0e-4
max_step = 1.0e-2
min_step = 1.0e-9
```

The parsed case reports lean moment storage and the automatic Vandermonde
solver policy. The `node_counts = 2, 2, 2` setting is the current default
configuration, not an optimality claim.

The moments written to the result files are determined by the configured CQMOM
case through:

```cpp
qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy()).required_indices
```

There is no hardcoded output-moment list in the implementation. Changing the
configured CQMOM hierarchy changes the written moment set.

## Physical Model

The continuous physical callback receives coordinates ordered as `{r, w, T}`
and returns rates ordered as:

```text
{dr/dt, dw/dt, dT/dt}
```

The implemented droplet model includes:

- buoyancy-corrected gravity;
- aerodynamic drag;
- a Reynolds-number-dependent drag coefficient;
- deformation correction using constant surface tension;
- condensation-driven radius change;
- Sherwood-number dependence of mass transfer;
- sensible and latent heat transfer;
- a base Nusselt correlation;
- temperature-dependent saturation pressure, water-vapor diffusivity, and
  latent heat where implemented.

Surface tension is constant in this example. The implementation deliberately
uses the base Nusselt correlation and omits the evaporation-oriented wet-bulb
correction from the reference paper.

The fixed ambient conditions in the implementation are:

```text
ambient.temperature = 313.15 K
ambient.relative_humidity = 1.0
ambient.vertical_air_velocity = 0.0 m/s
ambient.gravity = 9.81 m/s^2
```

The fixed liquid-water/model parameters are:

```text
water_density = 997.0 kg/m^3
water_heat_capacity = 4181.0 J/(kg K)
surface_tension = 0.072 N/m
water_molar_mass = 18.0153e-3 kg/mol
gas_constant = 8.3145 J/(mol K)
water_critical_temperature = 647.096 K
```

## Numerical Evolution

The time controls are local constants in both drivers:

```text
initial_time = 0.0
final_time = 1.0
time_step = 0.01
```

The output directory, enabled files, and output interval are read from
`[output]` in `configuration.qmom`.

### CQMOM Executable

`cqmom_condensing_droplets_example`:

1. samples the initial population;
2. computes empirical raw moments from that finite sample;
3. reconstructs the initial 3D CQMOM quadrature;
4. advances quadrature nodes directly in phase space;
5. calls `qmomkit::advance_cqmom_continuous(...)`;
6. uses adaptive RK23 through the configured node-advection options.

Quadrature weights remain unchanged because the evolution is continuous only.

### Monte Carlo Executable

`cqmom_condensing_droplets_mc_example`:

1. generates the same deterministic initial sample population;
2. advances every sampled droplet with the same physical callback;
3. calls `qmomkit::monte_carlo::evolve_ensemble(...)`;
4. uses RK23 options derived from the configured node-advection settings;
5. writes the same case-defined moment set as the CQMOM executable.

The Monte Carlo companion provides a sampled population evolved with the same
physical model and initial finite realization. Differences after the initial
time arise from the finite CQMOM representation and numerical evolution, not
from different sampled initial populations. Exact agreement after evolution is
not expected, and the default `2 x 2 x 2` topology is not claimed to be
converged or optimal.

## Build

This example is a standalone installed-package consumer. It requires both
QMOMKit and the optional QMOMKitMonteCarlo companion package.

From the repository root, install both packages into the shared
repository-local `install/` prefix:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build build
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

Then configure and build this example from its directory:

```bash
cd examples/cqmom_condensing_droplets

cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

This builds:

- `cqmom_condensing_droplets_example`
- `cqmom_condensing_droplets_mc_example`

CMake copies `configuration.qmom` into the example build directory.

## Run

Run from the example build directory so generated output stays there:

```bash
cd build
./cqmom_condensing_droplets_example
./cqmom_condensing_droplets_mc_example
```

The CQMOM executable prints phase-space and conditioning information,
configured node counts, transport and RK23 settings, sample count and seed,
initial and final moments, final quadrature, and output paths.

The Monte Carlo executable prints phase-space and RK23 information, sample
count and seed, initial and final moments, final ensemble size, and its output
path.

A representative startup excerpt is:

```text
phase-space dimensions: 3
phase-space ordering: (r,w,T)
node_counts = [2,2,2]
Node advection settings:
  scheme = adaptive_rk23
sample_count = 20000
initialization_seed = 20260199
```

## Result Files

The default output configuration writes:

- `output/initial_moments.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`
- `output/moment_history_mc.txt`

`initial_moments.txt` contains the case-defined initial CQMOM moments.
`moment_history.txt` contains the CQMOM moment evolution.
`quadrature_history.txt` contains terminal CQMOM nodes, weights, abscissas,
and conditioning paths. `moment_history_mc.txt` contains the matching Monte
Carlo moment evolution.

The CQMOM and Monte Carlo histories use the same case-defined moment columns.
Labels, ordering, and output times match. The initial rows agree to roundoff
because both derive from the same finite sample realization.

## Expected Behavior

For the default case:

- total population weight remains constant;
- droplet radii increase because of condensation;
- droplet temperatures approach the ambient temperature;
- vertical velocities evolve under gravity and drag;
- smaller and larger droplets exhibit different settling dynamics;
- all radii remain positive;
- default radii remain below the intended breakup range.

These are qualitative expectations for the shipped configuration, not a
quantitative accuracy or convergence claim.

## Limitations

This example covers isolated droplets in prescribed uniform ambient
conditions. It does not include feedback of condensation on the gas phase,
droplet-droplet interactions, breakup, coalescence, physical-space transport,
or a zero-radius population boundary treatment.

The heat-transfer model uses the base Nusselt correlation without the
wet-bulb correction. The default case is selected to remain in the
condensation regime.
