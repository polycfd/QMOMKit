# Five-Dimensional CQMOM With Non-Spherical Particle Dynamics And Binary Fragmentation

This example advances the same five-dimensional non-spherical particle model in
two ways: deterministic CQMOM transport and direct Monte Carlo ensemble
evolution. Both executables use the same sampled initial realization, continuous
particle dynamics, binary-fragmentation parameters, and configured output
controls. The example is a workflow and model demonstration, not a calibrated or
validated physical application.

The workflow combines continuous phase-space dynamics, a discontinuous binary
fragmentation operator, multidimensional CQMOM reconstruction, adaptive
branch-specific node counts, moment output, and quadrature output.

## Phase Space And Moments

The configured phase-space coordinate vector is:

```text
xi = (d, chi, u, v, w)
```

where:

- `d` is equivalent particle diameter `[m]`;
- `chi` is particle aspect ratio `[-]`;
- `u`, `v`, and `w` are particle velocity components `[m/s]`.

This order is used consistently by `configuration.qmom`, the continuous and
discontinuous callbacks, raw-moment exponents, sampled ensembles, and output
columns.

A raw moment is indexed by:

```text
alpha = (a_d, a_chi, a_u, a_v, a_w)
```

and represents:

```text
M_alpha =
int d^a_d chi^a_chi u^a_u v^a_v w^a_w
    f(d, chi, u, v, w) dxi
```

QMOMKit uses unnormalized raw moments: `M_(0,0,0,0,0)` is the total number
density, not one unless the configured distribution has unit number density.

## Continuous Particle Dynamics

The continuous operator is defined in `particle_model.hpp` and
`particle_model.cpp`. It advances the five phase-space coordinates with:

```text
dd/dt   = -k_attr |U_f(t) - U_p|^2
dchi/dt = (1 - chi) / tau_chi
du/dt   = lambda_drag (U_f,x(t) - u)
dv/dt   = lambda_drag (U_f,y(t) - v)
dw/dt   = lambda_drag (U_f,z(t) - w) - g
```

where `U_p = (u, v, w)`. The fluid velocity is prescribed as:

```text
U_f,x(t) = A_x cos(omega t)
U_f,y(t) = A_y sin(omega t)
U_f,z(t) = A_z cos(omega t)
```

The drag rate is:

```text
lambda_drag = 18 mu_f K(chi, Re) / (rho_p d^2)
Re          = rho_f |U_f(t) - U_p| d / mu_f
```

`K(chi, Re)` is the implemented non-spherical, oblate-spheroid drag correction,
with a separate spherical-limit branch near `chi = 1`.

Current example-local continuous-model constants are:

| Symbol | Source field | Value | Units |
|---|---:|---:|---|
| `rho_f` | `fluid_density` | `1.2` | `kg/m^3` |
| `mu_f` | `fluid_dynamic_viscosity` | `1.8e-5` | `Pa s` |
| `rho_p` | `particle_density` | `2000.0` | `kg/m^3` |
| `g` | `gravity` | `9.81` | `m/s^2` |
| `A_x` | `fluid_velocity_amplitude_x` | `2.0` | `m/s` |
| `A_y` | `fluid_velocity_amplitude_y` | `2.0` | `m/s` |
| `A_z` | `fluid_velocity_amplitude_z` | `0.0` | `m/s` |
| `omega` | `oscillation_angular_frequency` | `1.0` | `1/s` |
| `k_attr` | `attrition_rate` | `1.0e-4` | `s/m` |
| `tau_chi` | `aspect_ratio_relaxation_time` | `2.0` | `s` |

Diameter decreases with relative-velocity-squared attrition. The aspect ratio
relaxes toward one under the continuous operator. Velocity components are driven
toward the imposed fluid velocity by drag, while gravity acts only in the
vertical `w` equation.

## Binary Fragmentation Model

The discontinuous operator is defined in `breakage_model.hpp` and
`breakage_model.cpp`. It represents binary fragmentation: one parent particle is
removed and two daughter particles are created.

### Breakage Rate

The event rate is:

```text
beta(d) = k_frag d^p_frag,  d >= d_min
beta(d) = 0,                d < d_min
```

Current example-local fragmentation constants are:

| Symbol | Source field | Value | Units |
|---|---:|---:|---|
| `k_frag` | `fragmentation_rate_constant` | `1.0e8` | `m^-2 s^-1` for `p_frag = 2` |
| `p_frag` | `fragmentation_power` | `2.0` | `[-]` |
| `d_min` | `minimum_fragmentation_diameter` | `1.0e-6` | `m` |
| `c_chi` | `shape_relaxation_factor` | `0.5` | `[-]` |
| `sigma_kick` | `recoil_velocity_standard_deviation` | `0.5` | `m/s` |

In general, the units of `k_frag` depend on `p_frag` so that `beta` has units
of `s^-1`.

### Daughter Properties

Equal-volume binary fragmentation sets:

```text
d_daughter = d_parent / 2^(1/3)
```

so that:

```text
2 d_daughter^3 = d_parent^3
```

The daughter aspect ratio is:

```text
chi_daughter = c_chi chi_parent
```

With the current `c_chi = 0.5`, each daughter carries half the parent aspect
ratio. This is an example-local shape rule for the fragmentation event; it is
separate from the continuous aspect-ratio relaxation equation.

### Gaussian Recoil

Daughter velocities are centered on the parent velocity with Gaussian recoil.
In the Monte Carlo executable, one recoil vector is sampled for each event and
applied with opposite signs:

```text
U_daughter,1 = U_parent + xi
U_daughter,2 = U_parent - xi
xi_j ~ N(0, sigma_kick^2)
```

The paired kicks give zero net recoil for each daughter pair.

In the CQMOM source, recoil is represented analytically through shifted Gaussian
raw moments rather than sampled kicks:

```text
G_k(x, sigma) = E[(x + xi)^k]
```

The finite Gaussian moment expansion is evaluated for each requested velocity
exponent.

### Raw-Moment Source

For a parent monomial `P_alpha` and the expected one-daughter monomial
`D_alpha`, the deterministic CQMOM source is:

```text
S_alpha = sum_i w_i beta(d_i) [2 D_alpha,i - P_alpha,i]
```

The sum runs over active terminal nodes of the reconstructed stage quadrature.
Here `w_i` is the active terminal-node weight, `-P_alpha` removes one parent,
and `2 D_alpha` creates two statistically identical daughters.

The direct Monte Carlo source applies the same physical event model with
stochastic waiting times. Daughters re-enter the remaining source-stage interval,
so multiple fragmentation generations can occur during one discontinuous stage.

## Fragmentation-Stage Balances

The binary fragmentation operator has the following qualitative balances:

- particle number increases by one per event;
- equal-volume daughters conserve the diameter-cubed moment;
- paired recoil conserves diameter-cubed weighted momentum;
- Gaussian recoil increases velocity second moments;
- aspect-ratio moments change according to the factor `c_chi`.

For the discontinuous operator alone:

```text
S_(3,0,0,0,0) = 0
```

and paired recoil gives:

```text
S_(3,0,1,0,0) = 0
S_(3,0,0,1,0) = 0
S_(3,0,0,0,1) = 0
```

These are fragmentation-stage balances. They do not imply that the full
continuous-plus-discontinuous model conserves the same moments: attrition, drag,
gravity, and continuous aspect-ratio relaxation also contribute to moment
evolution.

## CQMOM And Monte Carlo Correspondence

| Feature | CQMOM executable | Monte Carlo executable |
|---|---|---|
| Target | `cqmom_non_spherical_particles_example` | `cqmom_non_spherical_particles_mc_example` |
| State representation | Recursive conditional quadrature state | Weighted particle ensemble |
| Continuous dynamics | Quadrature-node transport | Particle trajectory evolution |
| Fragmentation | Deterministic raw-moment source | Stochastic branching events |
| Recoil | Analytical Gaussian moments | Sampled paired kicks |
| Output | Moments and quadrature | Dense ensemble moments |

Both paths use the same initial sampled realization, phase-space ordering,
continuous particle model, fragmentation parameters, macro time step, final
time, and output interval. The CQMOM result is a finite-node deterministic
moment representation, while the Monte Carlo result is one fixed-seed stochastic
ensemble realization; exact one-to-one agreement is not implied.

## Initialization

Initialization is intentionally sample-based in both drivers:

1. `configuration.qmom` defines independent marginal distributions.
2. The example samples `10000` particles with `initialization_seed = 12345`.
3. Dense raw moments are computed from that same finite ensemble.
4. The CQMOM executable replaces the distribution-based `[initial]` section with
   sampled `[moments]` to initialize the CQMOM hierarchy.
5. The Monte Carlo executable evolves the sampled ensemble directly.

This means both approaches start from the same finite statistical realization
rather than from analytical moments in one path and sampled moments in the
other.

The sampled dense tensor extents are controlled by:

```text
maximum_orders = 3, 5, 3, 3, 3
```

The CQMOM hierarchy requests three branch-local nodes for the aspect-ratio
coordinate:

```text
node_counts = 2, 3, 2, 2, 2
```

The `chi` marginal is beta-distributed and physically bounded on `[0, 1]`.
`configuration.qmom` declares that interval in phase-space order:

```text
support_lower = none, 0.0, none, none, none
support_upper = none, 1.0, none, none, none
```

QMOMKit uses the interval as a bounded-support realizability criterion during
branch-local CQMOM reconstruction. It does not clip reconstructed aspect-ratio
nodes.

Monte Carlo sampled moments use the QMOMKit number-density convention:

```text
M_alpha = n0 / N sum_p product_r xi_(r,p)^alpha_r
```

with `n0 = 1.0` in the supplied configuration.

## Numerical Workflow

Continuous dynamics and discontinuous fragmentation are advanced through
operator splitting. In the CQMOM executable, the splitting scheme and continuous
node-advection settings are selected in `configuration.qmom`, not by editing the
driver source.

The supplied configuration selects:

```text
[time_integration]
scheme = strang

[time_integration.node_advection]
scheme = adaptive_rk23
```

with RK23 tolerances and step controls specified in
`[time_integration.node_advection.rk23]`. The CQMOM executable delegates the
configured transport step to QMOMKit.

The Monte Carlo executable applies the corresponding current sequence explicitly
to the ensemble: a continuous half-step, stochastic fragmentation over the full
macro step, and a second continuous half-step. Its fragmentation stage may
increase the ensemble size.

The executable time controls are fixed in the drivers:

```text
final_time = 1.0
time_step  = 0.01
```

## Configuration Summary

`configuration.qmom` controls the numerical CQMOM and output setup:

| Setting | Current value | Meaning |
|---|---|---|
| `phase_space.names` | `d, chi, u, v, w` | phase-space coordinate order |
| `node_counts` | `2, 2, 2, 2, 2` | requested CQMOM nodes per phase-space dimension |
| `conditioning_order` | `0, 1, 2, 3, 4` | recursive CQMOM conditioning order |
| `adaptive_branch_counts` | `true` | enables branch-specific active node counts |
| `type` | `independent_product` | independent marginal initialization |
| `number_density` | `1.0` | zeroth raw moment of the initial distribution |
| `maximum_orders` | `3, 3, 3, 3, 3` | dense sampled moment tensor extents |
| `scheme` | `strang` | configured CQMOM operator-splitting scheme |
| `node_advection.scheme` | `adaptive_rk23` | configured CQMOM continuous node-advection integrator |
| `output.directory` | `output` | output directory under the run directory |
| `output.interval` | `0.1` | history output interval |

The supplied marginal distributions are lognormal for `d`, beta for `chi`, and
normal for each velocity component. Continuous particle constants live in
`particle_model.hpp`; fragmentation constants live in `breakage_model.hpp`.

## Executables

- `cqmom_non_spherical_particles_example` samples the configured initialization,
  initializes CQMOM from the sampled raw moments, advances the CQMOM transport
  state, and writes moment and quadrature files.
- `cqmom_non_spherical_particles_mc_example` samples the same initialization,
  evolves the weighted Monte Carlo ensemble directly, and writes dense ensemble
  moment histories.

The two executables are companion workflows for the same model. Use their output
as a qualitative comparison unless a separate validation study defines a
quantitative acceptance criterion.

## Build And Run

This standalone example consumes installed QMOMKit and QMOMKitMonteCarlo
packages. From the repository root, install both packages into the shared local
prefix:

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

Configure and build the example from its directory:

```bash
cd examples/cqmom_non-spherical-particles
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
```

The build copies `configuration.qmom` into `build/`. Run from that directory:

```bash
cd build
./cqmom_non_spherical_particles_example
./cqmom_non_spherical_particles_mc_example
```

## Output

With the supplied output settings, the CQMOM executable writes:

```text
output/initial_moments.txt
output/moment_history.txt
output/initial_quadrature.txt
output/quadrature_history.txt
```

The CQMOM moment files use the complete hierarchy-required moment layout. The
quadrature files contain active CQMOM nodes, their conditioning paths, weights,
and abscissas.

The Monte Carlo executable writes:

```text
output/moment_history_mc.txt
```

The Monte Carlo moment file uses dense tensor order from `maximum_orders`.

Terminal summaries print only a small common subset of moments. That subset is
stable-sorted by total polynomial degree for readability, with the CQMOM
hierarchy order used as the within-degree tie-break. This terminal ordering does
not change the file-column ordering.

## Expected Qualitative Behavior

Under active fragmentation, particle number increases and characteristic
diameters decrease. Fragmentation alone conserves diameter-cubed content, while
the continuous attrition term changes diameter moments during the full coupled
run. The daughter aspect-ratio scaling changes aspect-ratio moments, Gaussian
recoil broadens velocity distributions, and continuous drag and gravity modify
velocity moments independently of fragmentation.
