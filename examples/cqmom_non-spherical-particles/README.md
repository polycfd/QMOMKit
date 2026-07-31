# Five-dimensional non-spherical particle CQMOM example

## Purpose

This standalone example advances a five-dimensional non-spherical particle
population with CQMOM and a direct Monte Carlo comparison. The phase-space
model includes drag in an oscillating carrier flow, gravity, regularized
diameter attrition, aspect-ratio relaxation, binary fragmentation, and recoil
velocity spread.

It includes three executables:

- `cqmom_non_spherical_particles_example`: the normal bounded-support CQMOM
  simulation.
- `cqmom_non_spherical_particles_unbounded_diagnostic`: an intentionally
  unbounded diagnostic that explains a conditional support-loss event.
- `cqmom_non_spherical_particles_mc_example`: a direct stochastic Monte Carlo
  comparison for the same particle model and sampled initialization.

The bounded CQMOM executable is the normal physical workflow. The unbounded
diagnostic is not a different physical model; it omits the aspect-ratio support
only to show why finite conditional moments can be Stieltjes-realizable while
failing bounded Hausdorff realizability on `[0,1]`.

## Phase-space model

The particle coordinate is:

```text
xi = (d, chi, u, v, w)
```

where:

- `d` is particle diameter, in meters.
- `chi` is particle aspect ratio, dimensionless.
- `u`, `v`, and `w` are particle-velocity components, in meters per second.

Continuous transport uses the phase-space velocity:

```text
dxi/dt = (dd/dt, dchi/dt, du/dt, dv/dt, dw/dt)
```

Binary fragmentation is treated separately as a discontinuous population source
acting on raw moments or, in the Monte Carlo executable, on sampled particles.
Raw moments are unnormalized number-density moments, so `M[0,0,0,0,0]` is the
total number density.

## Example-local physical parameters

These are the model values compiled into this example. They are not claimed to
be universal calibrated material constants.

| Symbol | Code member | Value | Units | Role |
|---|---|---:|---|---|
| `rho_f` | `fluid_density` | `1.2` | kg/m^3 | carrier-fluid density |
| `rho_p` | `particle_density` | `2000.0` | kg/m^3 | particle material density |
| `mu_f` | `fluid_dynamic_viscosity` | `1.8e-5` | Pa s | carrier-fluid dynamic viscosity |
| `g` | `gravity` | `9.81` | m/s^2 | downward gravitational acceleration |
| `A_x` | `fluid_velocity_amplitude_x` | `2.0` | m/s | x carrier-velocity amplitude |
| `A_y` | `fluid_velocity_amplitude_y` | `2.0` | m/s | y carrier-velocity amplitude |
| `A_z` | `fluid_velocity_amplitude_z` | `0.0` | m/s | z carrier-velocity amplitude |
| `omega` | `oscillation_angular_frequency` | `1.0` | 1/s | carrier-flow angular frequency |
| `k_attr` | `attrition_rate` | `1.0e-6` | s/m | diameter attrition coefficient |
| `U_sat` | `attrition_velocity_scale` | `5.0` | m/s | velocity scale for attrition saturation |
| `d_min` | `minimum_attrition_diameter` | `0.0` | m | lower scale in attrition taper |
| `d_transition` | `attrition_transition_diameter` | `1.0e-5` | m | attrition taper transition length |
| `tau_chi` | `aspect_ratio_relaxation_time` | `2.0` | s | continuous relaxation time toward `chi = 1` |
| `k_frag` | `fragmentation_rate_constant` | `1.0e8` | m^-2 s^-1 | coefficient in `beta(d)` for `p_frag = 2` |
| `p_frag` | `fragmentation_power` | `2.0` | dimensionless | diameter exponent in fragmentation rate |
| `d_frag,min` | `minimum_fragmentation_diameter` | `1.0e-6` | m | fragmentation cutoff diameter |
| `f_chi` | `shape_relaxation_factor` | `0.5` | dimensionless | daughter aspect-ratio multiplier |
| `sigma_recoil` | `recoil_velocity_standard_deviation` | `0.5` | m/s | standard deviation of each recoil component |

## Oscillating carrier-fluid velocity

The implemented carrier-fluid velocity is spatially uniform and time
dependent:

```text
U_f,x(t) = A_x cos(omega t)
U_f,y(t) = A_y sin(omega t)
U_f,z(t) = A_z cos(omega t)
```

The particle velocity is `U_p = (u, v, w)`. The relative velocity driving both
drag and attrition is:

```text
U_rel = U_f - U_p

U_rel,x = U_f,x - u
U_rel,y = U_f,y - v
U_rel,z = U_f,z - w

|U_rel| =
sqrt(
    (U_f,x - u)^2
  + (U_f,y - v)^2
  + (U_f,z - w)^2
)
```

## Reynolds number

The particle Reynolds number is:

```text
Re = rho_f |U_rel| d / mu_f
```

Here `rho_f` has units kg/m^3, `|U_rel|` has units m/s, `d` has units m, and
`mu_f` has units kg/(m s), so `Re` is dimensionless.

## Non-spherical drag model

The drag response uses an aspect-ratio-dependent correction for an oblate
spheroidal particle. For physical aspect ratios in `(0,1]`, define:

```text
e = sqrt(1 - chi^2)
```

If `|1 - chi| < 0.001`, the implementation uses the spherical-limit value:

```text
K0(chi) = 1
```

Otherwise it evaluates:

```text
K0(chi) =
    (8 / (3 chi^(1/3)))
    /
    (
        2 chi / (1 - chi^2)
      + 2 (1 - 2 chi^2) atan(e / chi) / (1 - chi^2)^(3/2)
    )
```

The complete shape factor includes empirical Reynolds-number corrections:

```text
K(chi, Re) =
    K0(chi)
  + 0.15 chi^95.91 Re^0.687
  + 0.2927 (1 - chi)^0.4374 Re^0.7512
```

The first term supplies the oblate-spheroid creeping-flow correction. The
second and third terms increase the drag response with Reynolds number and make
the correction depend on particle shape. The formula is evaluated directly by
the example callback; the callback does not clip `d` or `chi`.

## Drag and gravity

The implemented drag-rate coefficient is:

```text
lambda_drag =
    18 mu_f K(chi, Re)
    /
    (rho_p d^2)
```

It has units 1/s. This is algebraically equivalent to evaluating the drag force

```text
F_d = 1/2 rho_f A (24/Re) K Re |U_rel| U_rel
```

and dividing by the spherical particle mass. The velocity equations are:

```text
du/dt = lambda_drag (U_f,x - u)
dv/dt = lambda_drag (U_f,y - v)
dw/dt = lambda_drag (U_f,z - w) - g
```

Drag relaxes each particle-velocity component toward the carrier-fluid
velocity. Because `lambda_drag` scales like `1/d^2`, smaller particles respond
more rapidly to drag. Gravity acts only in the `w` equation, with positive `g`
subtracted from `dw/dt`; in this sign convention gravity accelerates particles
in the negative `w` direction.

## Regularized diameter attrition

Diameter attrition uses the same relative speed as drag, but with two
regularizations. First, the velocity-squared response is saturated:

```text
U_attr^2 =
    |U_rel|^2
    /
    (1 + |U_rel|^2 / U_sat^2)
```

Second, the recession rate is tapered near the configured lower diameter scale:

```text
dd/dt =
    -k_attr U_attr^2
    tanh((d - d_min) / d_transition)
```

For ordinary relative velocities, `U_attr^2` approaches the original quadratic
dependence `|U_rel|^2`. At extreme relative velocity, `U_attr^2` is bounded by
the saturation scale. Far above `d_min`, the `tanh` factor approaches one and
the original recession behavior is recovered. Near `d_min`, the taper makes the
boundary asymptotically attracting and avoids finite-time crossing.

The current minimum attrition diameter is `d_min = 0`. This regularization
changes the attrition law itself; it is not state clipping, and it does not
remove extreme velocity quadrature nodes.

## Aspect-ratio evolution

The continuous aspect-ratio equation is:

```text
dchi/dt = (1 - chi) / tau_chi
```

The relaxation time is `tau_chi = 2.0 s`. This continuous term drives `chi`
toward the spherical value `chi = 1`. It is distinct from the discontinuous
daughter aspect-ratio transformation used during fragmentation.

## Binary fragmentation rate

Binary fragmentation uses the diameter-dependent event rate:

```text
beta(d) =
    k_frag d^p_frag,  d >= d_frag,min
    0,                d < d_frag,min
```

with `k_frag = 1.0e8 m^-2 s^-1`, `p_frag = 2.0`, and
`d_frag,min = 1.0e-6 m`. The minimum fragmentation diameter stops further
fragmentation events for smaller particles, but it is not a lower bound on
particle diameter under continuous attrition.

## Daughter diameter and shape

One parent produces two daughters. Daughter diameter is:

```text
d_daughter = 2^(-1/3) d_parent
```

so:

```text
2 d_daughter^3 = d_parent^3
```

Binary fragmentation therefore conserves diameter-cubed content for each event.
The daughter aspect ratio is:

```text
chi_daughter = f_chi chi_parent
```

with `f_chi = 0.5`. This discontinuous multiplicative shape change is separate
from the continuous relaxation toward `chi = 1`.

## Recoil velocity

Each Monte Carlo fragmentation event samples one centered Gaussian recoil
vector:

```text
zeta = (zeta_u, zeta_v, zeta_w)
zeta_j ~ N(0, sigma_recoil^2)
```

with `sigma_recoil = 0.5 m/s`. The paired daughter velocities are:

```text
U_daughter,1 = U_parent + zeta
U_daughter,2 = U_parent - zeta
```

Equal and opposite kicks give zero net pair recoil. Because the daughters have
equal volume and the two daughter volumes sum to the parent volume, this paired
recoil preserves mass-weighted momentum for each fragmentation event. It also
increases velocity spread and velocity second moments.

The deterministic CQMOM source does not sample recoil. Instead, daughter
velocity moments use analytical shifted Gaussian moments:

```text
G_k(x, sigma) = E[(x + zeta)^k]
```

Only even centered-Gaussian powers contribute to the binomial expansion.
Representative low-order moments are:

```text
G_0 = 1
G_1 = x
G_2 = x^2 + sigma^2
G_3 = x^3 + 3 x sigma^2
```

## Fragmentation moment source

For a raw-moment multi-index `alpha = (a_d, a_chi, a_u, a_v, a_w)`, the parent
monomial is:

```text
P_alpha =
    d^a_d chi^a_chi u^a_u v^a_v w^a_w
```

The expected one-daughter monomial is:

```text
D_alpha =
    d_daughter^a_d
    chi_daughter^a_chi
    G_a_u(u, sigma_recoil)
    G_a_v(v, sigma_recoil)
    G_a_w(w, sigma_recoil)
```

The deterministic CQMOM source sums terminal-quadrature contributions as:

```text
S_alpha =
    sum_i w_i beta(d_i) (2 D_alpha,i - P_alpha,i)
```

The `-P_alpha` term removes one parent, `2 D_alpha` creates two statistically
identical daughters, `w_i` supplies the parent population weight, and the
daughter velocity factors include the Gaussian recoil expectation.

The Monte Carlo source applies the same event model stochastically. Within each
discontinuous-stage interval, daughter lineages re-enter with the remaining
stage time, so repeated fragmentation in one source interval is allowed.

## Physical balances

For the discontinuous fragmentation operator alone:

- particle number increases by one per binary event;
- diameter-cubed content is conserved by `2 d_daughter^3 = d_parent^3`;
- paired recoil preserves mass-weighted momentum eventwise for equal-volume
  daughters;
- recoil increases velocity spread and velocity second moments.

These are fragmentation-only balances. The complete model also includes
continuous drag, gravity, diameter attrition, and aspect-ratio relaxation.
Continuous drag exchanges momentum with the carrier flow, gravity accelerates
the vertical velocity component, and attrition changes particle volume and
diameter moments.

## Initialization shared by CQMOM and Monte Carlo

The CQMOM drivers parse the independent-product initialization, sample `10000`
particles with `initialization_seed = 12345`, and replace the distributional
initial condition with raw moments of that same finite ensemble. This makes the
bounded and unbounded CQMOM initial states identical except for the support
declaration.

The Monte Carlo executable uses the same sampled ensemble directly. Each
sample carries `n0 / N`, so the ensemble represents unnormalized raw moments of
the particle number-density distribution.

## Bounded-support CQMOM workflow

`configuration.qmom` is the normal bounded configuration. It declares the
physical aspect-ratio support in phase-space order:

```text
support_lower = none, 0.0, none, none, none
support_upper = none, 1.0, none, none, none
```

QMOMKit uses this interval as a support-aware recursive branch acceptance
check. The bounded workflow advances the full physical model, completes the
full CQMOM transport run, writes standard moment and quadrature histories, and
prints only concise support topology events.

## Unbounded diagnostic workflow

`configuration_unbounded.qmom` intentionally omits only the chi support
declaration while preserving the same phase-space hierarchy, initialization,
physical callbacks, integration settings, tolerances, and source model. Its
output directory is `output/unbounded_diagnostic`.

The diagnostic advances the same physics to the first source boundary where
ordinary unconstrained post-source inversion accepts a finite terminal node
with `chi < 0` or `chi > 1`. It then writes one self-contained Markdown report
and returns success before the next continuous stage can obscure the origin of
the event.

## Why both CQMOM workflows are included

Together, the bounded and unbounded workflows demonstrate:

- the practical five-dimensional non-spherical-particle model;
- conditional support loss in finite multidimensional moment closures;
- the distinction between Stieltjes realizability and bounded Hausdorff
  realizability;
- why positive classical Hankel-Hadamard determinants do not guarantee
  `chi <= 1`;
- how support-aware adaptive reconstruction resolves the issue without
  clipping node coordinates.

The unbounded workflow is not a different physical model. It is included to
show how a finite conditional moment sequence can leave the physical support
even though the fragmentation daughter law preserves `chi <= 1` when the parent
is inside `[0,1]`. The bounded workflow demonstrates the generic
support-aware response.

## Conditional support-loss mechanism

After the discontinuous source update, a branch-local conditional chi sequence
can remain realizable on the nonnegative real line while failing feasibility on
the bounded physical interval `[0,1]`. This is a higher-order conditional
moment representation/support issue. It is not random sampling noise, a
negative-weight quadrature, a Vandermonde failure, a Wheeler defect, or a
physical daughter rule that creates `chi > 1`.

## Classical Hankel-Hadamard and interval-localizing conditions

For a branch sequence `m0 ... m5`, the diagnostic reports classical
Hankel-Hadamard determinants:

```text
H(k,l) = [m_(k+i+j)] for i,j = 0,...,l
Delta(k,l) = det(H(k,l))
```

with `k = 0` and `k = 1`. These checks establish real-line and
nonnegative-support moment conditions.

For interval `[a,b]`, bounded support also requires the interval-localizing
matrices:

```text
Delta_lower(l) = det([m_(i+j+1) - a*m_(i+j)])
Delta_upper(l) = det([b*m_(i+j) - m_(i+j+1)])
```

For `[0,1]`, the upper condition is `det([m_(i+j) - m_(i+j+1)])`. The unbounded
diagnostic shows that the classical determinants can remain positive while the
upper localizing matrix fails. Wheeler can still be algebraically accurate in
this situation: it reconstructs the supplied unbounded-support sequence, but
one returned chi node lies outside the finite physical interval.

## Support-aware adaptive reconstruction

The bounded workflow rejects the support-infeasible six-moment three-node chi
candidate, accepts a support-valid lower-count branch, rebuilds descendants,
and canonicalizes the accepted state from the accepted quadrature.

This is not clipping. It is not projection of all six moments onto a bounded
moment cone. It is an adaptive reduction of represented branch order:
higher-order source-candidate information is replaced by the closure implied by
the accepted lower-count quadrature.

## Relationship between bounded and unbounded results

The unbounded workflow uses all six branch moments. Ordinary Wheeler
reconstruction returns three positive chi nodes, one above one, and reconstructs
all six supplied moments to roundoff.

The bounded workflow sees the same source-boundary event with the same sampled
realization and physical model. Because chi support is declared, it rejects the
infeasible three-node candidate and retains a lower-count support-valid branch.
The diagnostic report compares the unbounded accepted branch with the actual
bounded accepted branch from the bounded post-source quadrature.

## Executables

| Executable | Purpose |
|---|---|
| `cqmom_non_spherical_particles_example` | Normal bounded-support CQMOM simulation. |
| `cqmom_non_spherical_particles_unbounded_diagnostic` | Early-stopping diagnostic for the first unconstrained chi support-loss event. |
| `cqmom_non_spherical_particles_mc_example` | Direct stochastic Monte Carlo comparison. |

## Configuration-selected numerical choices

Physical constants are compiled into the example model. Numerical choices are
selected through `configuration.qmom` and, for the diagnostic, the matching
`configuration_unbounded.qmom`. These include node counts, conditioning order,
coordinate support, operator-splitting scheme, node-advection scheme, RK23
tolerances, adaptive branch-count behavior, and output controls.

The supplied CQMOM hierarchy and initialization settings are:

```text
node_counts = 3, 3, 3, 3, 3
conditioning_order = 0, 1, 2, 3, 4
maximum_orders = 5, 5, 5, 5, 5
adaptive_branch_counts = true
```

The supplied CQMOM transport settings are:

```text
[time_integration]
scheme = cd

[time_integration.node_advection]
scheme = adaptive_rk23
```

with RK23 tolerances and step limits in:

```text
[time_integration.node_advection.rk23]
relative_tolerance = 1.0e-6
absolute_tolerance = 1.0e-9
initial_step = 1.0e-3
max_step = 1.0e-2
min_step = 1.0e-9
max_steps = 2000000
```

The CQMOM driver time controls are:

```text
initial_time = 0.0
final_time = 1.0
time_step = 0.005
```

The Monte Carlo driver currently uses `time_step = 0.01`.

## Configuration files

`configuration.qmom` is the normal bounded-support configuration.
`configuration_unbounded.qmom` intentionally removes only the chi support
declaration and redirects diagnostic output to `output/unbounded_diagnostic`.

## Build and run

Install QMOMKit and the optional Monte Carlo companion package into a shared
prefix, then configure this example as a standalone consumer:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build

cmake -S examples/cqmom_non-spherical-particles \
    -B examples/cqmom_non-spherical-particles/build \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build examples/cqmom_non-spherical-particles/build
```

Run from the example build directory:

```bash
cd examples/cqmom_non-spherical-particles/build
./cqmom_non_spherical_particles_example
./cqmom_non_spherical_particles_unbounded_diagnostic
./cqmom_non_spherical_particles_mc_example
```

## Output

The bounded CQMOM executable writes:

```text
output/initial_moments.txt
output/moment_history.txt
output/initial_quadrature.txt
output/quadrature_history.txt
```

The unbounded diagnostic writes:

```text
output/unbounded_diagnostic/report.md
```

The Monte Carlo executable writes:

```text
output/moment_history_mc.txt
```

## Expected behavior

The bounded CQMOM run remains within the declared chi support and reaches
`final_time = 1.0`. It may report concise support-driven active-count
reductions and support-infeasible growth rejections.

The unbounded diagnostic is expected to stop early after writing its report. In
the supplied configuration it detects a small positive-weight chi node above
one immediately after an accepted post-source reconstruction.

The Monte Carlo executable is a direct stochastic comparison. It may take
substantially longer than CQMOM because it starts with `10000` particles,
fragmentation grows the ensemble, every particle follows continuous dynamics,
and dense moments are recomputed from the evolving ensemble. A long MC runtime
is expected behavior, not a defect.

## Interpretation and limitations

This example is a deterministic numerical workflow and diagnostic
demonstrator, not a calibrated physical application. Coordinate names are
example metadata; QMOMKit itself does not assign physical units or enforce
application physics. The bounded support declaration is the correct way for
this model to communicate the physical aspect-ratio interval to generic CQMOM
reconstruction.
