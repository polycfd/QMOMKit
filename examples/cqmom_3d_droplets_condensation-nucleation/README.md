# 3D CQMOM droplets: condensation and nucleation

## Purpose and workflows

This standalone example evolves a spatially uniform polydisperse population of
isolated water droplets in the three-dimensional phase space
$\boldsymbol{\xi}=(r,w,T)$, where $r$ is droplet radius in metres, $w$ is
vertical droplet velocity in m/s, and $T$ is droplet temperature in K. Positive
$w$ is upward, so falling droplets have negative vertical velocity.

The model contains phase-space evolution only. Existing droplets move through
$(r,w,T)$ by condensation, drag, and heat transfer. During a finite interval,
nucleation adds droplets at one specified birth state through a phase-space
boundary contribution.

The example contains two executables:

- `cqmom_3d_droplets_condensation_nucleation_example`: initializes CQMOM from empirical raw
  moments of a deterministic finite sampled population and advances the CQMOM
  quadrature with Strang-split droplet advection and a nucleation boundary
  contribution.
- `cqmom_3d_droplets_condensation_nucleation_mc_example`: evolves the same deterministic finite
  sampled population directly, using midpoint population injection for the
  same nucleation model.

Both workflows use the same [configuration.qmom](configuration.qmom),
distribution parameters, sample count, initialization seed, coordinate order,
and droplet model callback. The coordinate names and distributions are example
choices; the configuration and sampled-initialization machinery are generic
QMOMKit mechanisms.

The CQMOM driver loads the normal case first, then replaces only its initial
raw moments with the complete sampled tensor through
`with_initial_moments(...)`. The loaded hierarchy, support, and reconstruction
policies are retained; a sampler is only one possible source of the replacement
tensor.

## QMOMKit focus

Read this as a three-dimensional CQMOM phase-space evolution example before
focusing on the droplet correlations: it combines adaptive CQMOM
reconstruction, condensation-driven node advection, and a finite-time CQMOM
boundary birth contribution for nucleation. The direct Monte Carlo companion
evolves the same sampled initialization, providing a stochastic particle-based
companion/reference realization for comparing the two representations.

## Phase-space definition and raw moments

The configured phase-space order is:
$$(r\ [\mathrm{m}], w\ [\mathrm{m/s}], T\ [\mathrm{K}])$$

The same positional order is used by [configuration.qmom](configuration.qmom),
samples, the phase-space velocity callback, CQMOM moment exponents, quadrature
abscissas, and output files.

The radius coordinate and its moments use SI metres. Reconstruction does not
depend on the units of the coordinates, so no example-local phase-space
rescaling is needed and the physical callback receives the SI radius directly.

For a number density function $n(r,w,T)$, the raw unnormalized moments are

```math
\begin{aligned}
M_{\alpha_r,\alpha_w,\alpha_T}
&=
\iiint
r^{\alpha_r}
w^{\alpha_w}
T^{\alpha_T}
n(r,w,T)
\,\mathrm{d}r\,\mathrm{d}w\,\mathrm{d}T.
\end{aligned}
```

Here $\alpha_r$, $\alpha_w$, and $\alpha_T$ are nonnegative integer
moment orders in radius, vertical velocity, and temperature. The zeroth moment
$M_{0,0,0}$ is the total number density.

For a finite sampled population, the implementation estimates the same raw
moments as

$$
\begin{aligned}
M_{\alpha_r,\alpha_w,\alpha_T}^{(N)}
&=
\frac{n_0}{N}
\sum_{p=1}^{N}
r_p^{\alpha_r}
w_p^{\alpha_w}
T_p^{\alpha_T},
\end{aligned}
$$

where $N$ is the sample count, $n_0$ is the configured total number
density, $(r_p,w_p,T_p)$ is sample $p$, and $n_0/N$ is the constant
statistical weight carried by each sample.

The direct Monte Carlo executable evaluates the configured selected moment
request with QMOMKitMonteCarlo's native selective estimator, passing the
current total number density. Request order and duplicate entries are
preserved; no example-local monomial or per-particle-weight implementation is
required.

## Droplet model

The advection callback receives coordinates ordered as $(r,w,T)$ and
returns
$$\left(\frac{\mathrm{d}r}{\mathrm{d}t}, \frac{\mathrm{d}w}{\mathrm{d}t}, \frac{\mathrm{d}T}{\mathrm{d}t} \right).$$
The droplet model is a simplified version of the model used by [Denner (2026)](https://www.doi.org/10.1016/j.ijmultiphaseflow.2026.105696) for water droplets in air.
The condensation-model constants are compiled into [droplet_model.cpp](./droplet_model.cpp); the
nucleation birth parameters are declared in [droplet_model.hpp](./droplet_model.hpp).


### Ambient and material properties

The prescribed ambient state is:

| Symbol | Code member | Value | Units | Meaning |
|---|---|---:|---|---|
| $T_\infty$ | `ambient.temperature` | `313.15` | K | ambient air temperature |
| $\phi$ | `ambient.relative_humidity` | `1.0` | dimensionless | ambient relative humidity |
| $w_\infty$ | `ambient.vertical_air_velocity` | `0.0` | m/s | vertical air velocity |
| $g$ | `ambient.gravity` | `9.81` | m/s^2 | gravitational acceleration |

The liquid-water and vapor constants are:

| Symbol | Code member | Value | Units | Meaning |
|---|---|---:|---|---|
| $\rho_l$ | `water_density` | `997.0` | kg/m^3 | liquid-water density |
| $c_l$ | `water_heat_capacity` | `4181.0` | J/(kg K) | liquid-water heat capacity |
| $\sigma$ | `surface_tension` | `0.072` | N/m | constant surface tension |
| $M_w$ | `water_molar_mass` | `18.0153e-3` | kg/mol | water molar mass |
| $R_g$ | `gas_constant` | `8.3145` | J/(mol K) | universal gas constant |
| $T_c$ | `water_critical_temperature` | `647.096` | K | water critical temperature |

The implementation also uses:

| Symbol | Code value | Meaning |
|---|---:|---|
| $T_0$ | `273.15` K | Celsius-to-kelvin offset and minimum liquid-water correlation temperature |
| $\mathrm{Re}_{\min}$ | `1.0e-8` | low-Reynolds-number drag branch threshold |

The moist-air properties are evaluated from a shared polynomial form. With
$\theta=T_\infty-T_0$,

$$
\begin{aligned}
P(\theta,\phi)
&=
a_1
+\theta(b_1+b_2\phi)
+\theta^2(c_1+c_2\phi).
\end{aligned}
$$

The implementation applies this form as:

| Property | Symbol | $a_1$ | $b_1$ | $b_2$ | $c_1$ | $c_2$ |
|---|---|---:|---:|---:|---:|---:|
| air density | $\rho_a$ | `1.2932` | `-4.6688e-3` | `-8.6965e-4` | `1.6165e-5` | `2.2437e-5` |
| dynamic viscosity | $\mu_a$ | `1.7192e-5` | `4.9467e-8` | `-2.2425e-9` | `-4.2412e-11` | `-3.2389e-10` |
| thermal conductivity | $k_a$ | `2.4031e-2` | `7.5904e-5` | `-3.1199e-6` | `-4.8653e-8` | `-1.3015e-7` |
| heat capacity | $c_a$ | `1006.4` | `3.9654e-2` | `2.0101` | `4.15208e-4` | `-1.673e-1` |

The correlations hold for ambient air temperatures in $[273.15,373.15]$ K and
droplet temperatures in $[273.15,T_c)$, which the shipped case stays within.

### Vapor and heat-transfer correlations

The saturation pressure correlation is

$$
\begin{aligned}
p_{\mathrm{sat}}(T)
&=
\exp\left(
34.494-\frac{4924.99}{T-36.05}
\right)
(T-168.15)^{-1.57}.
\end{aligned}
$$

The water-vapor diffusivity is evaluated at the ambient temperature:

$$
\begin{aligned}
D_v
&=
2.178\times 10^{-5}
\left(\frac{T_\infty}{T_0}\right)^{1.81}.
\end{aligned}
$$

The latent heat correlation uses

$$
\begin{aligned}
\theta_d &= \frac{T_c-T}{T}, \\
\theta_c &= \frac{T_c-T}{T_c},
\end{aligned}
$$

and

```math
\begin{aligned}
L(T)
&=
6.8477
\frac{
\theta_d^{0.325}
+1.7383\,\theta_c^{0.835}
}{
1+1.3913\,\theta_c^{0.835}
}
\frac{R_g T_c}{M_w}.
\end{aligned}
```

The relative vertical velocity is

$$
\begin{aligned}
w_{\mathrm{rel}} &= w_\infty-w.
\end{aligned}
$$

Positive $w_{\mathrm{rel}}$ means that the ambient air velocity is upward
relative to the droplet velocity. The Reynolds and Weber numbers are

$$
\begin{aligned}
\mathrm{Re}
&=
\frac{2r\rho_a |w_{\mathrm{rel}}|}{\mu_a}, \\
\mathrm{We}
&=
\frac{2r\rho_a w_{\mathrm{rel}}^2}{\sigma}.
\end{aligned}
$$

The transfer-number helper is [[Beard & Pruppacher,
1971](https://www.doi.org/10.1175/1520-0469(1971)028<1455:AWTIOT>2.0.CO;2)]

```math
\begin{aligned}
\Psi(\mathrm{Re},\Pi)
&=
\begin{cases}
2+0.316\,\mathrm{Re}\,\Pi^{2/3},
& \sqrt{\mathrm{Re}}\,\Pi^{1/3} < 1.4,
\\ 1.56+0.616\,\sqrt{\mathrm{Re}}\,\Pi^{1/3},
& \sqrt{\mathrm{Re}}\,\Pi^{1/3} \ge 1.4.
\end{cases}
\end{aligned}
```

The Sherwood number is

$$
\begin{aligned}
\mathrm{Sh}
&=
\Psi\left(
\mathrm{Re},
\frac{\mu_a}{\rho_a D_v}
\right),
\end{aligned}
$$

and the base Nusselt number used by this example is

$$
\begin{aligned}
\mathrm{Nu}
&=
\Psi\left(
\mathrm{Re},
\frac{\mu_a c_a}{k_a}
\right).
\end{aligned}
$$

No wet-bulb correction is applied to $\mathrm{Nu}$ in this implementation.

### Condensation radius rate

The implemented radius rate is
```math
\begin{aligned}
\frac{\mathrm{d}r}{\mathrm{d}t}
&=
\frac{\mathrm{Sh}\,D_v\,M_w}{2r\rho_l R_g} \ \left(\phi\frac{p_{\mathrm{sat}}(T_\infty)}{T_\infty} - \frac{p_{\mathrm{sat}}(T)}{T}\right).
\end{aligned}
```

The default ambient relative humidity and temperatures produce
condensation for the supplied initial population.

### Nucleation boundary contribution

The example uses the paper model's prescribed constant birth flux, active only
for

$$
0.20\ \mathrm{s} \le t < 0.40\ \mathrm{s}.
$$

Let $M_{0,0,0}(0)$ be the reconstructed initial number density. The birth flux
is

$$
J = \frac{M_{0,0,0}(0)}{0.40\ \mathrm{s}-0.20\ \mathrm{s}},
$$

so nucleation adds one initial population over its active interval. $J$ has
units of number density per second. Every newborn droplet has the deterministic
phase-space state

$$
(r_n,w_n,T_n)=(50.0\times10^{-6}\ \mathrm{m},0\ \mathrm{m/s},298.15\ \mathrm{K}).
$$

For each configured raw moment $M_{a,b,c}$, the native
`CqmomMomentBoundaryContribution` returns

$$
\left.\frac{\mathrm{d}M_{a,b,c}}{\mathrm{d}t}\right|_{\mathrm{nuc}}
= J r_n^a w_n^b T_n^c.
$$

QMOMKit combines this rate with the zero volumetric moment source during the
moment-source stage and reconstructs the updated CQMOM state. The application
iterates the provided `CqmomSourceContext` entries, so it follows the moment
layout selected by the configured hierarchy rather than assuming a full
tensor.

### Momentum equation

If $w_{\mathrm{rel}}=0$, the drag acceleration is zero. Otherwise, for
$`\mathrm{Re}<\mathrm{Re}_{\min}`$, the low-Reynolds-number branch uses
```math
a_{\mathrm{drag}} = \frac{4.5\,\mu_a\,w_{\mathrm{rel}}}{\rho_l r^2}.
```
For $`\mathrm{Re}\ge \mathrm{Re}_{\min}`$, the deformation factor is
[[Lie et al., 1993](https://saemobilus.sae.org/papers/modeling-effects-drop-drag-breakup-fuel-sprays-930072);
[Hsiang & Faeth, 1992](https://doi.org/10.1016/0301-9322(92)90036-G)]
```math
F_{\mathrm{def}} = 1 + 2.632 \left[ 1 - \frac{1}{\left(1+0.19\sqrt{\mathrm{We}}\right)^2} \right],
```
the spherical drag coefficient correlation is
[[Clift & Gauvin, 1971](https://doi.org/10.1002/cjce.5450490403)]
```math
C_{d,0} = \frac{24}{\mathrm{Re}}
\left[1 +0.15\,\mathrm{Re}^{0.687} + \frac{0.0175\,\mathrm{Re}}{1+42500\,\mathrm{Re}^{-1.16}} \right],
```
and
$$C_d = C_{d,0}F_{\mathrm{def}}.$$

The general drag acceleration is
```math
\begin{aligned}
a_{\mathrm{drag}}
&=
\frac{
3\rho_a C_d |w_{\mathrm{rel}}| w_{\mathrm{rel}}
}{
8\rho_l r
}.
\end{aligned}
```

The vertical-velocity equation is

$$
\begin{aligned}
\frac{\mathrm{d}w}{\mathrm{d}t}
&=
-g\left(1-\frac{\rho_a}{\rho_l}\right)
+a_{\mathrm{drag}}.
\end{aligned}
$$

The gravity term is negative in the upward-positive sign convention, while the
drag term acts according to the sign of $w_{\mathrm{rel}}$.

### Temperature equation

The temperature rate combines sensible heat transfer and latent heat associated
with the radius change:

```math
\begin{aligned}
\frac{\mathrm{d}T}{\mathrm{d}t}
&=
\frac{3}{r c_l}
\left[
\frac{
\mathrm{Nu}\,k_a\,(T_\infty-T)
}{
2r\rho_l
}
+
L(T)\frac{\mathrm{d}r}{\mathrm{d}t}
\right].
\end{aligned}
```

The first term inside the brackets relaxes droplet temperature toward the
ambient temperature. The second term couples condensation or evaporation to the
droplet thermal state through the implemented latent heat correlation.

## Initialization shared by CQMOM and Monte Carlo

The supplied initialization is an independent product distribution with total
number density

$$
n_0 = 1.
$$

The configured marginals are:

- $r$ is gamma-distributed in metres with `shape = 16.0` and `scale = 18.75e-6`.
- $w$ is normal with `mean = -2.0` and `standard_deviation = 1.0`.
- $T$ is normal with `mean = 293.15` and `standard_deviation = 3.0`.

For the implemented gamma parameterization,

```math
\begin{aligned}
\mathbb{E}[r_{\mathrm{\mu m}}]
&=
k\theta_r
=
16.0 \times 18.75
=
300\ \mathrm{\mu m},
\end{aligned}
```

where $k$ is the gamma shape and $\theta_r$ is the gamma scale.

Both executables start from the same sample of $10000$ droplets.
The CQMOM executable computes empirical raw
moments from that deterministic finite sample and builds a sampled-moment
CQMOM case. The direct Monte Carlo executable evolves the same sampled
population.

## Numerical workflow and configuration

The numerical CQMOM hierarchy is selected in `configuration.qmom`:

```ini
[phase_space]
names = r, w, T
node_counts = 3, 3, 3
conditioning_order = 0, 2, 1

[qmom]
adaptive = true

[time_integration]
scheme = strang

[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-6
absolute_tolerance = 1.0e-9
initial_step = 1.0e-4
max_step = 1.0e-2
min_step = 1.0e-9
```

The `node_counts = 3, 3, 3` setting gives up to three nodes in $r$, three
conditional nodes in $w$, and three conditional nodes in $T$, for up to 27
terminal CQMOM nodes. The phase-space coordinate order remains $(r,w,T)$,
while `conditioning_order = 0, 2, 1` reconstructs $r$, then $T$, then $w$.
The CQMOM conditioning order is independent of the phase-space coordinate
order. With
`qmom.adaptive = true`, these configured node counts are maximum
branch capacities and the accepted runtime topology may contain fewer active
terminal nodes.

The configuration does not set `initial.maximum_orders`. When the sampling
helper does not receive explicit maximum orders, it uses the CQMOM hierarchy
moment extents, which are $2$ times the configured node count in each
dimension. For the default `3, 3, 3` hierarchy, that corresponds to sampled
orders $0,\ldots,5$ in $r$, $w$, and $T$. The result files are written
with the hierarchy-defined required CQMOM moment layout; there is no hardcoded
output-moment list in the example.

The driver time controls are compiled local constants:

```ini
initial_time = 0.0
final_time = 1.0
time_step = 0.01
```

The output settings are read from:

```ini
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
```

The CQMOM executable reconstructs the initial three-dimensional CQMOM
quadrature, creates a `qmomkit::CqmomState`, and uses
`qmomkit::advance_cqmom_phase_space(...)` while nucleation is active.
Configured Strang splitting applies advection for a half step, the boundary
birth contribution for a full source step, then advection for a half step. At
other times it uses the existing native advection workflow. Phase-space
advection moves quadrature abscissas while preserving terminal weights; the
boundary contribution changes raw moments and may change the accepted adaptive
topology after source reconstruction.

The Monte Carlo executable uses the same advection callback and adaptive-RK
tolerances. During the nucleation interval it advances the ensemble for a half
step, calls `qmomkit::monte_carlo::inject_population(...)` with $J\Delta t$,
then advances the enlarged ensemble for the second half step. The injection API
keeps a common representative sample weight and returns the increased total
number density.

## Build and run

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
cd examples/cqmom_3d_droplets_condensation-nucleation

cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

This builds:

- `cqmom_3d_droplets_condensation_nucleation_example`
- `cqmom_3d_droplets_condensation_nucleation_mc_example`

CMake copies `configuration.qmom` into the example build directory. Run from
that build directory so generated output stays there:

```bash
cd build
./cqmom_3d_droplets_condensation_nucleation_example
./cqmom_3d_droplets_condensation_nucleation_mc_example
```

Both executables print the same selection of initial and final raw moments, so
the two runs can be compared line by line: the number density $M_{0,0,0}$ and
the first and second moments of $r$, $w$ and $T$. The Monte Carlo executable
also prints the initial and final droplet counts.

## Output files

With the default output configuration, the CQMOM executable writes:

```text
output/moment_history.txt
output/quadrature_history.txt
```

`moment_history.txt` contains the CQMOM raw-moment history at the configured
output interval, starting with the initial state. It is a whitespace-delimited
`MomentWriter` table with a commented header containing `time` followed by
moment labels such as `M[0,0,0]`.

`quadrature_history.txt` contains terminal CQMOM nodes over time. It is a
`QuadratureWriter` table with time, public node index, conditioning path
entries, terminal weight, and abscissas ordered as $r$ in metres, $w$ in
m/s, and $T$ in K.

The direct Monte Carlo executable writes:

```text
output/moment_history_mc.txt
```

`moment_history_mc.txt` contains Monte Carlo estimated raw moments from the
evolving finite sample. Population injection preserves a common representative
sample weight while increasing both sample count and total number density. It
uses the same case-defined moment columns, labels, ordering, and output times
as the CQMOM moment history. The initial rows agree to roundoff because both
workflows start from the same deterministic finite sampled population.

## Plotting

The repository provides generic plotting helpers for QMOMKit-native text
output. From the repository root, after running the example, plot moment
histories with:

```bash
python tools/plotting/plot_moments.py examples/cqmom_3d_droplets_condensation-nucleation/build/output/moment_history.txt
python tools/plotting/plot_moments.py examples/cqmom_3d_droplets_condensation-nucleation/build/output/moment_history_mc.txt
```

Selected moment columns can be requested by label:

```bash
python tools/plotting/plot_moments.py examples/cqmom_3d_droplets_condensation-nucleation/build/output/moment_history.txt --moments 'M[0,0,0]' 'M[1,0,0]' 'M[0,1,0]' 'M[0,0,1]'
```

Plot the CQMOM quadrature history with:

```bash
python tools/plotting/plot_quadrature.py examples/cqmom_3d_droplets_condensation-nucleation/build/output/quadrature_history.txt
```

The plotting scripts read the existing whitespace-delimited output files and
write PNG files next to the input file.

## Expected qualitative behavior

For the default case, the vapor-drive term is positive for the supplied
ambient and initial temperatures, so seeded droplets grow by condensation.
From $0.20$ to $0.40$ s, nucleation increases total droplet number by one
initial population in CQMOM and by the corresponding Poisson realization in
Monte Carlo. The newborn droplets then also condense. Droplet temperatures
relax toward the prescribed ambient temperature while latent heat couples the
radius change back into the temperature equation. Vertical velocities evolve
under buoyancy-corrected gravity and drag, and droplets of different radii
respond differently because the Reynolds number, drag coefficient, heat
transfer, and mass transfer depend on $r$.

At $t = 1$ s, the CQMOM run gives $M_{0,0,0} = 2$ and the Monte Carlo run
$1.9947$ from 19,947 droplets. The mean temperature $M_{0,0,1}/M_{0,0,0}$ is
about 313 K in both runs, near the ambient 313.15 K, and the other printed
moments agree to within 0.3 percent.