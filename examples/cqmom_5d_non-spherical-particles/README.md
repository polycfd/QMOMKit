# Five-dimensional non-spherical particle CQMOM example

## Purpose and workflows

This standalone example advances a five-dimensional non-spherical particle
population with CQMOM and a direct Monte Carlo comparison. The model includes
drag in an oscillating carrier flow, gravity, regularized diameter attrition,
advection aspect-ratio relaxation, binary fragmentation, and recoil velocity
spread.

It includes three executables:

- `cqmom_5d_non_spherical_particles_example`: the normal bounded-support CQMOM
  simulation.
- `cqmom_5d_non_spherical_particles_unbounded_diagnostic`: an intentionally
  unbounded diagnostic that explains a conditional support-loss event.
- `cqmom_5d_non_spherical_particles_mc_example`: a direct stochastic Monte Carlo
  comparison for the same particle model and sampled initialization.

The bounded CQMOM executable is the normal physical workflow. The unbounded
diagnostic uses the same initial population and physical model, but omits the
aspect-ratio support declaration to show why finite conditional moments can be
Stieltjes-realizable while failing bounded Hausdorff realizability on $[0,1]$.
The Monte Carlo executable uses the same sampled initialization and particle
model to provide an event-driven stochastic reference workflow.

Each CQMOM workflow loads its configured case, then uses
`with_initial_moments(...)` to replace only the initialized tensor with the
complete sampled moments. The configured hierarchy, support policy, storage
selection, and reconstruction policies remain unchanged; the sampling utility
is only one possible tensor producer.

## Phase-space definition and raw moments

The particle phase-space coordinate is

$$
\boldsymbol{\xi} = (d, \chi, u, v, w),
$$

where $d$ is particle diameter in meters, $\chi$ is dimensionless particle
aspect ratio, and $u$, $v$, and $w$ are particle-velocity components in meters
per second.

Phase-space advection uses the phase-space velocity

$$
\frac{\mathrm{d}\boldsymbol{\xi}}{\mathrm{d}t}
=
\left(
\frac{\mathrm{d}d}{\mathrm{d}t},
\frac{\mathrm{d}\chi}{\mathrm{d}t},
\frac{\mathrm{d}u}{\mathrm{d}t},
\frac{\mathrm{d}v}{\mathrm{d}t},
\frac{\mathrm{d}w}{\mathrm{d}t}
\right).
$$

Binary fragmentation is treated separately as a moment-source population
source acting on raw moments or, in the Monte Carlo executable, on sampled
particles.

For a multi-index
$\boldsymbol{\alpha}=(\alpha_d,\alpha_\chi,\alpha_u,\alpha_v,\alpha_w)$, the
raw moment is

$$
M_{\boldsymbol{\alpha}}
=
\int
d^{\alpha_d}
\chi^{\alpha_\chi}
u^{\alpha_u}
v^{\alpha_v}
w^{\alpha_w}
f(d,\chi,u,v,w)
\,\mathrm{d}d\,\mathrm{d}\chi\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w.
$$

These moments are unnormalized number-density moments. In particular,
$M_{(0,0,0,0,0)}$ is total number density. Mean quantities are therefore
obtained by dividing first moments by $M_{(0,0,0,0,0)}$; the implementation
does not normalize raw moments internally.

## Complete particle model

The physical constants below are compiled into this example. They are
example-local model values, not universal calibrated material constants.

| Symbol | Code member | Value | Units | Role |
|---|---|---:|---|---|
| $\rho_f$ | `fluid_density` | `1.2` | kg/m^3 | carrier-fluid density |
| $\rho_p$ | `particle_density` | `2000.0` | kg/m^3 | particle material density |
| $\mu_f$ | `fluid_dynamic_viscosity` | `1.8e-5` | Pa s | carrier-fluid dynamic viscosity |
| $g$ | `gravity` | `9.81` | m/s^2 | downward gravitational acceleration |
| $A_x$ | `fluid_velocity_amplitude_x` | `2.0` | m/s | x carrier-velocity amplitude |
| $A_y$ | `fluid_velocity_amplitude_y` | `2.0` | m/s | y carrier-velocity amplitude |
| $A_z$ | `fluid_velocity_amplitude_z` | `0.0` | m/s | z carrier-velocity amplitude |
| $\omega$ | `oscillation_angular_frequency` | `1.0` | 1/s | carrier-flow angular frequency |
| $k_{\mathrm{attr}}$ | `attrition_rate` | `1.0e-6` | s/m | diameter attrition coefficient |
| $U_{\mathrm{sat}}$ | `attrition_velocity_scale` | `5.0` | m/s | velocity scale for attrition saturation |
| $d_{\min}$ | `minimum_attrition_diameter` | `0.0` | m | lower scale in attrition taper |
| $d_{\mathrm{transition}}$ | `attrition_transition_diameter` | `1.0e-5` | m | attrition taper transition length |
| $\tau_\chi$ | `aspect_ratio_relaxation_time` | `2.0` | s | advection relaxation time toward $\chi=1$ |
| $k_{\mathrm{frag}}$ | `fragmentation_rate_constant` | `1.0e8` | m^-2 s^-1 | coefficient in $\beta(d)$ for $p_{\mathrm{frag}}=2$ |
| $p_{\mathrm{frag}}$ | `fragmentation_power` | `2.0` | dimensionless | diameter exponent in fragmentation rate |
| $d_{\mathrm{frag,min}}$ | `minimum_fragmentation_diameter` | `1.0e-6` | m | fragmentation cutoff diameter |
| $f_\chi$ | `shape_relaxation_factor` | `0.5` | dimensionless | daughter aspect-ratio multiplier |
| $\sigma_{\mathrm{recoil}}$ | `recoil_velocity_standard_deviation` | `0.5` | m/s | standard deviation of each recoil component |

### Oscillating carrier-fluid velocity

The implemented carrier-fluid velocity is spatially uniform and time
dependent:

$$
\boldsymbol{U}_f(t)
=
\left(
A_x\cos(\omega t),
A_y\sin(\omega t),
A_z\cos(\omega t)
\right).
$$

Equivalently,

$$
\begin{aligned}
U_{f,x}(t) &= A_x\cos(\omega t),\\
U_{f,y}(t) &= A_y\sin(\omega t),\\
U_{f,z}(t) &= A_z\cos(\omega t).
\end{aligned}
$$

### Relative velocity and Reynolds number

The particle velocity is $\boldsymbol{U}_p=(u,v,w)$. The relative velocity
driving both drag and attrition is

$$
\boldsymbol{U}_{\mathrm{rel}}
=
\boldsymbol{U}_f-\boldsymbol{U}_p,
$$

with components

$$
\begin{aligned}
U_{\mathrm{rel},x} &= U_{f,x}-u,\\
U_{\mathrm{rel},y} &= U_{f,y}-v,\\
U_{\mathrm{rel},z} &= U_{f,z}-w.
\end{aligned}
$$

Its magnitude is

$$
\left|\boldsymbol{U}_{\mathrm{rel}}\right|
=
\sqrt{
(U_{f,x}-u)^2
+(U_{f,y}-v)^2
+(U_{f,z}-w)^2
}.
$$

The particle Reynolds number is

$$
\mathrm{Re}
=
\frac{
\rho_f\left|\boldsymbol{U}_{\mathrm{rel}}\right|d
}{
\mu_f
}.
$$

Here $\rho_f$ has units kg/m^3,
$\left|\boldsymbol{U}_{\mathrm{rel}}\right|$ has units m/s, $d$ has units m,
and $\mu_f$ has units kg/(m s), so $\mathrm{Re}$ is dimensionless.

### Non-spherical drag correction

The drag response uses an aspect-ratio-dependent correction for an oblate
spheroidal particle. For physical aspect ratios in $(0,1]$, define

$$
e=\sqrt{1-\chi^2}.
$$

If $\left|1-\chi\right|<0.001$, the implementation uses the spherical-limit
value

$$
K_0(\chi)=1.
$$

Otherwise it evaluates

$$
K_0(\chi)
=
\frac{
\frac{8}{3\chi^{1/3}}
}{
\frac{2\chi}{1-\chi^2}
+
\frac{
2(1-2\chi^2)\arctan(e/\chi)
}{
(1-\chi^2)^{3/2}
}
}.
$$

The complete shape factor includes empirical Reynolds-number corrections:

$$
K(\chi,\mathrm{Re})
=
K_0(\chi)
+0.15\,\chi^{95.91}\mathrm{Re}^{0.687}
+0.2927\,(1-\chi)^{0.4374}\mathrm{Re}^{0.7512}.
$$

The first term supplies the oblate-spheroid creeping-flow correction. The
second and third terms increase the drag response with Reynolds number and make
the correction depend on particle shape. The formula is evaluated directly by
the example callback; the callback does not clip $d$ or $\chi$.

### Drag acceleration and gravity

The implemented drag-rate coefficient is

$$
\lambda_{\mathrm{drag}}
=
\frac{
18\mu_f K(\chi,\mathrm{Re})
}{
\rho_p d^2
}.
$$

It has units 1/s. With projected spherical area $A=πd^2/4$, this is
algebraically equivalent to evaluating the drag force

$$
\boldsymbol{F}_d
=
\frac{1}{2}
\rho_f A
\left(
\frac{24K(\chi,\mathrm{Re})}{\mathrm{Re}}
\right)
\left|\boldsymbol{U}_{\mathrm{rel}}\right|
\boldsymbol{U}_{\mathrm{rel}}
$$

and dividing by the spherical particle mass $\rho_p\pi d^3/6$. Substituting
the Reynolds number gives
$\boldsymbol{F}_d/(\rho_p\pi d^3/6)=\lambda_{\mathrm{drag}}\boldsymbol{U}_{\mathrm{rel}}$.
The Reynolds-number dependence in $K(\chi,\mathrm{Re})$ remains part of the
model; no additional Reynolds-number factor multiplies the force. The velocity
equations are

$$
\begin{aligned}
\frac{\mathrm{d}u}{\mathrm{d}t}
&=
\lambda_{\mathrm{drag}}(U_{f,x}-u),\\
\frac{\mathrm{d}v}{\mathrm{d}t}
&=
\lambda_{\mathrm{drag}}(U_{f,y}-v),\\
\frac{\mathrm{d}w}{\mathrm{d}t}
&=
\lambda_{\mathrm{drag}}(U_{f,z}-w)-g.
\end{aligned}
$$

Drag relaxes each particle-velocity component toward the carrier-fluid
velocity. Because $\lambda_{\mathrm{drag}}$ scales like $1/d^2$, smaller
particles respond more rapidly to drag. Gravity acts only in the $w$ equation,
with positive $g$ subtracted from $\mathrm{d}w/\mathrm{d}t$; in this sign
convention gravity accelerates particles in the negative $w$ direction.

### Regularized diameter attrition

Diameter attrition uses the same relative speed as drag, but with two
regularizations. First, the velocity-squared response is saturated:

$$
U_{\mathrm{attr}}^2
=
\frac{
\left|\boldsymbol{U}_{\mathrm{rel}}\right|^2
}{
1+\left|\boldsymbol{U}_{\mathrm{rel}}\right|^2/U_{\mathrm{sat}}^2
}.
$$

Second, the recession rate is tapered near the configured lower diameter scale:

$$
\frac{\mathrm{d}d}{\mathrm{d}t}
=
-k_{\mathrm{attr}}U_{\mathrm{attr}}^2
\tanh\left(
\frac{d-d_{\min}}{d_{\mathrm{transition}}}
\right).
$$

For ordinary relative velocities, $U_{\mathrm{attr}}^2$ approaches the original
quadratic dependence
$\left|\boldsymbol{U}_{\mathrm{rel}}\right|^2$. At extreme relative velocity,
$U_{\mathrm{attr}}^2$ is bounded by the saturation scale. Far above $d_{\min}$,
the $\tanh$ factor approaches one and the original recession behavior is
recovered. Near $d_{\min}$, the taper makes the boundary asymptotically
attracting and avoids finite-time crossing.

The current minimum attrition diameter is $d_{\min}=0$. This regularization
changes the attrition law itself; it is not state clipping, and it does not
remove extreme velocity quadrature nodes.

### Advection aspect-ratio relaxation

The advection aspect-ratio equation is

$$
\frac{\mathrm{d}\chi}{\mathrm{d}t}
=
\frac{1-\chi}{\tau_\chi}.
$$

The relaxation time is $\tau_\chi=2.0$ s. This advection term drives $\chi$
toward the spherical value $\chi=1$. It is distinct from the moment-source
daughter aspect-ratio transformation used during fragmentation.

## Binary fragmentation

### Fragmentation rate

Binary fragmentation uses the diameter-dependent event rate

$$
\beta(d)
=
\begin{cases}
k_{\mathrm{frag}}d^{p_{\mathrm{frag}}},
& d\ge d_{\mathrm{frag,min}},\\
0,
& d<d_{\mathrm{frag,min}}.
\end{cases}
$$

with $k_{\mathrm{frag}}=1.0\times10^8$ m^-2 s^-1,
$p_{\mathrm{frag}}=2.0$, and $d_{\mathrm{frag,min}}=1.0\times10^{-6}$ m. The
minimum fragmentation diameter stops further fragmentation events for smaller
particles, but it is not a lower bound on particle diameter under advection
attrition.

### Daughter diameter and volume conservation

One parent produces two daughters. Daughter diameter is

$$
d_{\mathrm{daughter}}
=
2^{-1/3}d_{\mathrm{parent}},
$$

so

$$
2d_{\mathrm{daughter}}^3
=
d_{\mathrm{parent}}^3.
$$

Binary fragmentation therefore conserves diameter-cubed content for each
event.

### Daughter aspect ratio

The daughter aspect ratio is

$$
\chi_{\mathrm{daughter}}
=
f_\chi\chi_{\mathrm{parent}},
$$

with $f_\chi=0.5$. This moment-source multiplicative shape change is separate
from the advection relaxation toward $\chi=1$. If the parent is inside the
physical interval $[0,1]$, this daughter rule also stays inside $[0,1]$.

### Gaussian recoil velocity

Each Monte Carlo fragmentation event samples one centered Gaussian recoil
vector

$$
\boldsymbol{\zeta}
=
(\zeta_u,\zeta_v,\zeta_w),
\qquad
\zeta_j\sim\mathcal{N}(0,\sigma_{\mathrm{recoil}}^2),
$$

with $\sigma_{\mathrm{recoil}}=0.5$ m/s. The paired daughter velocities are

$$
\begin{aligned}
\boldsymbol{U}_{\mathrm{daughter},1}
&=
\boldsymbol{U}_{\mathrm{parent}}+\boldsymbol{\zeta},\\
\boldsymbol{U}_{\mathrm{daughter},2}
&=
\boldsymbol{U}_{\mathrm{parent}}-\boldsymbol{\zeta}.
\end{aligned}
$$

Equal and opposite kicks give zero net pair recoil. Because the daughters have
equal volume and the two daughter volumes sum to the parent volume, this paired
recoil preserves mass-weighted momentum for each fragmentation event. It also
increases velocity spread and velocity second moments.

The deterministic CQMOM source does not sample recoil. Instead, daughter
velocity moments use analytical shifted Gaussian moments

$$
G_k(x,\sigma)
=
\mathbb{E}\left[(x+\zeta)^k\right],
\qquad
\zeta\sim\mathcal{N}(0,\sigma^2).
$$

Only even centered-Gaussian powers contribute to the binomial expansion:

$$
G_k(x,\sigma)
=
\sum_{\substack{r=0\\r\ \mathrm{even}}}^{k}
\binom{k}{r}
x^{k-r}
\sigma^r
(r-1)!!.
$$

Representative low-order moments are

$$
\begin{aligned}
G_0 &= 1,\\
G_1 &= x,\\
G_2 &= x^2+\sigma^2,\\
G_3 &= x^3+3x\sigma^2.
\end{aligned}
$$

### Deterministic CQMOM moment source

For a raw-moment multi-index
$\boldsymbol{\alpha}=(\alpha_d,\alpha_\chi,\alpha_u,\alpha_v,\alpha_w)$, the
parent monomial is

$$
P_{\boldsymbol{\alpha}}
=
d^{\alpha_d}
\chi^{\alpha_\chi}
u^{\alpha_u}
v^{\alpha_v}
w^{\alpha_w}.
$$

The expected one-daughter monomial is

$$
D_{\boldsymbol{\alpha}}
=
d_{\mathrm{daughter}}^{\alpha_d}
\chi_{\mathrm{daughter}}^{\alpha_\chi}
G_{\alpha_u}(u,\sigma_{\mathrm{recoil}})
G_{\alpha_v}(v,\sigma_{\mathrm{recoil}})
G_{\alpha_w}(w,\sigma_{\mathrm{recoil}}).
$$

The deterministic CQMOM source sums terminal-quadrature contributions as

$$
S_{\boldsymbol{\alpha}}
=
\sum_i
W_i
\beta(d_i)
\left(
2D_{\boldsymbol{\alpha},i}
-P_{\boldsymbol{\alpha},i}
\right).
$$

The $-P_{\boldsymbol{\alpha}}$ term removes one parent,
$2D_{\boldsymbol{\alpha}}$ creates two statistically identical daughters,
$W_i$ supplies the parent population weight, and the daughter velocity factors
include the Gaussian recoil expectation.

### Monte Carlo fragmentation process

The Monte Carlo source applies the same event model stochastically. For a
particle with rate $\beta(d)$, the waiting time is exponentially distributed:

$$
T\sim\mathrm{Exp}(\beta(d)).
$$

If $T$ falls inside the remaining moment-source interval, the parent is
removed and two daughters are inserted with
$d_{\mathrm{daughter}}=2^{-1/3}d_{\mathrm{parent}}$,
$\chi_{\mathrm{daughter}}=f_\chi\chi_{\mathrm{parent}}$, and paired recoil
velocities. Daughter lineages re-enter the same moment-source stage with the remaining
stage time, so repeated fragmentation in one source interval is allowed.
Survivors and newly generated daughters keep the same macro-particle weight.

### Fragmentation-stage balances

For the moment-source fragmentation operator alone, particle number increases
by one per binary event:

$$
S_{(0,0,0,0,0)}
=
\sum_i W_i\beta(d_i).
$$

Diameter-cubed content is conserved by the daughter diameter rule:

$$
S_{(3,0,0,0,0)}=0.
$$

Paired recoil preserves diameter-cubed-weighted velocity moments eventwise for
equal-volume daughters:

$$
S_{(3,0,1,0,0)}
=
S_{(3,0,0,1,0)}
=
S_{(3,0,0,0,1)}
=
0.
$$

Recoil increases velocity spread and velocity second moments. These are
fragmentation-only balances. The complete model also includes advection drag,
gravity, diameter attrition, and aspect-ratio relaxation. Advection drag
exchanges momentum with the carrier flow, gravity accelerates the vertical
velocity component, and attrition changes particle volume and diameter moments.

## Initialization shared by CQMOM and Monte Carlo

The supplied initialization is an independent-product distribution with total
number density $n_0=1$. The coordinate distributions are configured in
`configuration.qmom`:

- $d$ is lognormal with `mu_log = -8.9596725473` and
  `sigma_log = 0.1529473098`.
- $\chi$ is beta-distributed with `alpha = 6.0` and `beta = 14.0`.
- $u$ uses the native CQMOM marginal-mixture syntax: two normal components,
  each with `number_density = 0.5` and `standard_deviation = 0.5`, centered at
  `mean = -1.0` and `mean = 1.0`, respectively.
- $v$ and $w$ are normal with `mean = 0.0` and
  `standard_deviation = 0.5`.

The $u$ marginal is therefore a bimodal distribution of two simultaneously
present velocity subpopulations. It provides a simple phase-space
representation of crossing velocity populations (and hence particle
trajectory crossing in a population sense); it does not assign persistent
individual identities or trajectories to CQMOM nodes. The configured
independent-product initialization combines this marginal mixture with the
other four marginals.

The CQMOM drivers parse the independent-product initialization, sample `10000`
particles with `initialization_seed = 12345`, and replace the distributional
initial condition with raw moments of that same finite ensemble. This makes the
bounded and unbounded CQMOM initial states identical except for the support
declaration.

The Monte Carlo executable uses the same sampled ensemble directly. Each
sample carries $n_0/N$, so the ensemble represents unnormalized raw moments of
the particle number-density distribution.

## Numerical workflow and configurable choices

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

Both CQMOM configurations use Strang splitting for phase-space evolution and
RK23 for node trajectories:

```text
[time_integration]
scheme = strang

[time_integration.node_advection]
scheme = rk23
```

with RK23 tolerances and step limits in:

```text
[time_integration.node_advection.adaptive]
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

The Monte Carlo driver uses `time_step = 0.01` and `final_time = 1.0`; the
CQMOM driver uses the smaller `time_step = 0.005`. Both apply the same Strang
operator order: advection half-step, fragmentation source stage, then
advection half-step. The CQMOM `rk23` setting controls only quadrature-node
trajectories within its advection stages. The moment-source fragmentation is
advanced through QMOMKit's CQMOM update path and through event-driven waiting
times for Monte Carlo, so the workflows remain a deterministic closure and a
stochastic reference rather than numerically identical integrations.

## Bounded-support CQMOM workflow

`configuration.qmom` is the normal bounded configuration. It declares the
physical aspect-ratio support in phase-space order:

```text
support_lower = none, 0.0, none, none, none
support_upper = none, 1.0, none, none, none
```

QMOMKit uses this interval as a support-aware recursive branch acceptance
check. The bounded workflow advances the full physical model, completes the
full CQMOM phase-space evolution run, writes standard moment and quadrature histories, and
prints only concise support topology events.

## Intentionally unbounded diagnostic workflow

`configuration_unbounded.qmom` preserves the same hierarchy, initialization,
physical callbacks, Strang splitting, RK23 settings, and source model, but
omits the $\chi$ support declaration. It also selects the diagnostic output
directory `output/unbounded_diagnostic` and a diagnostic-specific
`max_steps` guard; neither changes the initialized population or physical
model.

The diagnostic advances the same physics to the first source boundary where
ordinary unconstrained post-source inversion accepts a finite terminal node
with $\chi<0$ or $\chi>1$. It then writes one self-contained Markdown report and
returns success before the next advection stage can obscure the origin of the
event.

## Why both CQMOM workflows are included

Together, the bounded and unbounded workflows demonstrate:

- the practical five-dimensional non-spherical-particle model;
- conditional support loss in finite multidimensional moment closures;
- the distinction between Stieltjes realizability and bounded Hausdorff
  realizability;
- why positive classical Hankel-Hadamard determinants do not guarantee
  $\chi\le1$;
- how support-aware adaptive reconstruction resolves the issue without
  clipping node coordinates.

The unbounded workflow is not a different physical model or initialization. It
is included to show how a finite conditional moment sequence can leave the
physical support even though the fragmentation daughter law preserves
$\chi\le1$ when the parent is inside $[0,1]$. The bounded workflow demonstrates
the generic support-aware response.

## Conditional support-loss mechanism

After the moment-source update, a branch-local conditional $\chi$
sequence can remain realizable on the nonnegative real line while failing
feasibility on the bounded physical interval $[0,1]$. This is a higher-order
conditional moment representation/support issue. It is not random sampling
noise, a negative-weight quadrature, a Vandermonde failure, a Wheeler defect,
or a physical daughter rule that creates $\chi>1$.

In the unbounded workflow, ordinary Wheeler reconstruction can accept all
requested branch moments and return positive weights while placing one
terminal $\chi$ node above one. In the bounded workflow, the declared interval
turns the same candidate into a support-infeasible branch and lets adaptive
reconstruction accept a lower-count support-valid representation.

## Classical Hankel-Hadamard and support-localizing conditions

For a branch sequence $m_0,\ldots,m_5$, the diagnostic reports classical
Hankel-Hadamard matrices and classical Hankel-Hadamard determinants

$$
H_{k,\ell}
=
\left[m_{k+i+j}\right]_{i,j=0}^{\ell},
\qquad
\Delta_{k,\ell}
=
\det(H_{k,\ell}),
$$

with $k=0$ and $k=1$. These checks establish real-line and
nonnegative-support moment conditions.

For interval $[a,b]$, bounded support also requires support-localizing
matrices

$$
L_{a,\ell}
=
\left[m_{i+j+1}-a\,m_{i+j}\right]_{i,j=0}^{\ell},
\qquad
U_{b,\ell}
=
\left[b\,m_{i+j}-m_{i+j+1}\right]_{i,j=0}^{\ell}.
$$

For $[0,1]$, the upper-support localizing matrix is

$$
U_{1,\ell}
=
\left[m_{i+j}-m_{i+j+1}\right]_{i,j=0}^{\ell}.
$$

The diagnostic classifies these matrices by eigenvalues relative to an
effective tolerance,

$$
\tau_{\mathrm{eff}}
=
10^{-12}
\max\left(1,\left\|A\right\|_{\infty}\right),
$$

and uses determinant values as explanatory diagnostics. The unbounded
diagnostic shows that the classical Hankel-Hadamard determinants can remain
positive while the upper-support localizing matrix fails. Wheeler can still be
algebraically accurate in this situation: it reconstructs the supplied
unbounded-support sequence, but one returned $\chi$ node lies outside the
finite physical interval.

## Support-aware adaptive reconstruction

The bounded workflow rejects the support-infeasible six-moment three-node
$\chi$ candidate, accepts a support-valid lower-count branch, rebuilds
descendants, and canonicalizes the accepted state from the accepted quadrature.

This is not clipping. It is not projection of all six moments onto a bounded
moment cone. It is an adaptive reduction of represented branch order:
higher-order source-candidate information is replaced by the closure implied by
the accepted lower-count quadrature.

The unbounded workflow uses all six branch moments. Ordinary Wheeler
reconstruction returns three positive $\chi$ nodes, one above one, and
reconstructs all six supplied moments to roundoff. The bounded workflow sees
the same source-boundary event with the same sampled realization and physical
model, rejects the infeasible three-node candidate because $\chi$ support is
declared, and retains a lower-count support-valid branch. The diagnostic report
compares the unbounded accepted branch with the actual bounded accepted branch
from the bounded post-source quadrature.

## Executables

| Executable | Purpose |
|---|---|
| `cqmom_5d_non_spherical_particles_example` | Normal bounded-support CQMOM simulation. |
| `cqmom_5d_non_spherical_particles_unbounded_diagnostic` | Early-stopping diagnostic for the first unconstrained $\chi$ support-loss event. |
| `cqmom_5d_non_spherical_particles_mc_example` | Direct stochastic Monte Carlo comparison. |

## Configuration files

`configuration.qmom` is the normal bounded-support configuration. It declares
the physical $\chi$ support interval and writes CQMOM output under `output`.

`configuration_unbounded.qmom` intentionally removes only the $\chi$ support
declaration, redirects diagnostic output to `output/unbounded_diagnostic`, and
sets a diagnostic-specific `max_steps` guard. It keeps the same phase-space
hierarchy, initialization, Strang splitting, RK23 tolerances, and physical
source callbacks.

Both configuration files are copied into the example build directory by the
standalone example CMake project.

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

cmake -S examples/cqmom_5d_non-spherical-particles \
    -B examples/cqmom_5d_non-spherical-particles/build \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build examples/cqmom_5d_non-spherical-particles/build
```

Run from the example build directory:

```bash
cd examples/cqmom_5d_non-spherical-particles/build
./cqmom_5d_non_spherical_particles_example
./cqmom_5d_non_spherical_particles_unbounded_diagnostic
./cqmom_5d_non_spherical_particles_mc_example
```

## Output files

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

## Expected qualitative behavior

The bounded CQMOM run remains within the declared $\chi$ support and reaches
`final_time = 1.0`. It may report concise support-driven active-count
reductions and support-infeasible growth rejections.

The bounded five-dimensional CQMOM run can take noticeably longer than the
introductory examples because it repeatedly reconstructs an adaptive 5D
hierarchy.

The unbounded diagnostic is expected to stop early after writing its report. In
the supplied configuration it detects a small positive-weight $\chi$ node above
one immediately after an accepted post-source reconstruction.

The Monte Carlo executable is a direct stochastic comparison and is
substantially more expensive than the CQMOM run: it starts with `10000`
particles, fragmentation grows the ensemble, every particle follows advection
dynamics, and dense moments are recomputed from the evolving ensemble. A long
Monte Carlo runtime is expected behavior, not a defect.

## Interpretation and limitations

This example is a deterministic numerical workflow and diagnostic demonstrator,
not a calibrated physical application. Coordinate names are example metadata;
QMOMKit itself does not assign physical units or enforce application physics.
The bounded support declaration is the correct way for this model to
communicate the physical aspect-ratio interval to generic CQMOM reconstruction.

The advection model does not include spatial phase-space evolution, particle-particle
collision physics beyond the configured binary fragmentation source, calibrated
material parameters, or application boundary conditions. The unbounded
diagnostic is intentionally nonphysical with respect to $\chi$ support; it is
included to expose a moment-realizability issue and to document the bounded
workflow's support-aware response.
