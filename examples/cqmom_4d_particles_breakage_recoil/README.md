# 4D CQMOM particles: breakage and recoil

This standalone capability example evolves a spatially homogeneous particle
population in phase space $(d,u,v,w)$. It pairs a fixed CQMOM calculation with
a direct Monte Carlo (MC) companion. Both use deliberately simple constant-rate
binary breakage with recoil; the point is to demonstrate a nonlinear
parent-to-daughter source transformation evaluated from the source-stage
terminal quadrature.

The coordinate order is diameter $d$ in m followed by particle velocity
components $(u,v,w)$ in $\mathrm{m/s}$. Raw number-density moments are

$$
M_{p,q,r,s}=\int d^p u^q v^r w^s f(d,u,v,w)\,\mathrm{d}d\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w.
$$

In particular, $M_{0,0,0,0}$ is number density and $(\pi/6)M_{3,0,0,0}$ is
total particle volume. The configuration uses two CQMOM nodes in every
dimension, fixed branch counts, and conditioning order $(0,1,2,3)$.

## Initialization and continuous evolution

Each executable deterministically samples 10,000 particles from the same
independent-product distribution with number density one:

$$
d\sim\operatorname{LogNormal}(-9.210340371976184,0.12^2),\qquad
u,v,w\sim\mathcal N(0,0.5^2).
$$

The CQMOM calculation initializes from raw moments estimated from that
ensemble; the MC calculation advances the ensemble itself. Both verify the
sampled tensor before evolution, so they begin from the same finite-sample
population.

The carrier gas has $(\rho_g,\mu_g)=(1.2\ \mathrm{kg\,m^{-3}},1.8\times10^{-5}\
\mathrm{Pa\,s})$ and velocity $(2.0,1.5,0.8) \ \mathrm{m/s}$. Particles have
density $2000 \ \mathrm{kg/m}^{3}$. With slip
$\boldsymbol{s}=\boldsymbol{U}_g-\boldsymbol{u}$, the phase-space callback uses
the Schiller--Naumann form

$$
\mathrm{Re}=\frac{\rho_gd|\boldsymbol{s}|}{\mu_g},\qquad
\frac{\mathrm{d}\boldsymbol{u}}{\mathrm{d}t}=
\frac{3\rho_g}{4\rho_pd}\frac{24}{\mathrm{Re}}
\left(1+0.15\,\mathrm{Re}^{0.687}\right)|\boldsymbol{s}|\boldsymbol{s},
$$

with $\mathrm{d}d/\mathrm{d}t=0$.

## Breakage source and terminal quadrature

Each parent breaks at $\beta=10 \ \mathrm{s}^{-1}$ into two equal-volume
daughters,

$$
d'=2^{-1/3}d,
$$

with opposite recoil vectors. For an individual daughter, each recoil component
is normally distributed with zero mean and standard deviation

$$
\sigma_r(d)=0.25\frac{d}{10^{-4}\ \mathrm{m}}\ \mathrm{m\,s^{-1}}.
$$

The MC companion samples one recoil vector and assigns its positive and
negative forms to the two daughters. The CQMOM source uses the corresponding
conditional Gaussian moments and evaluates every requested birth-minus-loss
contribution through `CqmomSourceContext::for_each_terminal_node(...)`. This
visits the global terminal weights and $(d,u,v,w)$ coordinates of the
quadrature reconstructed from the current source stage.

Diameter-scaled recoil couples velocity powers to higher diameter powers. For
example, the source of the tracked $M_{0,3,0,0}$ contains the untracked mixed
moment $M_{2,1,0,0}$, because

$$
\mathbb E[(u+r_u)^3]=u^3+3u\,\sigma_r(d)^2.
$$

Applying the nonlinear daughter mapping directly to terminal nodes supplies
this reconstruction-dependent source closure naturally.

Breakage alone conserves total particle volume and volume-weighted momentum:
the equal-volume daughters preserve $d^3$, and their opposite recoils cancel in
the pair. It increases number density, while recoil increases velocity spread;
drag can change momentum during the coupled run.

## Build and run

This example consumes both QMOMKit and QMOMKitMonteCarlo, so install the two of
them to a common prefix first. Then build and run from this example directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_4d_particles_breakage_recoil_example
./cqmom_4d_particles_breakage_recoil_mc_example
```

The configured Strang-split CQMOM run and the MC companion both advance from
$0$ to $0.1$ s with an outer step of $0.001$ s. The MC path applies half drag
advection, stochastic breakage, then half drag advection. It uses fixed
initialization and breakage seeds; it remains a finite-sample reference rather
than an exact duplicate of the CQMOM result.

## Output files

The CQMOM executable writes `output/moment_history.txt` and `output/quadrature_history.txt` at the configured $0.01$ s cadence. The first contains time and the moments represented by the current CQMOM state; the second contains terminal-node weights and coordinates. The MC executable writes its matching selected raw moments to `output/moment_history_mc.txt`. Both executables print the initial and final tracked moments (up to 20 on the console).
