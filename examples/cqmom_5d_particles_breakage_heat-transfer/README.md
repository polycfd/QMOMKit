# 5D CQMOM particles: breakage and heat transfer

This standalone capability example evolves a spatially homogeneous particle
population in phase space $(V,T,u,v,w)$. It pairs fixed CQMOM with a direct
Monte Carlo (MC) companion and demonstrates the complementary source-closure
route to the 4D particle breakage and recoil example: the breakage source is
formulated in moment space and closed with integer and generalized moments,
alongside heat-transfer advection.

The coordinate order is volume $V$ in $\mathrm{m}^3$, temperature $T$ in K, and
velocity components $(u,v,w)$ in $\mathrm{m/s}$. The raw number-density moments
are

```math
M_{p,q,r,s,t}=\int V^pT^qu^rv^sw^t f(V,T,u,v,w)\,\mathrm{d}V\,\mathrm{d}T\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w.
```

Thus $M_{0,0,0,0,0}$ is number density and $M_{1,0,0,0,0}$ is total particle
volume. The configuration uses two fixed CQMOM nodes in every dimension
and conditioning order $(0,1,2,3,4)$.

## Initialization and continuous evolution

Each executable deterministically samples 10,000 particles from the same
independent-product distribution with number density one:

$$
V\sim\operatorname{LogNormal}(-28.278050699307204,0.36^2),\quad
T\sim\mathcal N(350,50^2)\ \mathrm{K},\quad
u,v,w\sim\mathcal N(0,0.5^2).
$$

The temperature distribution is centered at 350 K with a 50 K standard
deviation. The volume distribution is exactly the spherical transformation of
the 4D example's diameter distribution: for $V=\pi d^3/6$, its lognormal
parameters are $\ln(\pi/6)+3(-9.210340371976184)=-28.278050699307204$ and
$3(0.12)=0.36$.

The CQMOM path initializes from moments estimated from this ensemble; the MC
path advances the ensemble itself, so both begin from the same finite-sample
population.
The physical diameter used by the transport models is $d(V)=(6V/\pi)^{1/3}$.

The carrier gas has
$`(\rho_g,\mu_g)=(1.2\ \mathrm{kg\,m^{-3}},1.8\times10^{-5}\ \mathrm{Pa\,s})`$,
velocity $`(2.0,1.5,0.8)\ \mathrm{m\,s^{-1}}`$, thermal conductivity
$`0.026\ \mathrm{W\,m^{-1}\,K^{-1}}`$, Prandtl number $0.71$, and temperature
$T_g=300$ K. Particles have density $`2000\ \mathrm{kg\,m^{-3}}`$ and heat
capacity $`1000\ \mathrm{J\,kg^{-1}\,K^{-1}}`$. With
$\boldsymbol{s}=\boldsymbol{U}_g-\boldsymbol{u}$, the phase-space callback uses

```math
\mathrm{Re}=\frac{\rho_gd|\boldsymbol{s}|}{\mu_g},\qquad
\frac{\mathrm{d}\boldsymbol{u}}{\mathrm{d}t}=
\frac{3\rho_g}{4\rho_pd}\frac{24}{\mathrm{Re}}
\left(1+0.15\,\mathrm{Re}^{0.687}\right)|\boldsymbol{s}|\boldsymbol{s},
```

and Ranz--Marshall heat transfer,

```math
\mathrm{Nu}=2+0.6\,\mathrm{Re}^{1/2}\mathrm{Pr}^{1/3},\qquad
\frac{\mathrm{d}T}{\mathrm{d}t}=
\frac{6k_g\mathrm{Nu}}{\rho_pc_pd^2}(T_g-T).
```

Temperature is expressed directly in K in this equation and in all moment
indices; $\mathrm{d}V/\mathrm{d}t=0$.

## Breakage source and generalized moments

Each parent breaks at $\beta=10\ \mathrm{s^{-1}}$ into two identical,
equal-volume daughters with inherited velocity:

$$
V_d=\frac{V_p}{2},\qquad
T_d=T_p+\varepsilon\left(\frac{V_d}{V_{\mathrm{ref}}}\right)^{1/3}(T_g-T_p),\qquad
\varepsilon=0.01,
$$

where $V_{\mathrm{ref}}=\pi(10^{-4}\ \mathrm{m})^3/6$. The $V^{1/3}$ factor is
the deliberately simple surface-related temperature redistribution used by the
example. It also creates fractional volume powers when $T_d^q$ is expanded.

Integer volume orders are accessed with `CqmomSourceContext::moment(...)`;
fractional orders are evaluated from the source-stage quadrature by
`CqmomSourceContext::generalized_moment(...)`, without adding those generalized
moments to the transported state. For a concrete tracked moment, let
$a=\varepsilon(1/(2V_{\mathrm{ref}}))^{1/3}$. The current source for
$M_{1,1,0,0,0}$ is

$$
\dot M_{1,1,0,0,0}\big|_{\mathrm{breakage}}
=\beta a\left(T_gM_{4/3,0,0,0,0}-M_{4/3,1,0,0,0}\right).
$$

It therefore needs $M_{4/3,0,0,0,0}$ and $M_{4/3,1,0,0,0}$ to evaluate and
close the source of a tracked moment; it does not need them to reconstruct that
tracked moment. Their units are respectively number density times
$\mathrm{m^4}$ and number density times $`\mathrm{m^4\,K}`$, so the factor $a$
supplies the inverse-length scale needed by the source. For $q=0$, the source
reduces to $\beta(2^{1-p}-1)M_{p,0,r,s,t}$, preserving total volume.

## Build and run

This example consumes both QMOMKit and QMOMKitMonteCarlo, so install the two of
them to a common prefix first. Then build and run from this example directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_5d_particles_breakage_heat_transfer_example
./cqmom_5d_particles_breakage_heat_transfer_mc_example
```

The configured Strang-split CQMOM run and the MC companion advance from $0$ to
$0.1$ s with an outer step of $0.001$ s. The MC path applies half-step
heat-transfer advection, breakage, then another half step, with fixed
initialization and breakage seeds.

This advanced five-dimensional example can take longer than the introductory
examples. The direct Monte Carlo companion is substantially more expensive
than the CQMOM run because it evolves a finite population through breakage.

## Output files

The CQMOM executable writes `output/moment_history.txt` and
`output/quadrature_history.txt` at the first step that reaches each multiple of
the configured interval of $0.0025$ s. The first
contains time and the moments represented by the current CQMOM state; the
second contains terminal-node weights and coordinates. The MC executable writes
matching selected raw moments to `output/moment_history_mc.txt`. Both
executables print the same selection of initial and final raw moments, so the
two runs can be compared line by line: the number density, the first moments,
the second moment of each coordinate, and $M_{1,1,0,0,0}$, which size-dependent
heating couples. The MC executable also prints the initial and final particle
counts.

## What to expect

The equal-volume breakup transformation preserves total particle volume during
the breakage source while increasing particle number. Heat-transfer advection
changes the temperature moments, and the source evaluation exercises the
fractional-volume generalized moments required by the daughter-temperature
mapping. The Monte Carlo companion is a stochastic particle-based
companion/reference realization, so it is useful for comparing trends rather
than exact pointwise agreement.

At $t = 0.1$ s, the CQMOM run gives $M_{0,0,0,0,0} = 2.718$, which is $e$
because breakage at rate $10\ \mathrm{s^{-1}}$ multiplies the number density
by $e^{10 \cdot 0.1}$, and the Monte Carlo run $2.722$ from 27,224 particles.
The total volume $M_{1,0,0,0,0}$ is the same before and after in both runs.
The other printed moments agree to within 0.2 percent, except $M_{2,0,0,0,0}$,
which differs by 1.1 percent.
