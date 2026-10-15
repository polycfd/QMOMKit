# 5D CQMOM particles: breakage and heat transfer

This example evolves a spatially homogeneous particle population in the phase
space $(V,T,u,v,w)$ with CQMOM and fixed node counts, using
`advance_cqmom_phase_space(...)` for drag, heat transfer and breakage. Its
breakage source is written in moment space and needs moments of fractional
order in $V$, which it evaluates on the quadrature of the source stage with
`CqmomSourceContext::generalized_moment(...)`. The
[4D breakage and recoil example](../cqmom_4d_particles_breakage_recoil/README.md)
evaluates its source on the quadrature nodes instead. `main.cpp` runs QMOMKit,
`main_mc.cpp` the Monte Carlo reference.

## Model

The coordinates are the particle volume $V$ in m³, the particle temperature
$T$ in K and the particle velocity components $(u,v,w)$ in m/s. For the number
density $n(V,T,u,v,w)$, the raw moments are

```math
M_{p,q,r,s,t}=\int V^pT^qu^rv^sw^t\, n(V,T,u,v,w)\,\mathrm{d}V\,\mathrm{d}T\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w,
```

so $M_{0,0,0,0,0}$ is the number density and $M_{1,0,0,0,0}$ the total
particle volume.

The carrier gas has density $\rho_g=1.2$ kg/m³, viscosity
$\mu_g=1.8\times10^{-5}$ Pa s, velocity $\boldsymbol{U}_g=(2.0,1.5,0.8)$ m/s,
thermal conductivity $k_g=0.026$ W/(m K), Prandtl number $\mathrm{Pr}=0.71$
and temperature $T_g=300$ K. The particles have density $\rho_p=2000$ kg/m³,
heat capacity $c_p=1000$ J/(kg K) and diameter $d=(6V/\pi)^{1/3}$. With the
slip $\boldsymbol{s}=\boldsymbol{U}_g-\boldsymbol{u}$, they follow
Schiller–Naumann drag and Ranz–Marshall heat transfer,

```math
\begin{aligned}
\mathrm{Re}&=\frac{\rho_g d|\boldsymbol{s}|}{\mu_g},\qquad
\frac{\mathrm{d}\boldsymbol{u}}{\mathrm{d}t}=
\frac{3\rho_g}{4\rho_p d}\frac{24}{\mathrm{Re}}
\left(1+0.15\,\mathrm{Re}^{0.687}\right)|\boldsymbol{s}|\boldsymbol{s}, \\
\mathrm{Nu}&=2+0.6\,\mathrm{Re}^{1/2}\mathrm{Pr}^{1/3},\qquad
\frac{\mathrm{d}T}{\mathrm{d}t}=
\frac{6k_g\mathrm{Nu}}{\rho_p c_p d^2}(T_g-T),
\end{aligned}
```

with $\mathrm{d}V/\mathrm{d}t=0$.

Each parent breaks at the rate $\beta=10\ \mathrm{s^{-1}}$ into two
equal-volume daughters that keep its velocity and move their temperature
slightly toward the gas temperature:

```math
V_d=\frac{V_p}{2},\qquad
T_d=T_p+\varepsilon\left(\frac{V_d}{V_{\mathrm{ref}}}\right)^{1/3}(T_g-T_p),
```

with $\varepsilon=0.01$ (dimensionless) and
$V_{\mathrm{ref}}=\pi(10^{-4}\ \mathrm{m})^3/6$. Expanding $T_d^q$ gives
powers $V^{1/3}$, so the source of a moment needs moments of fractional order
in $V$. With $a=\varepsilon(2V_{\mathrm{ref}})^{-1/3}$, for example,

```math
\dot M_{1,1,0,0,0}\big|_{\mathrm{breakage}}
=\beta a\left(T_g M_{4/3,0,0,0,0}-M_{4/3,1,0,0,0}\right).
```

The source takes integer orders from `CqmomSourceContext::moment(...)` and
fractional orders from `CqmomSourceContext::generalized_moment(...)`; the
fractional moments are not added to the moments of the case. For $q=0$ the
source reduces to $\beta(2^{1-p}-1)M_{p,0,r,s,t}$, so breakage conserves the
total particle volume.

All model parameters are set in [physics.cpp](physics.cpp).
[configuration.qmom](configuration.qmom) selects two nodes per coordinate
(32 nodes) with `adaptive = false`, the identity conditioning order, Strang
splitting and adaptive RK23 node advection. Its `[initial]` section describes
an independent product with $M_{0,0,0,0,0}=1$: a log-normal volume, which is
the diameter distribution of the 4D example (median 100 µm) transformed to
$V=\pi d^3/6$, a normal temperature with mean 350 K and standard deviation
50 K, and normal velocity components with zero mean and standard deviation
0.5 m/s. Both executables sample the same 10000 particles from it, and the
QMOMKit run starts from the moments of that sample. Both advance from $t=0$ to
$t=0.1$ s with steps of 0.001 s; the Monte Carlo run applies half a step of
drag and heat transfer, stochastic breakage, and another half step.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_5d_particles_breakage_heat_transfer_example
./cqmom_5d_particles_breakage_heat_transfer_mc_example
```

The Monte Carlo run takes noticeably longer than the QMOMKit run, because it
follows every particle through breakage.

## Output

Both executables print the initial and final number density, first moments,
second moment of each coordinate, and $M_{1,1,0,0,0}$, which size-dependent
heating couples; the Monte Carlo executable also prints the initial and final
particle counts. The QMOMKit run writes `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`, the Monte Carlo run
`build/output/moment_history_mc.txt`. Each is written at $t=0$ and at the
first step that reaches each multiple of the configured interval of 0.0025 s.

## What to expect

Breakage at $10\ \mathrm{s^{-1}}$ multiplies the number density by
$e^{10 \cdot 0.1}$, so at $t=0.1$ s the QMOMKit run gives
$M_{0,0,0,0,0}=2.71828$, and the Monte Carlo run 2.7224 from 27224 particles.
The total volume $M_{1,0,0,0,0}$ is the same before and after in both runs.
The particles cool toward the gas temperature: the mean temperature
$M_{0,1,0,0,0}/M_{0,0,0,0,0}$ falls from 350 K to about 302 K. The other
printed moments of the two runs differ by less than 0.2 percent, except
$M_{2,0,0,0,0}$, which differs by 1.1 percent; the Monte Carlo moments carry
the sampling noise of a finite sample with stochastic breakage.
