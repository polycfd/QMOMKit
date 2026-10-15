# 4D CQMOM particles: breakage and recoil

This example evolves a spatially homogeneous particle population in the phase
space $(d,u,v,w)$ with CQMOM and fixed node counts, using
`advance_cqmom_phase_space(...)`. Its breakage source maps each parent to its
daughters nonlinearly, which it evaluates on the quadrature nodes of the
source stage with `CqmomSourceContext::for_each_terminal_node(...)`.
`main.cpp` runs QMOMKit, `main_mc.cpp` the Monte Carlo reference.

## Model

The coordinates are the particle diameter $d$ in m and the particle velocity
components $(u,v,w)$ in m/s. For the number density $n(d,u,v,w)$, the raw
moments are

```math
M_{p,q,r,s}=\int d^p u^q v^r w^s\, n(d,u,v,w)\,\mathrm{d}d\,\mathrm{d}u\,\mathrm{d}v\,\mathrm{d}w,
```

so $M_{0,0,0,0}$ is the number density and $(\pi/6)M_{3,0,0,0}$ the total
particle volume.

The carrier gas has density $\rho_g=1.2$ kg/m³, viscosity
$\mu_g=1.8\times10^{-5}$ Pa s and velocity
$\boldsymbol{U}_g=(2.0,1.5,0.8)$ m/s; the particles have density
$\rho_p=2000$ kg/m³. With the slip
$\boldsymbol{s}=\boldsymbol{U}_g-\boldsymbol{u}$, the particles follow Schiller–Naumann drag,

```math
\mathrm{Re}=\frac{\rho_g d|\boldsymbol{s}|}{\mu_g},\qquad
\frac{\mathrm{d}\boldsymbol{u}}{\mathrm{d}t}=
\frac{3\rho_g}{4\rho_p d}\frac{24}{\mathrm{Re}}
\left(1+0.15\,\mathrm{Re}^{0.687}\right)|\boldsymbol{s}|\boldsymbol{s},
```

with $\mathrm{d}d/\mathrm{d}t=0$.

Each parent breaks at the rate $\beta=10\ \mathrm{s^{-1}}$ into two
equal-volume daughters of diameter $d'=2^{-1/3}d$ with opposite recoil
velocities. Each recoil component is normally distributed with zero mean and
standard deviation

```math
\sigma_r(d)=0.25\frac{d}{10^{-4}\ \mathrm{m}}\ \mathrm{m/s}.
```

The source evaluates the birth-minus-loss contribution to every moment from
the conditional Gaussian moments of the recoil at each quadrature node.
Because the recoil scales with $d$, velocity moments couple to higher diameter
moments: since $\mathbb E[(u+r_u)^3]=u^3+3u\,\sigma_r(d)^2$, the source of
$M_{0,3,0,0}$ contains $M_{2,1,0,0}$, which is not among the moments of the
case. Evaluating the source on the nodes closes it. Breakage alone conserves
the total particle volume and the volume-weighted momentum, increases the
number density and widens the velocity distribution; drag changes the
momentum.

All model parameters are set in [physics.cpp](physics.cpp).
[configuration.qmom](configuration.qmom) selects two nodes per coordinate
(16 nodes) with `adaptive = false`, the identity conditioning order, Strang
splitting and adaptive RK23 node advection. Its `[initial]` section describes
an independent product with $M_{0,0,0,0}=1$: a log-normal diameter with median
100 µm (`sigma_log = 0.12`) and normal velocity components with zero mean and
standard deviation 0.5 m/s. Both executables sample the same 10000 particles
from it, and the QMOMKit run starts from the moments of that sample. Both
advance from $t=0$ to $t=0.1$ s with steps of 0.001 s. The Monte Carlo run
applies half a step of drag, stochastic breakage, and another half step; each
breakage event draws one recoil vector and gives it to the two daughters with
opposite signs.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_4d_particles_breakage_recoil_example
./cqmom_4d_particles_breakage_recoil_mc_example
```

## Output

Both executables print the initial and final number density, first moments,
second moment of each coordinate, and $M_{1,1,0,0}$, which drag couples; the
Monte Carlo executable also prints the initial and final particle counts. The
QMOMKit run writes `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`, the Monte Carlo run
`build/output/moment_history_mc.txt`, each every 0.01 s from $t=0$.

## What to expect

Both runs start from the same moments. Breakage at $10\ \mathrm{s^{-1}}$
multiplies the number density by $e^{10 \cdot 0.1}$, so at $t=0.1$ s the
QMOMKit run gives $M_{0,0,0,0}=2.71828$, and the Monte Carlo run 2.6897 from
26897 particles. Drag accelerates the particles toward the gas velocity, so
the first velocity moments grow from nearly zero. The other printed moments of
the two runs differ by less than 1.3 percent, comparable to the 1 percent
difference in number density that the stochastic breakage of the Monte Carlo
run produces.
