# 3D CQMOM droplets: condensation and nucleation

This example evolves a spatially uniform population of water droplets in the
three-dimensional phase space $(r,w,T)$ with adaptive CQMOM. Condensation,
drag and heat transfer move the droplets through phase space, and nucleation
adds droplets at a fixed birth state through a
`CqmomMomentBoundaryContribution`, composed with node advection by
`advance_cqmom_phase_space(...)`. The initial moments are those of a finite
sample of droplets, passed to the case with `config::with_initial_moments(...)`.
`main.cpp` runs QMOMKit, `main_mc.cpp` the Monte Carlo reference, which evolves
the same sample directly.

## Model

The coordinates are the droplet radius $r$ in m, the vertical droplet velocity
$w$ in m/s (positive upward, so falling droplets have $w<0$) and the droplet
temperature $T$ in K. For the number density $n(r,w,T)$, the raw moments are

```math
M_{\alpha_r,\alpha_w,\alpha_T}
= \iiint r^{\alpha_r} w^{\alpha_w} T^{\alpha_T}\, n(r,w,T)
\,\mathrm{d}r\,\mathrm{d}w\,\mathrm{d}T,
```

and $M_{0,0,0}$ is the total number density.

The droplet model is a simplified version of the model used by
[Denner (2026)](https://www.doi.org/10.1016/j.ijmultiphaseflow.2026.105696)
for water droplets in air. The droplets are isolated and exchange mass,
momentum and heat with a prescribed ambient air state, which they do not
change. With the relative velocity $w_{\mathrm{rel}}=w_\infty-w$, each droplet
follows

```math
\begin{aligned}
\frac{\mathrm{d}r}{\mathrm{d}t}
&= \frac{\mathrm{Sh}\,D_v\,M_w}{2r\rho_l R_g}
\left(\phi\frac{p_{\mathrm{sat}}(T_\infty)}{T_\infty}
- \frac{p_{\mathrm{sat}}(T)}{T}\right), \\
\frac{\mathrm{d}w}{\mathrm{d}t}
&= -g\left(1-\frac{\rho_a}{\rho_l}\right)
+ \frac{3\rho_a C_d\,|w_{\mathrm{rel}}|\,w_{\mathrm{rel}}}{8\rho_l r}, \\
\frac{\mathrm{d}T}{\mathrm{d}t}
&= \frac{3}{r c_l}\left[
\frac{\mathrm{Nu}\,k_a\,(T_\infty-T)}{2r\rho_l}
+ L(T)\frac{\mathrm{d}r}{\mathrm{d}t}\right].
\end{aligned}
```

The radius rate is driven by the vapor difference between the ambient air and
the droplet surface, the velocity rate by buoyancy-corrected gravity and drag,
and the temperature rate by heat transfer to the ambient air and the latent
heat of condensation. The ambient state and the material constants are set in
[droplet_model.cpp](droplet_model.cpp):

| Symbol | Value | Unit | Meaning |
|---|---:|---|---|
| $T_\infty$ | 313.15 | K | ambient air temperature |
| $\phi$ | 1.0 | dimensionless | ambient relative humidity |
| $w_\infty$ | 0.0 | m/s | vertical air velocity |
| $g$ | 9.81 | m/s² | gravitational acceleration |
| $\rho_l$ | 997.0 | kg/m³ | liquid-water density |
| $c_l$ | 4181.0 | J/(kg K) | liquid-water heat capacity |
| $\sigma$ | 0.072 | N/m | surface tension |
| $M_w$ | 18.0153e-3 | kg/mol | water molar mass |
| $R_g$ | 8.3145 | J/(mol K) | universal gas constant |
| $T_c$ | 647.096 | K | water critical temperature |

The remaining quantities are correlations, implemented in
[droplet_model.cpp](droplet_model.cpp):

- moist-air density $\rho_a$, dynamic viscosity $\mu_a$, thermal conductivity
  $k_a$ and heat capacity $c_a$: quadratic polynomials in $T_\infty$ and
  $\phi$, evaluated once at the ambient state;
- saturation vapor pressure $p_{\mathrm{sat}}(T)$ in Pa, water-vapor
  diffusivity $D_v$ in m²/s at $T_\infty$, and latent heat $L(T)$ in J/kg;
- Sherwood number $\mathrm{Sh}$ and Nusselt number $\mathrm{Nu}$, from the
  transfer-number correlation of
  [Beard & Pruppacher (1971)](https://www.doi.org/10.1175/1520-0469(1971)028<1455:AWTIOT>2.0.CO;2)
  with the Reynolds number $\mathrm{Re}=2r\rho_a|w_{\mathrm{rel}}|/\mu_a$ and
  no wet-bulb correction;
- drag coefficient $C_d=C_{d,0}F_{\mathrm{def}}$, with the sphere drag
  $C_{d,0}$ of [Clift & Gauvin (1971)](https://www.doi.org/10.1002/cjce.5450490403)
  and the deformation factor $F_{\mathrm{def}}$ of
  [Lie et al. (1993)](https://saemobilus.sae.org/papers/modeling-effects-drop-drag-breakup-fuel-sprays-930072)
  and [Hsiang & Faeth (1992)](https://www.doi.org/10.1016/0301-9322(92)90036-G),
  a function of the Weber number
  $\mathrm{We}=2r\rho_a w_{\mathrm{rel}}^2/\sigma$;
  below $\mathrm{Re}=10^{-8}$ Stokes drag is used.

The correlations hold for ambient air temperatures in $[273.15,373.15]$ K and
droplet temperatures in $[273.15,T_c)$ K, which this case stays within.

Nucleation is active for $0.20\ \mathrm{s} \le t < 0.40\ \mathrm{s}$ with the
constant birth flux

```math
J = \frac{M_{0,0,0}(0)}{0.40\ \mathrm{s}-0.20\ \mathrm{s}},
```

in number density per second, so it adds one initial population. Every newborn
droplet has the state
$(r_n,w_n,T_n)=(50\times10^{-6}\ \mathrm{m},\ 0\ \mathrm{m/s},\ 298.15\ \mathrm{K})$,
so the boundary contribution to each moment is

```math
\left.\frac{\mathrm{d}M_{a,b,c}}{\mathrm{d}t}\right|_{\mathrm{nuc}}
= J\, r_n^a w_n^b T_n^c.
```

The nucleation parameters are set in [droplet_model.hpp](droplet_model.hpp).

[configuration.qmom](configuration.qmom) selects up to 3 nodes per coordinate
(`node_counts = 3, 3, 3`, at most 27 nodes) with `adaptive = true`, the
conditioning order $r$, $T$, $w$ (`conditioning_order = 0, 2, 1`), Strang
splitting and adaptive RK23 node advection. Its `[initial]` section describes
the initial distribution with $M_{0,0,0}=1$: a gamma distribution in $r$ with
mean 300 µm, and normal distributions in $w$ (mean $-2$ m/s, standard deviation
1 m/s) and $T$ (mean 293.15 K, standard deviation 3 K). Both executables draw
the same sample of 10000 droplets from it ([sampling.hpp](sampling.hpp)), and
the QMOMKit run starts from the moments of that sample.

`main.cpp` advances from $t=0$ to $t=1$ s with steps of 0.01 s, cut at 0.20 s
and 0.40 s so that no step crosses the start or end of nucleation. While
nucleation is active, each step is `advance_cqmom_phase_space(...)`, which
composes node advection with the boundary contribution and a zero volumetric
source; at all other times it is `advance_cqmom_advection(...)`. The Monte
Carlo run takes the same steps and, while nucleation is active, adds
$J\,\Delta t$ of number density between two advection half-steps with
`monte_carlo::inject_population(...)`, which draws a Poisson-distributed
number of new droplets.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_3d_droplets_condensation_nucleation_example
./cqmom_3d_droplets_condensation_nucleation_mc_example
```

## Output

Both executables print the initial and final values of $M_{0,0,0}$ and of the
first and second moments of $r$, $w$ and $T$; the Monte Carlo executable also
prints the initial and final droplet counts. The QMOMKit run writes
`build/output/moment_history.txt` and `build/output/quadrature_history.txt`,
the Monte Carlo run `build/output/moment_history_mc.txt`, each every 0.1 s
from $t=0$. The moment files have the same columns, and their initial rows
agree to round-off because both runs start from the same sample.
[tools/plotting](../../tools/plotting/README.md) plots these files.

## What to expect

The ambient air is saturated and warmer than the droplets, so the droplets grow
by condensation. Their temperatures relax toward the ambient 313.15 K while the
latent heat couples the radius change back into the temperature. Their
velocities evolve under gravity and drag, which depend on $r$. Between 0.20 s
and 0.40 s nucleation doubles the number density, and the newborn droplets
condense as well.

At $t=1$ s the QMOMKit run gives $M_{0,0,0}=2$ and the Monte Carlo run 1.9947
from 19947 droplets. The mean temperature $M_{0,0,1}/M_{0,0,0}$ is about 313 K
in both runs, and the other printed moments differ by less than 0.3 percent,
comparable to the difference in number density that the Poisson-distributed
births of the Monte Carlo run produce.
