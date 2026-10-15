# 2D CQMOM bubble cavitation

This example, based on
[Bryngelson et al. (2023)](https://www.doi.org/10.1016/j.jcp.2023.111917),
moves a spatially uniform population of oscillating bubbles through the phase
space of bubble radius and radial velocity by nonlinear CQMOM node advection.
`main.cpp` loads a `qmomkit::CqmomPreparedCase`, replaces its initial
moments with those of a finite sample through
`qmomkit::config::with_initial_moments(...)` and advances them with
`qmomkit::advance_cqmom_advection(...)`; `main_mc.cpp`, the Monte Carlo
reference, moves the sampled bubbles themselves with
`qmomkit::monte_carlo::evolve_ensemble(...)`.

## Model

The phase-space coordinates are, in this order, the bubble radius $R$ and the
radial velocity $U = \mathrm{d}R/\mathrm{d}t$, and $n(R,U)$ is the number
density; $R$, $U$, time $t$ and all parameters are dimensionless. The raw
moments are $M_{i,j} = \iint R^i U^j\, n(R,U)\,\mathrm{d}R\,\mathrm{d}U$. Each
bubble follows

```math
\frac{\mathrm{d}R}{\mathrm{d}t} = U,
\qquad
\frac{\mathrm{d}U}{\mathrm{d}t}
= \frac{1}{R}\left[
\left(\frac{R_o}{R}\right)^{3\gamma}
- \frac{1}{C_p}
- \frac{2}{\mathrm{We}\,R_o}\left(\frac{R_o}{R}
- \left(\frac{R_o}{R}\right)^{3\gamma}\right)
- \frac{3}{2}U^2
- \frac{4}{\mathrm{Re}}\frac{U}{R}
\right],
```

with the constants set in [physics.hpp](./physics.hpp):

| Symbol | Meaning | Value |
|---|---|---:|
| $R_o$ | reference radius | 1.0 |
| $C_p$ | pressure coefficient | 0.3 |
| $\gamma$ | heat-capacity ratio | 1.4 |
| $\mathrm{Re}$ | Reynolds number | 10.0 |
| $\mathrm{We}$ | Weber number | 1.39 |

The velocity callback is in [physics.cpp](./physics.cpp). Both drivers advance
from $t = 0$ to $t = 13.9$ in steps of $0.01$. Bubbles are only moved, so the
node advection keeps the weights and moves the abscissas.

[configuration.qmom](./configuration.qmom) defines the initial number density
as the total number density $n_0 = 1$ times the product of a log-normal
distribution of $R$ with `mu_log = -0.02` and `sigma_log = 0.2`, whose mean is
$\exp(-0.02 + 0.2^2/2) = 1$, and a normal distribution of $U$ with mean $0$ and
standard deviation $0.2$. It sets `node_counts = 2, 2`, two nodes in $R$ and,
conditioned on each, two in $U$; `adaptive = false` in `[qmom]`, so all four
nodes are kept; and `maximum_orders = 3, 3`, so the sample provides $M_{i,j}$
for $i, j \le 3$. Both executables draw the same $N = 10000$ samples
$(R_p, U_p)$ with a fixed seed ([sampling.hpp](./sampling.hpp)) and give each
the weight $n_0/N$, and the Monte Carlo companion estimates

```math
M_{i,j} = \frac{n_0}{N}\sum_{p=1}^{N} R_p^i U_p^j.
```

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md); this example also needs the
Monte Carlo companion, see its [README](../../tools/monte_carlo/README.md).
From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_2d_bubbles_cavitation_example
./cqmom_2d_bubbles_cavitation_mc_example
```

## Output

Both executables print the initial and final raw moments $M_{0,0}$,
$M_{1,0}$, $M_{0,1}$, $M_{2,0}$, $M_{1,1}$ and $M_{0,2}$. Every $0.1$ in time,
the QMOMKit run writes `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`, and the Monte Carlo run
`build/output/moment_history_mc.txt`. The scripts in
[tools/plotting](../../tools/plotting/README.md) plot these files.

## What to expect

Both runs start from the same sample, so their initial moments agree, and
$M_{0,0}$ stays 1 because bubbles are only moved. At $t = 13.9$, the mean
radius $M_{1,0}$ has fallen from about 1.000 to 0.839 in both runs, and
$M_{2,0}$ agrees to within 0.1 percent. The velocity moments have decayed by
two orders of magnitude, and at that level the runs differ more: $M_{0,2}$ is
$2.8 \times 10^{-4}$ in the QMOMKit run and $3.2 \times 10^{-4}$ in the Monte
Carlo run, while $M_{0,1}$ and $M_{1,1}$ stay near zero in both. The Monte
Carlo run is a stochastic reference with sampling noise, not an exact solution.
