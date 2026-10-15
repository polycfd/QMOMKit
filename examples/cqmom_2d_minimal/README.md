# Minimal 2D CQMOM evolution

This example is the smallest CQMOM time loop: it builds a `CqmomCase` in code,
prepares it as a `CqmomPreparedCase`, creates a `CqmomState` with
`initial_state()`, and advances it with `advance_cqmom_phase_space(...)`, given
a phase-space velocity, a moment source and a boundary contribution. Everything
is set in [main.cpp](./main.cpp).

## Model

The phase-space coordinates $x_0$ and $x_1$ and the time $t$ are
dimensionless, and the raw moments of the number density $n(x_0, x_1)$ are

```math
M_{\boldsymbol{\alpha}}
= \int\!\int x_0^{\alpha_0} x_1^{\alpha_1}\,n(x_0, x_1)\,
\mathrm{d}x_0\,\mathrm{d}x_1.
```

The hierarchy has two nodes per coordinate and conditioning order $\{0, 1\}$.
The initial moment tensor $M_{i,j}$, $i, j = 0, \dots, 3$, holds the moments of
four weighted points of weight $0.25$ at $x_0 \in \{1.5, 2.5\}$ and
$x_1 \in \{-0.5, 0.5\}$. CQMOM reconstructs from ten of its sixteen entries:
`create_cqmom_moment_layout(...)` lists them, and `extract_cqmom_moments(...)`
selects them from the tensor; see
[Moments and quadrature](../../docs/user/4_moments_and_quadrature.md#from-your-own-moments).

The nodes move with

```math
\frac{\mathrm{d}x_0}{\mathrm{d}t} = 1-x_0, \qquad
\frac{\mathrm{d}x_1}{\mathrm{d}t} = -x_1.
```

Each element breaks at the rate $0.5$ into two with half its $x_0$ and the
same $x_1$, and the boundary contribution adds $0.1\,(0.1)^{\alpha_0}$,
independent of $\alpha_1$, which injects number density at the rate $0.1$ at
the point $(x_0, x_1) = (0.1, 1)$:

```math
\frac{\mathrm{d}M_{\boldsymbol{\alpha}}}{\mathrm{d}t}
= 0.5\left(2^{1-\alpha_0}-1\right)M_{\boldsymbol{\alpha}}
+ 0.1\,(0.1)^{\alpha_0}.
```

The run goes from $t = 0$ to $1$ with steps of $0.01$, Strang splitting and
SSP-RK3 node advection, as described in
[Evolving in time](../../docs/user/5_evolving_in_time.md).

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_2d_minimal_example
```

## Output

The example prints three lines: that the run is complete, and the names of the
two files it wrote. `MomentWriter` and `QuadratureWriter` write the moments and
the quadrature every $0.1$ to `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`; see
[Output files](../../docs/user/7_reference.md#output-files). The moment history
holds the ten moments the state carries; `reconstruct_raw_moments(...)`
recovers any other moment from `state.quadrature()`.

## What to expect

$M_{0,0}$ follows $1.2\,e^{0.5t} - 0.2$, from 1 to about 1.78 at $t = 1$.
$M_{0,1}$ starts at 0 and turns positive, to about 0.079 at $t = 1$, because
the injected elements sit at $x_1 = 1$. At $t = 1$ the two $x_0$ nodes lie at
about 0.57 and 1.35.
