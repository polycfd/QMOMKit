# 2D CQMOM adaptive node counts

This example shows adaptive node counts in a two-dimensional CQMOM state
evolved with `advance_cqmom_phase_space(...)`: each velocity branch loses a
node when its velocity variance collapses, grows one back when a source pulse
injects variance, and loses it again afterwards. It prints every change as
QMOMKit reports it in `result.diagnostics.cqmom_node_count_transitions`. See
[Adaptive node counts](../../docs/user/5_evolving_in_time.md#adaptive-node-counts)
for the rules.

## Model

The phase-space coordinates are the particle size $d$ and the particle
velocity $u$, in that order; $d$, $u$ and the time $t$ are dimensionless. The
raw moments are $M_{p,q} = \iint d^p u^q\,n(d,u)\,\mathrm{d}d\,\mathrm{d}u$.
Node advection keeps the sizes fixed and relaxes the velocities,

```math
\frac{\mathrm{d}d}{\mathrm{d}t} = 0, \qquad
\frac{\mathrm{d}u}{\mathrm{d}t} = -a\left(u-u_\mathrm{eq}\right),
```

so the velocity variance decays as $\exp(-2at)$. As it approaches zero
relative to $u_\mathrm{eq}^2$, the reconstruction after the source reduces
each velocity branch to one node. The nonzero $u_\mathrm{eq}$ matters:
relaxation toward zero would only shrink the velocity distribution without
changing its shape, and adaptive reduction does not depend on the unit of the
coordinate, so the branches would keep their nodes.

During the pulse $t_\mathrm{start} \le t < t_\mathrm{end}$ the source splits
each particle into two of the same size with velocities $u \pm \Delta u$ at
the rate $\lambda$,

```math
(d,u) \to \tfrac{1}{2}(d,u-\Delta u) + \tfrac{1}{2}(d,u+\Delta u),
```

which gives

```math
\frac{\mathrm{d}M_{p,q}}{\mathrm{d}t}
= \lambda \iint d^p\left(\frac{(u-\Delta u)^q+(u+\Delta u)^q}{2}-u^q\right)
  n\,\mathrm{d}d\,\mathrm{d}u
= \lambda \sum_{r = 2, 4, \ldots}^{q} \binom{q}{r}\,\Delta u^r M_{p,q-r}.
```

The source keeps $M_{p,0}$ and $M_{p,1}$, hence the number density, the size
distribution and the mean velocity, and adds velocity variance through
$\mathrm{d}M_{p,2}/\mathrm{d}t = \lambda\,\Delta u^2 M_{p,0}$. The
parameters $a = 5$, $u_\mathrm{eq} = 0.5$, $\lambda = 20$, $\Delta u = 0.10$
and the pulse $[t_\mathrm{start}, t_\mathrm{end}) = [4, 4.05)$ are set in
`VarianceInjectionParameters` in [physics.hpp](./physics.hpp); the callbacks
are defined in [physics.cpp](./physics.cpp).

[configuration.qmom](./configuration.qmom) starts from the independent
product $d \sim \mathcal{N}(2, 0.5^2)$ and $u \sim \mathcal{N}(1, 1)$ with
number density 1, on a hierarchy with `node_counts = 2, 2` and
`conditioning_order = 0, 1`, so velocity is conditioned on size. With
`adaptive = true` the node counts are maximums. [main.cpp](./main.cpp)
advances from $t = 0$ to $8$ with steps of $0.002$.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_2d_adaptive_node_counts_example
```

## Output

The example prints each node-count change of a velocity branch: the time, the
parent path of the branch, such as `[0]` or `[1]`, and the active node counts
before and after. It writes the moments and the terminal nodes of every branch
every $0.5$ to `build/output/moment_history.txt` and
`build/output/quadrature_history.txt`, as set in the `[output]` section of
`configuration.qmom`.

## What to expect

```text
Node-count changes of the velocity branches:
  t = 3.016  branch [0]  2 -> 1 nodes (reduced)
  t = 3.080  branch [1]  2 -> 1 nodes (reduced)
  t = 4.000  branch [0]  1 -> 2 nodes (grew)
  t = 4.000  branch [1]  1 -> 2 nodes (grew)
  t = 6.618  branch [0]  2 -> 1 nodes (reduced)
  t = 6.630  branch [1]  2 -> 1 nodes (reduced)
```

Relaxation reduces both velocity branches to one node shortly after $t = 3$,
the pulse at $t = 4$ restores the second node of each, and relaxation reduces
them again after $t = 6.6$.
