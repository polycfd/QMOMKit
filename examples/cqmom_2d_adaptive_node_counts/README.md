# CQMOM: 2D Adaptive Node Counts

This deterministic example evolves a two-dimensional CQMOM phase space ordered
as $(d, u)$, where $d$ is size and $u$ is velocity. It demonstrates adaptive
branch-local node counts: velocity branches reduce as variance collapses,
persist at reduced counts during phase-space advection, recover by
at most one node per accepted moment-source stage, and reduce again after
the source pulse.

Variance collapse and symmetric variance injection are the numerical mechanisms
used to exercise that lifecycle. They are not the sole purpose of the example.

With the default two-velocity-node configuration, each velocity branch follows:

```text
2 -> 1 -> 2 -> 1
```

Exact transition times depend on the configured parameters, reconstruction
tolerances, and output interval. The implementation supports other conditioning
orders and configured maximum counts; focused regression coverage for those
invariants lives in the root QMOMKit test suite.

## Configuration

The key setup is:

```ini
[phase_space]
names = d, u
node_counts = 2, 2
conditioning_order = 0, 1

[qmom]
adaptive = true
```

`load_cqmom_prepared(...)` uses the configured independent-product distribution to
compute the required analytical raw moments. `node_counts = 2, 2` requests two
size nodes and two velocity nodes; with adaptive reconstruction, those counts
are maximums.

Post-source reconstruction automatically tests controlled branch-local growth
back toward the configured maximum. Growth is checked per velocity branch and
may increase by at most one active node per accepted moment-source stage.

The shipped configuration uses `conditioning_order = 0, 1`, so velocity is
conditional on size. The physical velocity coordinate remains phase-space
dimension `u`; conditioning order only changes recursive reconstruction order.

## Evolution

The phase-space velocity in [physics.cpp](./physics.cpp) is:

$$
\begin{aligned}
\frac{\mathrm{d}d}{\mathrm{d}t} &= 0, \\
\frac{\mathrm{d}u}{\mathrm{d}t} &= -a\left(u-u_{\mathrm{eq}}\right).
\end{aligned}
$$

with $a=5$ and $u_{\mathrm{eq}}=0.5$ by default. Therefore the size distribution
remains fixed, while the velocities collapse onto $u_{\mathrm{eq}}$ and their
variance decays as:

```math
\operatorname{Var}[u](t)
= \operatorname{Var}[u](0)\exp(-2at).
```

As the velocity variance approaches zero relative to $u_{\mathrm{eq}}^2$,
adaptive Wheeler reduces active velocity branches toward one node. The nonzero
equilibrium matters: relaxation toward $u_{\mathrm{eq}}=0$ would only shrink the
velocity distribution without changing its shape, and adaptive reduction does
not depend on the unit of the coordinate, so the branches would keep their
nodes.

The moment-source is active only over a fixed pulse interval. For each
parent state $(d,u)$, it applies the deliberately simple splitting event:

$$
(d,u) \to \frac{1}{2}(d,u-\Delta u) + \frac{1}{2}(d,u+\Delta u).
$$

Here $\Delta u$ is the configured velocity displacement and $\lambda$ is the
moment-source rate.

For the raw, unnormalized mixed moment

```math
M_{p,q} = \int\!\int d^p u^q n(d,u)\,\mathrm{d}d\,\mathrm{d}u,
```

the source callback evaluates:

```math
S_{p,q}
= \lambda\,\mathbb{E}\!\left[
d^p\left(\frac{(u-du)^q+(u+du)^q}{2}-u^q\right)
\right].
```

This gives $S_{p,0}=0$ and $S_{p,1}=0$, so the source preserves total weight,
all velocity-independent moments, and velocity mean. For $q=2$, it injects
positive velocity variance through:

```math
S_{p,2} = \lambda\,du^2 M_{p,0}.
```

After the pulse ends, the same advection relaxation again collapses the
velocity variance and adaptive Wheeler reduces recovered branches back toward
one node. This repeated collapse and regrowth is the active-node-count behavior
demonstrated by the example.

## Build and run

From the repository root, install QMOMKit, then configure the example as a
standalone package consumer:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cd examples/cqmom_2d_adaptive_node_counts
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_2d_adaptive_node_counts_example
```

The executable performs concise runtime consistency checks for the demonstrated
adaptive lifecycle and exits with failure if that lifecycle is not obtained.

## Output files

The executable writes the standard deterministic CQMOM files when output is
enabled:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

It also prints output-interval samples with time, source on/off state, velocity
mean, velocity variance, branch-identified active velocity counts, and terminal
paths. Active velocity counts use:

```text
parent-path: active/configured
```

For a marginal velocity branch, the parent path is shown as `root`. Conditional
velocity branches use their actual parent paths, such as `[0]` or `[1]`.

The transition table lists the time, velocity-branch parent path, count change,
and configured maximum.

The implementation is generic in conditioning order and branch identity. The
repository tests verify additional conditioning orders and maximum counts.
