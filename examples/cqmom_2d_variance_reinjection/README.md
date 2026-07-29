# CQMOM: 2D Variance Reinjection

This deterministic example evolves a two-dimensional CQMOM phase space ordered
as `(d, u)`, where `d` is size and `u` is velocity. It is based closely on
`examples/cqmom_2d_vanishing_variance`, but adds a finite-time symmetric
velocity-redistribution source so adaptive velocity branches reduce, regrow
one node per accepted source step, and reduce again. With the default
two-velocity-node configuration, each velocity branch follows:

```text
2 -> 1 -> 2 -> 1
```

It's also possible to change the conditioning order or increase the maximum node count.

## Configuration

The key setup is:

```text
[phase_space]
names = d, u
node_counts = 2, 2
conditioning_order = 0, 1
adaptive_branch_counts = true

[time_integration.cqmom]
source_node_count_evolution = allow_growth_to_configured_maximum

[qmom]
adaptive = true
```

`load_cqmom(...)` uses the configured independent-product distribution to
compute the required analytical raw moments. `node_counts = 2, 2` requests two
size nodes and two velocity nodes; with adaptive branch counts enabled, those
counts are maximums. The example test configurations also exercise
`node_counts = 2, 3`.

The source-stage policy lets post-source reconstruction test controlled
branch-local growth back toward the configured maximum. Growth is validated per
velocity branch and may increase by at most one active node per accepted source
step.

Two conditioning orders are supported for the 2D problem:

```text
conditioning_order = 0, 1
conditioning_order = 1, 0
```

The physical velocity coordinate remains phase-space dimension `u` in either
case; only the recursive conditioning order changes.

## Evolution

The continuous phase-space velocity in `physics.cpp` is:

```text
dd/dt = 0
du/dt = -a * (u - u_eq)
```

with `a = 5` and `u_eq = 0` by default. Therefore the size distribution remains
fixed, while the velocity variance decays as:

```text
Var[u](t) = Var[u](0) exp(-2*a*t)
```

As the velocity variance approaches zero, adaptive Wheeler reduces active
velocity branches toward one node.

The discontinuous source is active only over a fixed pulse interval. For each
parent state `(d,u)`, it applies the deliberately simple numerical validation
event:

```text
(d,u) -> 0.5 * (d,u-du) + 0.5 * (d,u+du)
```

For a raw mixed moment `M[p,q]`, the source callback evaluates:

```text
S[p,q] = lambda * E[d^p * (((u-du)^q + (u+du)^q)/2 - u^q)]
```

This gives `S[p,0] = 0` and `S[p,1] = 0`, so the source preserves total weight,
all velocity-independent moments, and velocity mean. For `q = 2`, it injects
positive velocity variance through:

```text
S[p,2] = lambda * du^2 * M[p,0]
```

After the pulse ends, the same continuous relaxation again collapses the
velocity variance and adaptive Wheeler reduces recovered branches back toward
one node.

## Build And Run

From the repository root, build the example with the in-tree target:

```bash
cmake -S . -B build -DQMOMKIT_BUILD_EXAMPLES=ON
cmake --build build --target cqmom_2d_variance_reinjection_example
cd build/examples/cqmom_2d_variance_reinjection
./cqmom_2d_variance_reinjection_example
```

To run the example validation variants through CTest:

```bash
cmake -S . -B build -DQMOMKIT_BUILD_TESTS=ON -DQMOMKIT_BUILD_EXAMPLES=ON
cmake --build build --target cqmom_2d_variance_reinjection_example
ctest --test-dir build -R cqmom_2d_variance_reinjection --output-on-failure
```

The standalone installed-package workflow also works:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cd examples/cqmom_2d_variance_reinjection
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_2d_variance_reinjection_example
```

## Output

The executable writes the standard deterministic CQMOM files when output is
enabled:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

It also prints output-interval samples with time, source on/off state, velocity
mean, velocity variance, branch-identified active velocity counts, and terminal
paths. Active-count transition lines include conditioning level, phase-space
dimension, parent path, count change, configured maximum, and a centered
rank-indicator diagnostic when velocity moments through fourth order are
available.

The implementation is generic in conditioning order and branch identity. The
repository tests currently validate velocity maxima two and three.
