# CQMOM: 2D Vanishing Velocity Variance

This deterministic example evolves a two-dimensional CQMOM phase space ordered
as `(d, u)`, where `d` is size and `u` is velocity. The initial distributions
are defined in `configuration.qmom`; QMOMKit reads those analytical
distributions and generates the CQMOM moments needed by the configured
hierarchy.

## Configuration

The key setup is:

```text
[phase_space]
names = d, u
node_counts = 2, 2
conditioning_order = 0, 1
adaptive_branch_counts = true

[time_integration.cqmom]
source_node_count_evolution = preserve_active_counts

[qmom]
adaptive = true

[initial]
type = independent_product
number_density = 1.0

[initial.dimension.0]
family = normal
mean = 2.0
standard_deviation = 0.5

[initial.dimension.1]
family = normal
mean = 1.0
standard_deviation = 1.0
```

`load_cqmom(...)` uses the configured independent-product distribution to
compute the required analytical raw moments. `node_counts = 2, 2` requests two
size nodes and two velocity nodes; with adaptive branch counts enabled, those
counts are maximums.

To check the alternate conditioning order, change only:

```text
conditioning_order = 1, 0
```

## Evolution

The continuous phase-space velocity in `physics.cpp` is:

```text
dd/dt = 0
du/dt = -a*u
```

with $a=5$ by default. Therefore the size distribution remains fixed, while:

```text
E[u](t) = E[u](0) exp(-a*t)
Var[u](t) = Var[u](0) exp(-2*a*t)
```

As the velocity variance approaches zero, adaptive Wheeler reduces the
velocity quadrature, while `preserve_active_counts` prevents that no-op source stage from growing a reduced branch back to its configured maximum.

For `conditioning_order = 0, 1`, the final paths should show two size roots
with one velocity child each:

```text
[0,0], [1,0]
```

For `conditioning_order = 1, 0`, the final paths should show one velocity root
with two size children:

```text
[0,0], [0,1]
```

Both structures represent the same physical result: velocity has collapsed to
one active node and size remains represented by two nodes.

The example can, of course, also be run with larger initial node counts, the end result is the same.

## Build And Run

Install QMOMKit into the repository-local prefix, then configure the example:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cd examples/cqmom_2d_vanishing_variance
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build 
```

Run from the build directory:

```bash
cd examples/cqmom_2d_vanishing_variance/build
./cqmom_2d_vanishing_variance_example
```

## Output

The executable writes the standard deterministic CQMOM files when output is
enabled:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

It also prints selected initial and final moments, a compact size/velocity
variance history reconstructed from the quadrature history, and the final
terminal paths. `quadrature_history.txt` is the most direct file for observing
the adaptive reduction because it has one row per active terminal node at each
output time.
