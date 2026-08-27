# CQMOM Independent Product Initialization

This example loads an independent-product CQMOM setup from configuration and
shows the complete normal initialization lifecycle:

```text
configuration
    -> load_cqmom(...)
    -> CqmomCase
    -> compute_cqmom(...)
    -> quadrature nodes and weights

CqmomCase
    -> CqmomCase::initial_state()
    -> CqmomState
```

The shipped configuration is two-dimensional and orders its phase space as
`(size, velocity)`. It requests a fixed `2 x 2` CQMOM hierarchy with
`conditioning_order = 0, 1`: two size nodes and two conditional velocity
nodes. The initial number density is `2.5`, with a gamma distribution in
`size` (shape `3.0`, scale `0.2`) and a normal distribution in `velocity`
(mean `-0.5`, standard deviation `1.25`).

For the number density $n(\mathrm{size},\mathrm{velocity})$, the raw multivariate moments are

$$
M_{i,j} = \int\!\int \mathrm{size}^i\mathrm{velocity}^j\,
n(\mathrm{size},\mathrm{velocity})\,
\mathrm{d}\mathrm{size}\,\mathrm{d}\mathrm{velocity}.
$$

They are unnormalized, so $M_{0,0}$ is the total number density `2.5`.
`compute_cqmom(...)` converts the configured moment tensor into the initial
conditional quadrature. The executable then creates a `CqmomState` to show
the state lifecycle but does not add phase-space advection, a moment source,
or a boundary contribution.

`CqmomCase` is the immutable configuration and problem definition.
`CqmomState` is the runtime state that a phase-space evolution workflow evolves later.
This example creates the state but does not add phase-space evolution physics.

Requires: QMOMKit.

Required: `configuration.qmom`, `main.cpp`, and an installed QMOMKit package.
Optional: adapt the distribution, `node_counts`, or `conditioning_order` in
`configuration.qmom`.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_independent_product_initialization_example
```

The executable writes no files. It prints the hierarchy, conditioning order,
and reconstructed node weights and abscissas. The `x0` and `x1` labels in the
terminal output are example-specific formatting for the shipped two-dimensional
case: `x0` is `size` and `x1` is `velocity`. They are not a limitation of
QMOMKit's multidimensional representation.

## What to expect

The shipped fixed `2 x 2` hierarchy produces four terminal CQMOM nodes. Their
global weights sum to the configured total number density `2.5`, while their
coordinates represent the independent `size` and `velocity` initialization.
