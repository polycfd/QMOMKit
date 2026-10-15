# CQMOM initialization from an independent product

This example loads a two-dimensional CQMOM case with `load_cqmom_prepared(...)`
from [configuration.qmom](./configuration.qmom), creates its `CqmomState` with
`initial_state()`, and prints the quadrature from `state.quadrature()`.

## Model

The phase space is `(size, velocity)`; the names are labels, and both
coordinates are dimensionless. `type = independent_product` makes the number
density $n(\mathrm{size}, \mathrm{velocity})$ the product of a gamma
distribution in size (shape $3.0$, scale $0.2$) and a normal distribution in
velocity (mean $-0.5$, standard deviation $1.25$), scaled once by
`number_density = 2.5`, so the raw moment $M_{0,0} = 2.5$.
`node_counts = 2, 2` asks for two size nodes and two velocity nodes, and
`conditioning_order = 0, 1` reconstructs size first and velocity conditioned
on it; see
[Conditioning order](../../docs/user/4_moments_and_quadrature.md#conditioning-order).
With `conditioning_order = 1, 0`, velocity is reconstructed first: for this
independent product the same four weighted points result, listed with velocity
as the outer branch, and `x0` still denotes size and `x1` velocity.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./cqmom_independent_product_initialization_example
```

## Output

The example prints the weight and coordinates of each node, with `x0` the size
and `x1` the velocity. It writes no files.

## What to expect

Four nodes: size nodes at 0.4 and 1.2, each paired with velocity nodes at
-1.75 and 0.75. The weights are 0.9375 at size 0.4 and 0.3125 at size 1.2, and
sum to $M_{0,0} = 2.5$.
