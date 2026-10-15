# QMOM initialization from a mixture

This example loads a QMOM case with `load_qmom(...)` from
[configuration.qmom](./configuration.qmom), in which QMOMKit generates the raw
moments of a two-component mixture, and reconstructs a three-node quadrature
with `compute_qmom(...)`.

## Model

The number density $n(\xi)$ of a dimensionless coordinate $\xi$ is the sum of

- a gamma distribution with number density $1.5$, shape $2.0$ and scale $0.4$;
- a log-normal distribution with number density $0.8$, logarithmic mean
  $-0.2$ and logarithmic standard deviation $0.5$.

`type = mixture` selects the mixture, and each
`[initial.mixture.component.K]` section sets one component. Component number
densities are absolute, not fractions, so each component's moments are scaled
by its number density and added, and $m_0 = 1.5 + 0.8 = 2.3$. Three nodes
(`node_counts = 3`) need $m_0, \dots, m_5$; `maximum_order = 5` requests
them. See the [Reference](../../docs/user/7_reference.md#initial) for the
keys.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_mixture_initialization_example
```

## Output

The example prints $m_0$ and the weight and abscissa of each node. It writes no
files.

## What to expect

$m_0 = 2.3$. The quadrature has weights of about 1.37, 0.90 and 0.03 at
abscissas of about 0.44, 1.38 and 3.34, and the weights sum to $m_0$.
