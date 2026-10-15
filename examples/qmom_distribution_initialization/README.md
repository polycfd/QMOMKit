# QMOM initialization from a distribution

This example generates the raw moments of a log-normal distribution with
`compute_raw_moments(...)` and reconstructs a two-node QMOM quadrature from
them with `compute_qmom(...)`, as described in
[Moments and quadrature](../../docs/user/4_moments_and_quadrature.md#from-distributions).

## Model

The number density $n(\xi)$ of a dimensionless coordinate $\xi$ is the
log-normal distribution `LogNormalDistribution{-0.2, 0.4}`, with logarithmic
mean $-0.2$ and logarithmic standard deviation $0.4$, scaled to the number
density $2.5$. The raw moments are $m_k = \int \xi^k\,n(\xi)\,\mathrm{d}\xi$,
so $m_0 = 2.5$. A two-node quadrature needs $m_0, \dots, m_3$, so
[main.cpp](./main.cpp) requests the moments through order
$2 \cdot 2 - 1 = 3$.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_distribution_initialization_example
```

## Output

The example prints the number density, the moments $m_0, \dots, m_3$, and the
weight and abscissa of each node. It writes no files.

## What to expect

The moments are about 2.5, 2.22, 2.31 and 2.82. The quadrature has weights of
about 1.94 and 0.56 at abscissas of about 0.69 and 1.57; the weights sum to
$m_0 = 2.5$.
