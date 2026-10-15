# QMOM from explicit moments

This example reconstructs a two-node QMOM quadrature with `compute_qmom(...)`
from given raw moments: $m_0 = 1$, $m_1 = 0.8$, $m_2 = 2.8$ and $m_3 = 4.4$,
the moments of the number density
$n(\xi) = 0.4\,\delta(\xi+1) + 0.6\,\delta(\xi-2)$ of a dimensionless
coordinate $\xi$, that is, of weight 0.4 at $\xi = -1$ and weight 0.6 at
$\xi = 2$.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_explicit_moments_example
```

## Output

The example prints $m_0$, the total weight, and the weight and abscissa of each
node. It writes no files.

## What to expect

The quadrature recovers the two weighted points: weights 0.4 and 0.6 at
abscissas -1 and 2.
