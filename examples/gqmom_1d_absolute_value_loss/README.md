# 1D GQMOM loss proportional to |ξ|

This example evolves the same moments with GQMOM, the `generalized` QMOM
family [[Fox, Laurent & Passalacqua,
2023](https://doi.org/10.1016/j.jaerosci.2022.106096)], and with classical
QMOM, under a source that removes the elements of a population at a rate
proportional to $|\xi|$, and compares both with the exact solution. The source
is evaluated on the nodes of the quadrature, where the many nodes of GQMOM
represent the number density far better than the few nodes of classical QMOM.

## Model

The phase-space coordinate $\xi$ and time $t$ are dimensionless, and $\xi$
takes any real value. In a spatially homogeneous population, each element is
lost at the rate $c\,|\xi|$, with $c = 1$, so the number density $n(t, \xi)$
obeys $\partial_t n = -c\,|\xi|\,n$ and its moments change as

```math
\frac{\mathrm{d}m_k}{\mathrm{d}t}
= -c\int |\xi|\,\xi^k\,n(t, \xi)\,\mathrm{d}\xi
\approx -c\sum_i w_i\,|\xi_i|\,\xi_i^k,
```

where $w_i$ and $\xi_i$ are the weights and abscissas of the quadrature that
QMOMKit reconstructs from the moments of each Runge–Kutta stage. $|\xi|$ is not
a polynomial, so the right-hand side is not a combination of the moments the
state holds. The exact solution is $n(t, \xi) = n(0, \xi)\,e^{-c|\xi|t}$:
elements far from $\xi = 0$ are lost first, so both the number density $m_0$
and the mean of $\xi^2$, $m_2/m_0$, fall. For a normal $n(0, \xi)$, its moments
follow in closed form. [physics.cpp](./physics.cpp) defines $c$, the source
and these exact moments.

[configuration.qmom](./configuration.qmom) defines the standard normal
distribution as $n(0, \xi)$, with node count $N = 2$, so that the state holds
$m_0, \ldots, m_4$, and the `generalized` family with $N_G = 20$ nodes
(`generalized_node_counts`). [main.cpp](./main.cpp) evolves this case and the
same case with the `classical` family, which has $N$ nodes, with source steps
only, from $t = 0$ to $2$ with steps of $0.01$. The loss changes the number of
elements, not their coordinates, so there is no node advection.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./gqmom_1d_absolute_value_loss_example
```

## Output

The example prints $m_0$ and $m_2/m_0$ at $t = 0$, and at $t = 2$ those of the
exact solution, of GQMOM and of classical QMOM. It writes the moments and the
quadratures of both states every $0.5$ to `build/output/gqmom_moments.txt`,
`qmom_moments.txt`, `gqmom_quadrature.txt` and `qmom_quadrature.txt`, as set in
the `[output]` section of `configuration.qmom`.

## What to expect

At $t = 2$, the example prints:

```text
          m0          m2/m0
  exact   0.336204    0.253569
  GQMOM   0.320518    0.264196
  QMOM    0.135335    1
```

GQMOM is within 5% of the exact $m_0$ and $m_2/m_0$. Classical QMOM ends with
60% fewer elements than the exact solution, and its $m_2/m_0$ does not change
at all: its two nodes lie at $\xi = \pm 1$, have the same $|\xi|$, and so lose
their weights at the same rate. The 20 nodes of GQMOM resolve the elements
near $\xi = 0$, which survive the longest. Raising `generalized_node_counts` in
`build/configuration.qmom` brings GQMOM closer to the exact values, at a higher
cost per step.
