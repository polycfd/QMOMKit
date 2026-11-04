# HyQMOM characteristic speeds

This example computes the characteristic speeds of HyQMOM, the `hyperbolic`
QMOM family, with `compute_hyqmom_characteristic_speeds(...)`, for the moments
of section 5 and Fig. 1 of [Fox & Laurent
(2022)](https://www.doi.org/10.1137/21M1406143). When the moments are
transported in physical space, with $\xi$ the velocity, the characteristic
speeds are the speeds at which changes in the moments travel. A code that
transports them needs the slowest and the fastest, for example for the fluxes
and the time step of a finite-volume scheme, as in section 6 of the paper. That
transport belongs to the host code, so it is not part of this example.
Everything is set in [main.cpp](./main.cpp).

## Moments

The velocity distributions have $m_0 = 1$, mean 0 and variance 1, so their raw
moments are the standardized moments $S_k$ of the paper. For each node count
$n$, $S_3, \ldots, S_{2n-1}$ are fixed and $m_{2n}$ varies with
$H_{2n} = \det(m_{i+j})_{i,j=0}^{n}$, which measures how far the moments lie
from the boundary of the realizable set, on which $H_{2n} = 0$. The ranges are
those plotted in Fig. 1 of the paper:

| $n$ | $S_3, \ldots, S_{2n-1}$ | Plotted range of $H_{2n}$ |
| --- | --- | --- |
| 2 | $-1$ | 0 to 20 |
| 3 | $-1, 5, -8$ | 0 to 700 |
| 4 | $-1, 5, -8, 67.3, -100$ | 0 to $10^6$ |
| 5 | $-1, 5, -8, 67.3, -100, 3000, 0$ | 0 to $10^{11}$ |

On the boundary, $m_{2n}$ is that of the classical quadrature of $n$ nodes of
$m_0, \ldots, m_{2n-1}$, which `compute_qmom(...)` in fixed mode and
`reconstruct_raw_moments(...)` give. Away from it,

```math
m_{2n} = m_{2n}^{\mathrm{boundary}} + \frac{H_{2n}}{H_{2n-2}},
\qquad
H_{2n-2} = \prod_{k=1}^{n-1}\,\prod_{j=1}^{k} b_j,
```

with $b_j$ the recurrence coefficients of `compute_wheeler_recurrence(...)`.
For 101 values of $H_{2n}$ on a uniform grid, the example creates the state of
$m_0, \ldots, m_{2n}$ with `state_candidate(...)` of a `QmomCase` of the
`hyperbolic` family with node count $n$, as a host code does for moments it
computed itself, and computes its speeds.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./hyqmom_characteristic_speeds_example
```

## Output

For each $n$, the example prints $m_{2n}$ on the boundary, $H_{2n-2}$, and the
node count and speeds of the states at $H_{2n} = 0$ and at the largest
$H_{2n}$. It writes `build/output/speeds_n2.txt` to `speeds_n5.txt`,
whitespace-delimited tables with a commented header line: $H_{2n}$, the node
count of the state, and its $2n + 1$ speeds in ascending order. Fig. 1 of the
paper plots the speeds of moments with mean 0 and variance 1, so plotting the
speeds of a file against its first column reproduces a panel of the figure, for
example with `numpy` and `matplotlib`, from `build`:

```bash
python3 -c "
import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt('output/speeds_n3.txt')
plt.plot(d[:, 2:], d[:, 0], 'k.', markersize=2)
plt.xlabel('characteristic speed'); plt.ylabel('H_6')
plt.savefig('speeds_n3.png')"
```

## What to expect

Every state has $2n + 1$ speeds. For $H_{2n} > 0$, its quadrature has $n + 1$
nodes. The speeds in even positions, the second, the fourth and so on, are
then the abscissas of the classical quadrature of $n$ nodes of $m_0, \ldots, m_{2n-1}$, which do not depend on
$H_{2n}$, and the other $n + 1$ speeds lie between and beyond them (Theorem 4.8
of the paper):

| $n$ | Abscissas of the classical quadrature | Slowest and fastest speeds at the largest $H_{2n}$ |
| --- | --- | --- |
| 2 | $-1.6180$, $0.61803$ | $-7.8914$, $6.4013$ |
| 3 | $-2.4732$, $0.077968$, $1.7286$ | $-13.537$, $13.635$ |
| 4 | $-3.9691$, $-1.3286$, $0.59381$, $4.2335$ | $-21.931$, $22.022$ |
| 5 | $-7.1591$, $-2.2215$, $0.15969$, $1.7186$, $8.6833$ | $-22.465$, $24.346$ |

As $H_{2n}$ grows, the other speeds spread out; as it falls to 0, they join the
abscissas. On the boundary, the speeds are these abscissas, each twice, and
their mean (Remark 4.9 of the paper). There, $m_{2n}$ is rounded, so the state
holds either the classical quadrature of $n$ nodes, as in this run, or the
HyQMOM quadrature of $n + 1$ nodes, one of them with almost no weight; the
speeds are the same either way, to rounding.
