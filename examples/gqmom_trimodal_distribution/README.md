# GQMOM of a trimodal distribution

This example reconstructs the trimodal number density of Section 4.1 of
[Fox, Laurent & Passalacqua (2023)](https://doi.org/10.1016/j.jaerosci.2022.106096)
with GQMOM, the `generalized` QMOM family, and with classical QMOM, the case of
Fig. 5 of the paper. It loads the case with `load_qmom(...)` from
[configuration.qmom](./configuration.qmom) and reconstructs both quadratures
with `compute_qmom(...)`.

## Model

The number density $n(\xi)$ of a dimensionless coordinate $\xi$ is the sum of
three normal distributions of equal number density $1/3$, with means $-3$, $0$
and $4$ and standard deviations $1/10$, $1/20$ and $1/30$, so $m_0 = 1$. With
node count $N$ (`node_counts`), QMOMKit generates $m_0, \ldots, m_{2N}$.

- **GQMOM** reconstructs $N_G = 101$ nodes (`generalized_node_counts`) that
  reproduce $m_0, \ldots, m_{2N}$, with the default `normal` continuation and
  $\nu = 1$.
- **Classical QMOM** of the same case, with `family = classical`, reconstructs
  $N$ nodes that reproduce $m_0, \ldots, m_{2N-1}$.

With $N = 3$, two of the GQMOM modes remain nearly points, and the middle one
becomes a continuous distribution of many nodes. With $N = 6$ and $9$, the
GQMOM nodes cluster around all three modes.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./gqmom_trimodal_distribution_example
```

Set `node_counts` in `build/configuration.qmom` to 6 or 9 and run the example
again for the other node counts of Fig. 5.

## Output

The example prints $m_0$, the number of GQMOM nodes with their base node
count and range, and the classical nodes. With `write_quadrature = true` in
`[output]`, it writes `build/output/gqmom_quadrature.txt` and
`build/output/qmom_quadrature.txt`, one row per node with its weight and
abscissa `x0`. Plotting the weights against the abscissas gives a panel of
Fig. 5, for example with `numpy` and `matplotlib`, from `build`:

```bash
python3 -c "
import numpy as np, matplotlib.pyplot as plt
for name in ('gqmom', 'qmom'):
    d = np.loadtxt(f'output/{name}_quadrature.txt', ndmin=2)
    plt.plot(d[:, 3], d[:, 2], 'o-', markersize=3, label=name)
plt.xlabel('xi'); plt.ylabel('w'); plt.legend()
plt.savefig('quadrature.png')"
```

## What to expect

For $N = 3$, the example prints:

```text
m0 = 1
GQMOM: 101 nodes on a base of 3, from xi = -3.01284 to 4.00192
QMOM: 3 nodes
  xi = -3.00896, weight = 0.331611
  xi = -0.008205, weight = 0.334914
  xi = 3.99996, weight = 0.333475
```

The first and the last GQMOM node carry weights of $0.330$ and $0.333$, at
$\xi = -3.013$ and $4.002$. The other 99 nodes carry the middle mode, a smooth
distribution of weights that peaks near $\xi = -0.12$ and falls below
$10^{-6}$ outside $-0.5 < \xi < 1.1$. GQMOM and classical QMOM share
$m_0, \ldots, m_5$, and GQMOM also reproduces $m_6$.
