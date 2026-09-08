# QMOM Explicit Moments Example

This example reconstructs a one-dimensional QMOM quadrature directly from an
explicit sequence of raw, unnormalized moments.

The phase-space coordinate is $\xi$. For a number density $n(\xi)$, the raw
moment of order $k$ is

$$
m_k = \int_{\mathbb{R}} \xi^k n(\xi)\,\mathrm{d}\xi.
$$

The source uses the manufactured two-node discrete measure

$$
n(\xi) = 0.4\,\delta(\xi+1) + 0.6\,\delta(\xi-2),
$$

which gives $m_0=1$, $m_1=0.8$, $m_2=2.8$, and $m_3=4.4$. QMOMKit is
asked for two nodes, so the expected reconstruction is weights 0.4 and
0.6 at abscissas -1 and 2, up to numerical roundoff.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_explicit_moments_example
```

The executable writes no files. It prints the total weight and reconstructed
node weights and abscissas to the terminal.
