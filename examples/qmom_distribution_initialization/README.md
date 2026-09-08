# QMOM Distribution Initialization Example

This example generates raw moments from an analytical log-normal distribution
and reconstructs a one-dimensional QMOM quadrature.

The phase-space coordinate is the scalar distribution variable $\xi$. The
implemented distribution is `LogNormalDistribution{-0.2, 0.4}`, where the
parameters are the logarithmic mean and logarithmic standard deviation. The
total number density is `2.5`. QMOMKit requests moments through order
`2 * node_count - 1 = 3` for the requested two-node quadrature.

The generated values are raw, unnormalized moments. In particular, $m_0$ is
the total number density `2.5`; it is not normalized away before QMOM
reconstruction. QMOMKit evaluates the analytical distribution moments, then
`compute_qmom(...)` converts those moments into two weights and two
abscissas. The printed weights and abscissas are the finite QMOM
representation of the log-normal distribution, not samples from it.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_distribution_initialization_example
```

The executable writes no files. It prints the generated raw moments, including
$m_0$, followed by the reconstructed quadrature weights and abscissas.
