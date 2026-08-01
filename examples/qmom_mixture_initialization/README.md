# QMOM Mixture Initialization Example

This example generates raw moments from an absolute-density distribution
mixture and reconstructs a one-dimensional QMOM quadrature.

The phase-space coordinate is the scalar distribution variable $\xi$. The
mixture contains:

- a gamma component with number density `1.5`, shape `2.0`, and scale `0.4`;
- a log-normal component with number density `0.8`, logarithmic mean `-0.2`,
  and logarithmic standard deviation `0.5`.

The component densities are absolute number densities, not normalized mixture
fractions. Therefore the total number density is
$m_0 = 1.5 + 0.8 = 2.3$. For raw, unnormalized moments, each component's
analytical moment is multiplied by its absolute number density and the
components are added. QMOMKit requests moments through order `5` for the
requested three-node quadrature.

Requires: QMOMKit.

Required: mixture components in `main.cpp` and an installed QMOMKit package.
Optional: use this after `qmom_distribution_initialization` when one
analytical family is not enough.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_mixture_initialization_example
```

The executable writes no files. It prints the component densities and $m_0$,
then prints the reconstructed QMOM weights and abscissas. The quadrature is a
three-node representation of the combined mixture.
