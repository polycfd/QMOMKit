# QMOM Mixture Initialization Example

This example generates raw moments from an absolute-density distribution
mixture and reconstructs a one-dimensional QMOM quadrature.

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

The executable prints the component densities and reconstructed quadrature.
