# QMOM Distribution Initialization Example

This example generates raw moments from an analytical log-normal distribution
and reconstructs a one-dimensional QMOM quadrature.

Requires: QMOMKit.

Required: distribution parameters in `main.cpp` and an installed QMOMKit
package. Optional: change the distribution family, number density, or node
count.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_distribution_initialization_example
```

The executable prints the generated raw moments and reconstructed quadrature.
