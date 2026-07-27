# QMOM Explicit Moments Example

This example reconstructs a one-dimensional QMOM quadrature directly from an
explicit raw moment sequence.

Requires: QMOMKit.

Required: raw moments in `main.cpp` and an installed QMOMKit package.
Optional: replace the moment values with moments from your application.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_explicit_moments_example
```

The executable prints the total weight and reconstructed node weights and
abscissas.
