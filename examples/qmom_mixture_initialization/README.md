# QMOM Mixture Initialization Example

This example loads an absolute-density distribution mixture from
[configuration.qmom](./configuration.qmom), lets QMOMKit generate its raw moments, and reconstructs a
one-dimensional QMOM quadrature.

The phase-space coordinate is the scalar distribution variable $\xi$. The
mixture contains:

- a gamma component with number density $1.5$, shape $2.0$, and scale $0.4$;
- a log-normal component with number density $0.8$, logarithmic mean $-0.2$,
  and logarithmic standard deviation $0.5$.

The component densities are absolute number densities, not normalized mixture
fractions. Therefore the total number density is
$m_0 = 1.5 + 0.8 = 2.3$. For raw, unnormalized moments, each component's
analytical moment is multiplied by its absolute number density and the
components are added. QMOMKit requests moments through order $5$ for the
requested three-node quadrature.

The native configuration syntax is `type = mixture` with contiguous
`[initial.mixture.component.K]` sections. Each component supplies its
`number_density` in [configuration.qmom](./configuration.qmom); QMOMKit performs the mixture moment
generation during `load_qmom(...)`.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_mixture_initialization_example
```

The executable writes no files. It prints $m_0$, then the reconstructed QMOM
weights and abscissas. The quadrature is a three-node representation of the
combined mixture.
