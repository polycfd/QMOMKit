# 1D QMOM Aerosol Evaporation

This example, inspired by [Feng et al.
(2023)](https://www.doi.org/10.1016/j.ijmultiphaseflow.2023.104500), evolves a
one-dimensional aerosol diameter distribution through evaporation with a fixed
six-node QMOM quadrature. The phase-space coordinate is the particle diameter,
$d$, in micrometres.

## Model

The initial distribution is a three-component log-normal mixture defined in
[configuration.qmom](./configuration.qmom). QMOMKit generates the raw moments
$M_0,...,M_{11}$ from that analytical definition.

The evaporation law is the $d^2$-law,

$$
\frac{\mathrm{d}(d^2)}{\mathrm{d}t} = -500,
$$

or equivalently

$$
\frac{\mathrm{d}d}{\mathrm{d}t} = -\frac{500}{2d}.
$$

The particle diameter is defined in micrometres and time is defined in seconds. The negative sign moves
the QMOM abscissas toward smaller diameters. A residual diameter of
$`0.8\,\mathrm{\mu m}`$ is retained after each node-advection step.

[main.cpp](./main.cpp) shows the complete QMOMKit workflow:

```text
configuration 
    -> raw moments 
    -> six-node QMOM reconstruction
    -> native node advection 
    -> updated raw moments and output
```

[physics.cpp](./physics.cpp) declares the callback function for the evaporation
law and residual-diameter condition.

The advection of the quadrature nodes preserves all six weights and QMOMKit
reconstructs the raw moments from the moved quadrature.

## Build and run

From the repository root, install QMOMKit first:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cd examples/qmom_1d_aerosol_evaporation
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_1d_aerosol_evaporation_example
```

The executable prints the initial and final total number density and mean
diameter, and the final six-node quadrature. It also writes QMOMKit-native
moment and quadrature histories:

- `output/moments.txt`
- `output/quadrature.txt`

The output directory and cadence are configured in the `[output]` section of
`configuration.qmom`.

## What to expect

This is pure node advection, so the quadrature weights and total particle
number remain unchanged. Diameters move monotonically downward under
evaporation, and the residual-diameter rule keeps every node at or above
$0.8 \ \mu\mathrm{m}$.
