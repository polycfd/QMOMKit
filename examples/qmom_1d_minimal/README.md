# Minimal 1D QMOM transport

This is the minimal QMOMKit QMOM API example. It shows the smallest complete
transport loop used by an application:

```text
moments
   |
   v
advance_transport_qmom(...)
   |
   v
moments
```

The example keeps all application logic in `main.cpp`. The user provides:

- initial raw moments;
- a QMOM node count;
- time controls;
- continuous phase-space physics through `qmomkit::PhaseSpaceVelocity`;
- discontinuous moment physics through `qmomkit::MomentSource`.

Those callbacks are required by this example because its model includes both
continuous motion and fragmentation. Other QMOMKit transport workflows can use
only the callbacks needed by their physics; see `docs/user/transport.md` for
the lower-level node-advection and moment-source routines.

QMOMKit handles the quadrature representation, moment reconstruction, time
integration, and transport splitting inside `advance_transport_qmom(...)`.

The continuous physics is

```text
dd/dt = 1 - d
```

The discontinuous source is binary fragmentation,

```text
d -> d/2
```

with the same raw-moment source formulation used by the full 1D continuous
breakage example.

Requires: QMOMKit.

Required for this example: initial raw moments, node count, velocity callback,
source callback, and an installed QMOMKit package. Optional: change time
controls or transport options in `main.cpp`.

Build and run from this example directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_1d_minimal_example
```

The executable prints a concise console summary with the initial raw moments,
final raw moments, and final quadrature nodes and weights. It also writes two
QMOMKit-native output files:

- `output/moment_history.txt`
- `output/quadrature_history.txt`

`moment_history.txt` records the raw-moment evolution, which is the primary
distribution history for many moment-method workflows. `quadrature_history.txt`
records the reconstructed node weights and abscissas, which are useful for
inspecting the numerical representation. Both files are whitespace-delimited
text tables written through QMOMKit-native writers.
