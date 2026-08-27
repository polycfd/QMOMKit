# Minimal 1D QMOM phase-space evolution

This is the minimal QMOMKit QMOM API example. It shows the smallest complete
phase-space evolution loop used by an application:

```text
moments
   |
   v
QmomCase::initial_state()
   |
   v
QmomState
   |
   v
advance_qmom_phase_space(...)
   |
   v
QmomState
```

The driver keeps top-level orchestration in `main.cpp`, while the example-local
physics callbacks live in `physics.cpp`. The user provides:

- initial raw moments;
- a QMOM node count;
- time controls;
- phase-space advection physics through `qmomkit::PhaseSpaceVelocity`;
- moment-source physics through `qmomkit::MomentSource`.

Those callbacks are required by this example because its model includes both
advection motion and fragmentation. Other QMOMKit phase-space evolution workflows can use
only the callbacks needed by their physics; see `docs/user/phase_space_evolution.md` for
the lower-level node-advection and moment-source routines.

QMOMKit handles the quadrature representation, moment reconstruction, time
integration, and phase-space evolution splitting inside `advance_qmom_phase_space(...)`.
The accepted `QmomState` carries the canonical moments and quadrature from one
step to the next. Users copy, move, and advance accepted states; they do not
construct them directly from independent moments and quadrature.

The phase-space evolution is

$$
\frac{\mathrm{d}d}{\mathrm{d}t} = 1-d.
$$

The moment-source is binary fragmentation,

$$
d \to \frac{d}{2},
$$

with the same raw-moment source formulation used by the full 1D advection
breakage example.

Requires: QMOMKit.

Required for this example: initial raw moments, node count, velocity callback,
source callback, and an installed QMOMKit package. Optional: change time
controls or phase-space evolution options in `main.cpp`.

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

## What to expect

The binary-fragmentation source increases total number density, while the
advection law moves each abscissa toward `d = 1`. The final moments and
quadrature therefore reflect both the increasing population count and the
relaxing particle-diameter coordinates.
