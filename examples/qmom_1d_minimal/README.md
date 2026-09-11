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

The example includes:
- phase-space advection through `qmomkit::PhaseSpaceVelocity`;
- a moment source through `qmomkit::MomentSource`.

QMOMKit handles the quadrature representation, moment reconstruction, time
integration, and phase-space evolution splitting inside `advance_qmom_phase_space(...)`.
The accepted `QmomState` carries the canonical moments and quadrature from one
step to the next. 

The phase-space advection is defined as

$$
\frac{\mathrm{d}d}{\mathrm{d}t} = 1-d.
$$

The moment-source is binary fragmentation,

$$
d \to \frac{d}{2}.
$$

Build and run from this example directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_1d_minimal_example
```

The executable prints a concise console summary with the 
final moments and the final quadrature. It also writes two
QMOMKit-native output files:

- `output/moment_history.txt`
- `output/quadrature_history.txt`

`moment_history.txt` records the raw-moment evolution and `quadrature_history.txt`
records the reconstructed node weights and abscissas. 
Both files are whitespace-delimited text tables written through QMOMKit-native writers.

## What to expect

The binary-fragmentation source increases total number density, while the
advection law moves each abscissa toward $d = 1$. The final moments and
quadrature therefore reflect both the increasing population count and the
relaxing particle-diameter coordinates.
