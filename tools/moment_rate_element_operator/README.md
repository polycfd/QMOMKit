# Moment-Rate Element Operator Prototype

This directory is a standalone downstream prototype for QMOMKit's public local
moment-rate API. It is a manufactured element operator, not a complete CG or DG
solver.

The prototype demonstrates this external workflow:

```text
element coefficients
-> interpolate local moments at spatial quadrature points
-> call QMOMKit local moment-rate evaluators
-> project rates back to element coefficients
-> apply the inverse element mass matrix
```

All spatial interpolation, projection, and mass-matrix logic lives in this
standalone project. The code consumes only the installed QMOMKit package and
public headers.

## Build

Install QMOMKit first, using a Release build for timings:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build
```

Configure this downstream project against that install prefix:

```sh
cmake -S tools/moment_rate_element_operator \
  -B tools/moment_rate_element_operator/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/install"

cmake --build tools/moment_rate_element_operator/build
```

## Correctness Check

Run:

```sh
tools/moment_rate_element_operator/build/moment_rate_element_operator_check
```

The check covers:

- QMOM continuous local rates with an analytical constant-velocity reference;
- QMOM discontinuous local rates from a deterministic raw-moment source;
- CQMOM continuous local-rate composition through the public API;
- CQMOM discontinuous source-context ordering and rate ordering;
- a local coefficient-level SSP-RK3 demonstration using continuous and
  discontinuous RHS evaluators;
- failure diagnostics for a rejected local continuous rate call.

The CQMOM continuous case validates downstream composition and status handling;
it is not an independent validation of CQMOM inversion mathematics.

## Timing

Run a short smoke benchmark:

```sh
tools/moment_rate_element_operator/build/moment_rate_element_operator_benchmark --short
```

Run a longer local timing pass:

```sh
tools/moment_rate_element_operator/build/moment_rate_element_operator_benchmark
```

The benchmark prints a human-readable table by default. Use CSV when collecting
results for a spreadsheet or plotting script:

```sh
tools/moment_rate_element_operator/build/moment_rate_element_operator_benchmark --csv > moment-rate-element-operator.csv
```

Useful controls include:

```sh
tools/moment_rate_element_operator/build/moment_rate_element_operator_benchmark --samples 9 --target-ms 50
tools/moment_rate_element_operator/build/moment_rate_element_operator_benchmark --filter cqmom_3x2 --csv
```

Each case performs an untimed fixture validation, an untimed warm-up, an
automatic iteration calibration to reach the requested per-sample target
duration, and repeated measured samples. The reported median is the primary
timing value; minimum and maximum sample timings are printed as a simple spread
indicator. Every measured result contributes to a checksum consumed after
timing so the compiler cannot discard the work. Console output is outside the
timed loops.

Categories include:

- spatial interpolation only;
- projection and inverse-mass application only;
- public QMOM inversion only;
- public CQMOM inversion only;
- QMOM continuous local-rate evaluation;
- QMOM discontinuous local-rate evaluation;
- CQMOM continuous local-rate evaluation;
- CQMOM discontinuous local-rate evaluation;
- complete element RHS;
- complete external SSP-RK3 coefficient step.

The matrix includes quadratic and cubic manufactured elements, single-element
and batched serial loops, one- and two-node QMOM cases, identity and
nonidentity CQMOM conditioning order, a larger fixed CQMOM hierarchy, an
adaptive branch-specific CQMOM fixture, no and detailed CQMOM diagnostics, zero
and cheap affine velocities, a modest nonlinear velocity, and nearly free or
modest arithmetic discontinuous callbacks.

The inversion-only rows call the public `compute_qmom(...)` and
`compute_cqmom(...)` entry points with the same local moment values used by the
continuous local-rate rows. They are included to compare full continuous-rate
timing with public inversion timing. Do not treat full-call minus
inversion-only timing as a precise component measurement; the two measurements
have different surrounding work, result shapes, and noise.

The CSV output reports compiler/build metadata, target sample duration, case
metadata, operation category, basis size, spatial quadrature count, transported
moment count, phase-space dimension, node configuration, conditioning-order
label, diagnostic level, callback label, element batch size, local call counts,
sample count, calibrated iterations per sample, median total time, median time
per element, median time per local call, minimum sample time, maximum sample
time, and checksum. CSV timing fields are nanoseconds.

Results are serial, machine-specific, compiler-specific, and build-specific.
Use a Release build and compare trends within a single run before drawing
conclusions across machines.

A useful first profiling case is the largest fixed CQMOM continuous local-rate
case with minimal diagnostics and cheap affine velocity:

```sh
tools/moment_rate_element_operator/build/moment_rate_element_operator_benchmark \
  --filter cqmom_3x2_identity_none_affine_velocity \
  --target-ms 100
```

On macOS, run that command under Instruments Time Profiler. On Linux, use a
similar filtered run under `perf`. The filter isolates repeated public
continuous local-rate evaluation better than a full matrix run and keeps
callback cost modest.

The benchmark intentionally does not add workspaces, caches, batching, or
pre-inverted overloads. It is a first measurement scaffold for deciding what,
if anything, deserves optimization later.
