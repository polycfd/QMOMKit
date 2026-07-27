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
tools/moment_rate_element_operator/build/moment_rate_element_operator_benchmark --repetitions 1000
```

The timing executable prints CSV rows. Categories include:

- spatial interpolation only;
- QMOM continuous local-rate calls only;
- CQMOM continuous local-rate calls only;
- QMOM discontinuous local-rate calls only;
- CQMOM discontinuous local-rate calls only;
- projection and inverse-mass application only;
- complete element RHS;
- complete external SSP-RK3 coefficient step.

Continuous CQMOM local-rate calls perform CQMOM inversion before velocity
evaluation and rate reconstruction. Native QMOMKit node integration avoids that
repeated public local-rate inversion in compatible pointwise workflows, so it
can remain the cheaper path for those formulations. Measure before drawing
optimization conclusions.

The benchmark intentionally does not add workspaces, caches, batching, or
pre-inverted overloads. It is a first measurement scaffold for deciding what,
if anything, deserves optimization later.
