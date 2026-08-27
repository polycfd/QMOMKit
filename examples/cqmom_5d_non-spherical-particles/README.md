# 5D CQMOM non-spherical-particle evolution

This example evolves a non-spherical particle number-density distribution with
CQMOM and provides an independent stochastic particle-based Monte Carlo
reference. Both executables load the supplied configuration path and start
from the same deterministic sampled population. The CQMOM driver advances one
`CqmomState` with QMOMKit's native `advance_cqmom_phase_space(...)` routine.

The phase-space coordinate order is `(d, chi, u, v, w)`:

- `d` is particle diameter.
- `chi` is a shape coordinate.
- `u`, `v`, and `w` are velocity components.

The particle callback supplies deterministic phase-space velocity. Binary
breakage removes a parent and creates two equal-volume daughters, relaxes
`chi`, and applies Gaussian recoil to velocity. The CQMOM source and MC event
model use those same rules. Initial raw moments and the MC ensemble come from a
reproducible finite independent-product sample of the configured distributions:
10,000 samples with seed `12345`.

## Configuration

[`configuration.qmom`](configuration.qmom) declares three nodes in each
dimension, adaptive CQMOM branch counts, Strang splitting, RK23 node advection,
and output options. Its support declaration

```text
support_lower = none, 0.0, none, none, none
support_upper = none, 1.0, none, none, none
```

states that `chi` must remain in `[0, 1]`. QMOMKit applies that support during
reconstruction. The example neither clips coordinates nor adds a separate
`chi` check or stopping rule. Configuration diagnostics are disabled with
`diagnostic_level = none`.

## Build and run

Install QMOMKit and its optional Monte Carlo companion first; the companion
provides the deterministic sampled-initialization utility used here. Configure
this standalone consumer with both package prefixes available to CMake:

```sh
cmake -S examples/cqmom_5d_non-spherical-particles -B /tmp/qmomkit-cqmom-example \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<qmomkit-install>;<qmomkit-monte-carlo-install>"
cmake --build /tmp/qmomkit-cqmom-example
/tmp/qmomkit-cqmom-example/cqmom_5d_non_spherical_particles_example \
  examples/cqmom_5d_non-spherical-particles/configuration.qmom
/tmp/qmomkit-cqmom-example/cqmom_5d_non_spherical_particles_mc_example \
  examples/cqmom_5d_non-spherical-particles/configuration.qmom
```

The configuration path is required for both executables. Neither copies nor
searches for configuration files. The CQMOM run prints selected initial and
final moments and the final quadrature. The MC run prints corresponding sampled
reference moments and its final particle count. With the supplied configuration
they write CQMOM output files and `moment_history_mc.txt`, respectively, to the
configured `output` directory relative to the launch directory. The MC driver
does not apply CQMOM reconstruction support as particle clipping.

## Why the support matters

For an intentional failure experiment, make a temporary copy of the supplied
configuration and remove both `chi` support lines. Do not modify the source
configuration. The unbounded evolution can eventually pass `chi > 1`; the
particle model's physical callback then rejects the state through QMOMKit's
normal non-finite-velocity failure path. This demonstrates why the `[0, 1]`
support belongs in the CQMOM configuration.
