# 2D Minimal CQMOM Example

## Purpose

This is the main QMOMKit example for CQMOM phase-space evolution. It initializes a CQMOM
case through the standard configuration path, creates the runtime state, then
evolves that state with user-supplied advection and moment-source physics
callbacks. The example also includes an optional phase-space boundary exchange
callback to demonstrate nucleation-style moment production.

The supplied configuration is two-dimensional so the output is easy to read:

The phase-space ordering is $(d,u)$,

where `d` is particle size and `u` is velocity. The implementation is written
against coordinate vectors and CQMOM source contexts, so the physics callbacks
do not require `D == 2`.

The workflow is:

```text
configuration.qmom
        |
        v
load_cqmom(...) -> CqmomCase
        |
        v
compute_cqmom(...) -> initial quadrature
        |
        v
CqmomCase::initial_state() -> CqmomState
        |
        v
physics callbacks
        |
        v
advance_cqmom_phase_space(...)
        |
        v
updated state and quadrature
```

`CqmomCase` stores the configuration and problem definition.
`CqmomState` is the evolving runtime state passed to
`advance_cqmom_phase_space(...)`.

The callbacks in this example are one representative combination: advection
velocity, moment source, and optional boundary exchange. QMOMKit does not store
these physical models in the case or state, and other workflows can provide
only the callbacks required by their model.

## Configuration

Edit `configuration.qmom` to change the distribution and CQMOM hierarchy. The
example uses an independent-product initial distribution:

$$
f(d,u) = f_d(d)f_u(u).
$$

The default configuration uses:

$$
d \sim \mathcal{N}(2, 0.5^2),
\qquad
u \sim \mathcal{N}(0, 0.5^2).
$$

For number density $n(d,u)$, the raw, unnormalized CQMOM moments are

$$
M_{\alpha_0,\alpha_1}
= \int\!\int d^{\alpha_0}u^{\alpha_1}n(d,u)\,
\mathrm{d}d\,\mathrm{d}u.
$$

The zeroth moment $M_{0,0}$ is the configured total number density.

QMOMKit reads this configuration with `qmomkit::config::load_cqmom(...)`,
which returns the user-facing CQMOM case consumed by `qmomkit::compute_cqmom(...)`
and `qmomkit::CqmomCase::initial_state()`.

The same file also uses `[time_integration.node_advection]` to select the
advection node-trajectory integrator. The shipped configuration chooses
`scheme = ssp_rk3` for fixed-stage, predictable node motion; `scheme = rk23`
is the default adaptive Bogacki-Shampine 3(2) choice, and `scheme = rk45`
selects Dormand-Prince RK5(4)7M. Embedded-RK controls can be supplied under
`[time_integration.node_advection.adaptive]`; omitted values retain the
`AdaptiveRkOptions` defaults, including `initial_step = 1.0e-6 s` for RK23
and RK45. QMOMKit does not estimate the initial step from the state,
derivatives, or tolerances. The executable prints the effective
node-advection settings at startup. This choice does not change the
Strang/`advsrc`/`srcadv` phase-space evolution ordering or the SSP-RK3 moment-source update.

The default CQMOM hierarchy is `node_counts = 2, 2`, meaning:

- two quadrature nodes in `d`;
- two conditional quadrature nodes in `u`;
- four terminal CQMOM nodes.

Adaptive branch counts are also supported in this normal `CqmomState`
phase-space evolution workflow. When `phase_space.adaptive_branch_counts = true`,
`node_counts` are treated as maximum branch counts and CQMOM automatically uses
adaptive Wheeler for marginal and conditional inversions. The reconstructed
quadrature may contain fewer active terminal nodes. Output should therefore be
read through terminal weights and abscissas rather than assuming a fixed
rectangular node grid. The standalone QMOM adaptive option remains independent.

To adapt the example to more dimensions, update `names`, `node_counts`,
`conditioning_order`, and add one `[initial.dimension.N]` section per
phase-space dimension. The normal workflow still uses `load_cqmom(...)`.

## Advection Velocity

The user-defined velocity is in `physics.cpp`. The default model treats
coordinate `0` as the breakage/growth coordinate and relaxes every later
coordinate toward zero:

$$
\begin{aligned}
\frac{\mathrm{d}x_0}{\mathrm{d}t} &= 1-x_0, \\
\frac{\mathrm{d}x_k}{\mathrm{d}t} &= -x_k, \qquad k>0.
\end{aligned}
$$

For the shipped two-dimensional configuration, this is
$\mathrm{d}d/\mathrm{d}t = 1-d$ and $\mathrm{d}u/\mathrm{d}t = -u$. The phase-space evolution helper advances the terminal CQMOM nodes during
the advection substeps.

## Moment Source

The moment-source is also in `physics.cpp`. It models binary
fragmentation with constant rate along coordinate `0`:

$$
x_0 \to \frac{x_0}{2}.
$$

All other coordinates are unchanged by fragmentation, so the source acts on
each CQMOM moment as:

$$
\frac{\mathrm{d}M_{\boldsymbol{\alpha}}}{\mathrm{d}t}
= \lambda\left(2^{1-\alpha_0}-1\right)M_{\boldsymbol{\alpha}}.
$$

The example also defines a small `CqmomMomentBoundaryContribution` representing
creation of new samples through the boundary of the first phase-space
dimension. Boundary contributions provide additive raw-moment rates from
phase-space exchange; in this example, the injected size is prescribed in the
first coordinate. The callback depends only on the exponent of that injected
coordinate and returns a unit raw-moment factor for every unspecified
coordinate dimension. This is an example-specific assumption; QMOMKit does not
choose or infer a missing injected-coordinate distribution.

For a moment multi-index $\boldsymbol{\alpha}$, the boundary exchange contribution is:

$$
\left.\frac{\mathrm{d}M_{\boldsymbol{\alpha}}}{\mathrm{d}t}\right|_{\mathrm{boundary}}
= J x_b^{\alpha_0}.
$$

where $J$ is the nucleation or injection rate and $x_b$ is the injected size in
the first phase-space dimension. The callback receives the CQMOM moment
context, evaluates the raw-moment rate for each context entry, and returns one
rate value per required moment. Applying the contribution consistently across
the affected moment set lets an application add its phase-space exchange model
without modifying QMOMKit.

This is not a mesh boundary condition. It does not prescribe mesh faces,
normals, or reflection behavior; it is an additive raw-moment rate supplied by
a user-defined physics callback. To run the same example without boundary
exchange, remove the `boundary` callback from `main.cpp` and call the
source-only `advance_cqmom_phase_space(...)` overload.

The example advances advection motion, the source, and the optional boundary
contribution with `qmomkit::advance_cqmom_phase_space(...)` and
`qmomkit::CqmomState`. The operator-splitting scheme remains Strang splitting:

```text
advection half-step
        |
        v
moment-source update, including optional boundary contribution
        |
        v
post-source CQMOM reconstruction
        |
        v
advection half-step
```

The time controls are local constants in `main.cpp`:

- `initial_time = 0.0`
- `final_time = 1.0`
- `time_step = 0.01`

The output directory, enabled files, and output cadence are read from
`[output]` in `configuration.qmom` through `qmomkit::OutputOptions`.

The CQMOM source and boundary contribution callbacks receive a
`qmomkit::CqmomSourceContext`. For entry `i`, `context.index(i)` is the moment
multi-index and `context.value(i)` is the current value of that moment. The
multi-index follows phase-space order, not conditioning order. Return one rate
value per context entry in the same order: `rate[i]` is
$\mathrm{d}M_{\boldsymbol{\alpha}}/\mathrm{d}t$ for `context.index(i)`.

Most users adapting this example change:

- `configuration.qmom`: phase-space names, node counts, conditioning order,
  initial distributions, and `[output]` settings;
- `physics.cpp`: `advection_velocity()`, `binary_fragmentation_source()`, and
  optional boundary exchange callbacks;
- `main.cpp`: time controls and `CqmomPhaseSpaceEvolutionOptions`.

## Output

Requires: QMOMKit.

Required: `configuration.qmom`, `main.cpp`, `physics.cpp`, `utilities.cpp`,
and an installed QMOMKit package. Optional: output-file analysis after the run.

Build from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

The build copies `configuration.qmom` into `build/`. Run from the build
directory so generated output stays there:

```bash
cd build
./cqmom_2d_minimal_example
```

Running the executable writes:

- `output/initial_moments.txt`
- `output/initial_quadrature.txt`
- `output/moment_history.txt`
- `output/quadrature_history.txt`

It also prints a concise console summary with selected initial and final raw
moments and the final terminal quadrature nodes and weights.

`initial_moments.txt` records the raw moments reconstructed from the initial
CQMOM quadrature. `moment_history.txt` records the raw-moment evolution at the
configured output cadence. These moment files are the primary distribution
history for many workflows.

`initial_quadrature.txt` and `quadrature_history.txt` record the terminal node
weights and abscissas. These quadrature files are useful for inspecting the
representation and diagnosing reconstruction behavior. All four files are
written through QMOMKit-native writers as whitespace-delimited text tables with
commented headers.

## What to expect

The terminal nodes relax under the configured phase-space advection law. The
fragmentation source and additive boundary birth both change the raw-moment
history, so the output reflects their combined number-density and
coordinate-moment effects rather than advection alone. Inspect the terminal
quadrature history together with the moment history when adapting the callback
combination to another model.
