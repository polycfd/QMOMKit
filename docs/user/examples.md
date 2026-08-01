# Choosing An Example

QMOMKit examples are small standalone applications that consume an installed
QMOMKit package. Use them as starting points for common workflows, not as
complete applications or solvers.

If you have not built the examples before, start with `installation.md`. Most
examples require only QMOMKit. The Monte Carlo reference examples additionally
require the optional QMOMKitMonteCarlo companion package.

Examples that write QMOM or CQMOM moments and quadrature use QMOMKit-native
writers and `.txt` files. These files are whitespace-delimited text tables
with commented header lines.

## Recommended Onboarding Path

Use these examples in order when you want the shortest path from first QMOM
inversion to CQMOM transport.

```text
Primary onboarding
  QMOM
    qmom_explicit_moments
    qmom_distribution_initialization
    qmom_1d_minimal
    qmom_1d_continuous_breakage
  CQMOM
    cqmom_independent_product_initialization
    cqmom_transport

Reference workflows
  qmom_mixture_initialization
  realizability_check
  qmom_1d_continuous_breakage_mc_initialization
  cqmom_2d_continuous_breakage_mc_initialization
  cqmom_condensing_droplets
  cqmom_non-spherical-particles
```

For physics callbacks, start with `qmom_1d_minimal`,
`qmom_1d_continuous_breakage`, or `cqmom_transport`. The callback contract is
explained in `transport.md`.

### Beginner QMOM Workflow

1. `qmom_explicit_moments`

   Start here when another calculation, file, or application already provides
   raw moments. This is the direct `MomentSet -> compute_qmom(...)` path.

2. `qmom_distribution_initialization`

   Use this next when QMOMKit should generate raw moments from an analytical
   distribution before calling `compute_qmom(...)`.

3. `qmom_1d_minimal`

   This is the smallest complete QMOM transport loop using `QmomState` and
   `advance_transport_qmom(...)` with user-supplied velocity and source
   callbacks.

4. `qmom_1d_continuous_breakage`

   Use this when you want the full readable version with configuration,
   explicit Strang stages, physics callbacks, and file output.

### Beginner CQMOM Workflow

5. `cqmom_independent_product_initialization`

   Start here for multidimensional CQMOM. It shows
   `load_cqmom(...) -> CqmomCase -> compute_cqmom(...) -> CqmomQuadrature`.
   When transport is needed, the same case can then create `CqmomState`.

6. `cqmom_transport`

   Use this canonical CQMOM transport example when your workflow needs custom
   continuous velocity, discontinuous moment sources, transport, and
   quadrature output.

## Initialization Workflows

Use these examples when your first question is how to produce the raw moments
that QMOM or CQMOM will consume.

### `qmom_explicit_moments`

Purpose: reconstruct QMOM directly from an explicit raw moment sequence.

Use this example when another code, file, or calculation already provides
`m0, m1, ...` and you want QMOMKit to compute quadrature nodes and weights.

Example README:
[`examples/qmom_explicit_moments/README.md`](../../examples/qmom_explicit_moments/README.md)

### `qmom_distribution_initialization`

Purpose: generate raw moments from an analytical log-normal distribution, then
compute a one-dimensional QMOM quadrature.

Use this example when your initial condition is naturally described by a
supported analytical distribution.

Example README:
[`examples/qmom_distribution_initialization/README.md`](../../examples/qmom_distribution_initialization/README.md)

### `qmom_mixture_initialization`

Purpose: generate raw moments from an absolute-density mixture of analytical
distributions, then compute QMOM.

Use this example when one distribution family is not enough, such as a
population with multiple peaks or regimes.

Example README:
[`examples/qmom_mixture_initialization/README.md`](../../examples/qmom_mixture_initialization/README.md)

## CQMOM Initialization Workflow

Use this example when you want the normal multidimensional CQMOM setup path.

### `cqmom_independent_product_initialization`

Purpose: load a multidimensional CQMOM case from configuration and compute its
initial quadrature. If the workflow later enters transport, the same case can
create the initial CQMOM state.

Use this example when you want the normal CQMOM onboarding path without
choosing a storage representation.

Example README:
[`examples/cqmom_independent_product_initialization/README.md`](../../examples/cqmom_independent_product_initialization/README.md)

## Configuration And Diagnostics

Use these examples when you want to load input data from `.qmom` files or
check moment validity before using a state.

### Configuration snippets

Purpose: provide small `.qmom` files for explicit moments, mixtures, and
runtime realizability policy settings.

Use these files when you want to see the configuration syntax without a full
application around it.

Configuration directory:
[`examples/configuration/`](../../examples/configuration/)

### `realizability_check`

Purpose: run an explicit one-dimensional realizability check and print the
reported status.

Use this example when you want to check whether a raw moment sequence is
consistent with a physical support such as the nonnegative real line.

Example README:
[`examples/realizability_check/README.md`](../../examples/realizability_check/README.md)

## Transport And Evolution

Use these examples when you want to evolve a quadrature-based representation
through continuous node motion and moment source terms.

### `qmom_1d_minimal`

Purpose: advance one-dimensional accepted QMOM state through the high-level
`advance_transport_qmom(...)` helper.

Use this example when you want the smallest complete QMOM transport loop with
user-supplied velocity and source callbacks.

Example README:
[`examples/qmom_1d_minimal/README.md`](../../examples/qmom_1d_minimal/README.md)

### `qmom_1d_explicit_strang_splitting`

Purpose: advance an adaptive one-dimensional QMOM state with application-owned
Strang splitting by calling the public continuous and discontinuous stage APIs.

Use this example when the surrounding application needs to own the
`C(dt/2) -> D(dt) -> C(dt/2)` sequence while QMOMKit still owns state
canonicalization and adaptive post-source reconstruction.

Example README:
[`examples/qmom_1d_explicit_strang_splitting/README.md`](../../examples/qmom_1d_explicit_strang_splitting/README.md)

### `qmom_1d_continuous_breakage`

Purpose: initialize and evolve a one-dimensional QMOM representation with
continuous node motion and a discontinuous fragmentation source.

Use this example when you want a readable QMOM transport workflow that keeps
the Strang sequence visible: continuous half-step, moment source update, QMOM
reconstruction, continuous half-step.

Example README:
[`examples/qmom_1d_continuous_breakage/README.md`](../../examples/qmom_1d_continuous_breakage/README.md)

### `cqmom_transport`

Purpose: initialize and evolve a CQMOM state with continuous phase-space
velocity and a discontinuous fragmentation source. The shipped configuration
is two-dimensional, while the physics callbacks are written in a
dimension-ready style.

Use this example when your workflow needs multidimensional moments, a CQMOM
hierarchy, and terminal quadrature output over time.

Example README:
[`examples/cqmom_transport/README.md`](../../examples/cqmom_transport/README.md)

## Monte Carlo Reference

These examples require the optional QMOMKitMonteCarlo package. They are
sampled-initialization and reference workflows, not required for normal QMOM
or CQMOM use.

Use these examples when you want one of the Monte Carlo companion tool's two
roles:

- sample realization: draw a finite sample ensemble from a distribution and
  evolve or inspect that ensemble directly;
- independent reference generation: generate a stochastic Monte Carlo
  raw-moment history that can later be compared against a deterministic
  QMOMKit moment history.

The sampled ensemble is finite, so its moments include statistical sampling
uncertainty. Increasing the sample count generally reduces that uncertainty.

### `qmom_1d_continuous_breakage_mc_initialization`

Purpose: initialize a one-dimensional QMOM workflow from a fixed-seed Monte
Carlo sample and provide a companion direct Monte Carlo evolution executable.

Use this example when you want to compare analytical initialization with
sample-derived initialization, or when you want a one-dimensional sampled
sample realization and reference workflow.

Example README:
[`examples/qmom_1d_continuous_breakage_mc_initialization/README.md`](../../examples/qmom_1d_continuous_breakage_mc_initialization/README.md)

### `cqmom_2d_continuous_breakage_mc_initialization`

Purpose: initialize a two-dimensional CQMOM workflow from a fixed-seed Monte
Carlo sample and provide a companion direct Monte Carlo evolution executable.

Use this example when you want a multidimensional sampled-initialization
workflow with shared phase-space ordering, velocity callbacks, output
conventions, and a direct sampled-ensemble reference path.

Example README:
[`examples/cqmom_2d_continuous_breakage_mc_initialization/README.md`](../../examples/cqmom_2d_continuous_breakage_mc_initialization/README.md)

### `cqmom_condensing_droplets`

Purpose: initialize and evolve a three-dimensional CQMOM representation of a
condensing droplet population with a companion direct Monte Carlo evolution
executable.

Directory: `examples/cqmom_condensing_droplets`

Phase space: `(r,w,T)`, where `r` is radius, `w` is vertical velocity, and `T`
is droplet temperature.

Executables:

- `cqmom_condensing_droplets_example`
- `cqmom_condensing_droplets_mc_example`

Use this example when you want a 3D sampled-initialization workflow where
CQMOM and Monte Carlo start from the same finite sample realization and evolve
under the same continuous-only droplet phase-space callback. Both paths use
RK23; the CQMOM driver calls `advance_cqmom_continuous(...)`, and the Monte
Carlo driver calls `qmomkit::monte_carlo::evolve_ensemble(...)`.

Build and run it with the normal standalone example pattern:

```bash
cd examples/cqmom_condensing_droplets
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_condensing_droplets_example
./cqmom_condensing_droplets_mc_example
```

The default run writes `output/initial_moments.txt`,
`output/moment_history.txt`, `output/quadrature_history.txt`, and
`output/moment_history_mc.txt`. The CQMOM and Monte Carlo moment histories use
the same case-defined moment columns, labels, ordering, and output times.

Example README:
[`examples/cqmom_condensing_droplets/README.md`](../../examples/cqmom_condensing_droplets/README.md)

For broader Monte Carlo guidance, see `monte_carlo.md`.

## Suggested Path

If you are new to QMOMKit, a practical order is:

1. `qmom_explicit_moments`
2. `qmom_distribution_initialization`
3. `qmom_1d_minimal`
4. `qmom_1d_continuous_breakage`
5. `cqmom_independent_product_initialization`
6. `cqmom_transport`
7. Optional initialization examples only when needed:
   `qmom_mixture_initialization`.
8. Reference examples only when needed: `realizability_check` and Monte Carlo
   initialization examples.
