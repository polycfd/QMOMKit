# Choosing An Example

QMOMKit examples are small standalone applications that consume an installed
QMOMKit package. They demonstrate numerical workflows rather than complete
applications or solvers. Start with the stages relevant to your work; the full
suite is not a mandatory linear sequence.

Install QMOMKit before building an example. Examples with Monte Carlo (MC)
companions also require the optional QMOMKitMonteCarlo package; see
[installation.md](installation.md) and [monte_carlo.md](monte_carlo.md).

Examples that write moments or quadrature use whitespace-delimited QMOMKit
`.txt` files with commented headers. Individual READMEs describe their specific
build commands, configuration, output, and expected behavior.

## How To Use This Guide

The progression below starts with 1D QMOM reconstruction, then introduces
phase-space evolution, CQMOM, numerical controls, and specialist workflows.
Each stage states the prerequisite knowledge it assumes and identifies optional
branches.

For a concise repository-level map and generic build pattern, see
[`examples/README.md`](../../examples/README.md).

## Stage 1: Raw Moments And Initialization

Start here before evolving a QMOM or CQMOM state. These examples establish that
QMOMKit consumes raw, unnormalized number-density moments: the zeroth moment is
the total number density and is not silently normalized.

### `qmom_explicit_moments`

Use this first when another application already provides `m0, m1, ...`. It
teaches the direct `MomentSet -> compute_qmom(...)` reconstruction path and
uses a manufactured two-node measure with an exact expected quadrature.

Example README:
[`examples/qmom_explicit_moments/README.md`](../../examples/qmom_explicit_moments/README.md)

### `qmom_distribution_initialization`

Use this next when an analytical distribution, rather than an existing moment
vector, defines the initial condition. It generates log-normal raw moments and
then reconstructs their two-node QMOM representation.

Example README:
[`examples/qmom_distribution_initialization/README.md`](../../examples/qmom_distribution_initialization/README.md)

### Optional: `qmom_mixture_initialization`

Use this when one distribution family is insufficient. It shows that mixture
component `number_density` values are absolute densities, not normalized
fractions, and that the component moments add to form the QMOM input moments.

Example README:
[`examples/qmom_mixture_initialization/README.md`](../../examples/qmom_mixture_initialization/README.md)

### Optional: `realizability_check`

Use this diagnostic branch when you need to classify a raw moment sequence
under a declared support assumption before reconstructing it. It is not an
evolution or inversion workflow.

Example README:
[`examples/realizability_check/README.md`](../../examples/realizability_check/README.md)

For configuration syntax, including explicit moments, distributions, mixtures,
and realizability policy, see [configuration.md](configuration.md). Complete
`.qmom` examples live beside the executable workflows that use them.

## Stage 2: Basic QMOM Evolution

Prerequisite: understand the Stage 1 raw-moment and quadrature concepts. These
examples separate the main 1D evolution mechanisms before combining them.

### `qmom_1d_minimal`

Start here for the smallest complete `QmomCase`/`QmomState` evolution loop. It
uses `advance_qmom_phase_space(...)` with a velocity callback and a
fragmentation source, showing the normal library-managed state workflow.

Example README:
[`examples/qmom_1d_minimal/README.md`](../../examples/qmom_1d_minimal/README.md)

### `qmom_1d_aerosol_evaporation`

Use this next to isolate node advection in a physical setting. A
configuration-defined six-node aerosol distribution evolves by evaporation;
weights remain fixed and no raw-moment source is applied.

Example README:
[`examples/qmom_1d_aerosol_evaporation/README.md`](../../examples/qmom_1d_aerosol_evaporation/README.md)

### `qmom_1d_breakage_advection`

Use this after the minimal workflow when the model needs both node motion and a
raw-moment source. It is the readable normal example for library-managed Strang
splitting: advection acts on nodes and binary fragmentation acts on moments.

Example README:
[`examples/qmom_1d_breakage_advection/README.md`](../../examples/qmom_1d_breakage_advection/README.md)

## Stage 3: Application-Owned QMOM Workflows

Prerequisite: complete a normal library-managed QMOM evolution workflow first.
These are advanced integration patterns for applications that own more of the
time-integration machinery.

### `qmom_1d_explicit_strang_splitting`

Use this when the surrounding application must own the
`C(dt/2) -> D(dt) -> C(dt/2)` sequence. QMOMKit still performs the individual
advection and moment-source stages, state canonicalization, and adaptive
reconstruction; the application decides how to compose successful stages.

Example README:
[`examples/qmom_1d_explicit_strang_splitting/README.md`](../../examples/qmom_1d_explicit_strang_splitting/README.md)

### `qmom_1d_breakage_advection_local-rates_with-mc`

Use this host-solver branch when the application owns the raw-moment state and
its SSP-RK3 integration. Unlike the native breakage example, QMOMKit supplies
instantaneous local advection and source moment rates rather than advancing the
whole state. This is useful when placing QMOMKit operations inside a residual,
finite-element, or finite-volume time-integration workflow.

Example README:
[`examples/qmom_1d_breakage_advection_local-rates_with-mc/README.md`](../../examples/qmom_1d_breakage_advection_local-rates_with-mc/README.md)

## Stage 4: QMOM Monte Carlo Companion/Reference Workflow

Prerequisite: understand the normal 1D breakage/advection workflow. MC examples
draw a finite sampled population and evolve particles directly under the same
physical model. They are stochastic particle-based companion/reference
realizations, not automatic formal validation tests or exact solutions. Initial
sampling can be deliberately shared, but exact pointwise agreement later is not
expected because finite sampling and finite quadrature are different
representations.

### `qmom_1d_breakage_advection_with-mc`

Use this to compare sampled initialization and normal QMOM evolution with a
direct particle companion. It is the natural MC branch from
`qmom_1d_breakage_advection`.

Example README:
[`examples/qmom_1d_breakage_advection_with-mc/README.md`](../../examples/qmom_1d_breakage_advection_with-mc/README.md)

The Stage 3 local-rate example also includes an MC companion; choose it instead
when application-owned moment integration is the point of interest.

## Stage 5: CQMOM Fundamentals

Prerequisite: understand Stage 1 and the distinction between raw moments and a
quadrature. These are the primary CQMOM starting points.

### `cqmom_independent_product_initialization`

Start here for multidimensional CQMOM. It demonstrates
`load_cqmom(...) -> CqmomCase -> compute_cqmom(...)`, an independent-product
initial distribution, and the relationship between phase-space order and
conditioning order. It creates a `CqmomState` but does not evolve it.

Example README:
[`examples/cqmom_independent_product_initialization/README.md`](../../examples/cqmom_independent_product_initialization/README.md)

### `cqmom_2d_minimal`

Use this next for the canonical CQMOM phase-space-evolution workflow. It adds
`CqmomState`, velocity, moment-source, and boundary-contribution callbacks to a
two-dimensional conditional hierarchy. It is the right starting point before
adaptive or physics-specialist CQMOM examples.

Example README:
[`examples/cqmom_2d_minimal/README.md`](../../examples/cqmom_2d_minimal/README.md)

## Stage 6: CQMOM Numerical Controls

Prerequisite: understand the hierarchy and evolution workflow in Stage 5. These
examples demonstrate different numerical capabilities and are not
interchangeable.

### `cqmom_2d_adaptive_node_counts`

Use this to see adaptive branch-local topology in an evolving CQMOM state.
Conditional velocity branches collapse as variance disappears, persist during
advection, and recover by controlled growth when the source pulse restores
moment information.

Example README:
[`examples/cqmom_2d_adaptive_node_counts/README.md`](../../examples/cqmom_2d_adaptive_node_counts/README.md)

### `cqmom_2d_borderline_conditioning`

Use this static diagnostic study when reconstruction becomes nearly degenerate.
It separates conditioning-order effects, solver/Vandermonde policy behavior,
failed-stage interpretation, and successful lower-node adaptive reconstruction.
It is not a phase-space physical application.

Example README:
[`examples/cqmom_2d_borderline_conditioning/README.md`](../../examples/cqmom_2d_borderline_conditioning/README.md)

## Stage 7: CQMOM Sampled Physical Workflows

Prerequisite: Stage 5; Stage 6 is useful when interpreting reconstruction
topology. These are application-motivated CQMOM workflows, each with a distinct
numerical lesson. The MC companions follow the convention stated in Stage 4.

### `cqmom_2d_breakage_advection_with-mc`

Use this for sampled multidimensional initialization followed by CQMOM
advection and fragmentation. It compares the resulting finite CQMOM closure
with a direct particle-based companion using the same sampled initialization.

Example README:
[`examples/cqmom_2d_breakage_advection_with-mc/README.md`](../../examples/cqmom_2d_breakage_advection_with-mc/README.md)

### `cqmom_2d_bubbles_cavitation`

Use this for nonlinear two-dimensional node advection in a cavitation model.
The focus is the CQMOM representation of an evolving bubble population, not a
source term or a complete cavitating-flow solver.

Example README:
[`examples/cqmom_2d_bubbles_cavitation/README.md`](../../examples/cqmom_2d_bubbles_cavitation/README.md)

### `cqmom_3d_droplets_condensation-nucleation`

Use this higher-dimensional example after the basic CQMOM workflow to combine
adaptive reconstruction, condensation-driven node advection, and a finite-time
boundary birth contribution for nucleation. Its direct MC companion begins from
the same sampled population.

Example README:
[`examples/cqmom_3d_droplets_condensation-nucleation/README.md`](../../examples/cqmom_3d_droplets_condensation-nucleation/README.md)

## Stage 8: Specialist High-Dimensional CQMOM Workflows

These optional advanced branches assume comfort with CQMOM evolution,
conditioning, and source callbacks. They demonstrate three distinct closure or
support mechanisms. Each also includes a direct stochastic particle-based
companion/reference realization; apply the Stage 4 comparison convention rather
than treating it as an exact validation oracle.

### `cqmom_4d_particles_breakage_recoil`

Use this when a daughter-state transformation needs information not directly
available in the tracked moment set. The source traverses source-stage terminal
quadrature nodes to evaluate a nonlinear four-dimensional recoil transformation.

Example README:
[`examples/cqmom_4d_particles_breakage_recoil/README.md`](../../examples/cqmom_4d_particles_breakage_recoil/README.md)

### `cqmom_5d_particles_breakage_heat-transfer`

Use this complementary source-closure example when a moment-space formulation
requires fractional orders. It accesses them with
`CqmomSourceContext::generalized_moment(...)` while keeping only the configured
integer raw moments in the evolved state.

Example README:
[`examples/cqmom_5d_particles_breakage_heat-transfer/README.md`](../../examples/cqmom_5d_particles_breakage_heat-transfer/README.md)

### `cqmom_5d_non-spherical-particles`

Use this for bounded-support CQMOM in a five-dimensional adaptive hierarchy. It
combines a multimodal `u`-velocity marginal, an aspect-ratio support interval,
adaptive branch reduction, and paired bounded/unbounded workflows that diagnose
conditional support loss without clipping nodes.

Example README:
[`examples/cqmom_5d_non-spherical-particles/README.md`](../../examples/cqmom_5d_non-spherical-particles/README.md)

## Next Steps

The [phase-space evolution guide](phase_space_evolution.md) defines callback
contracts and splitting behavior. The [configuration guide](configuration.md)
defines `.qmom` syntax and mapping rules. Individual example READMEs remain the
authoritative place for model-specific configuration, commands, expected output,
and limitations.
