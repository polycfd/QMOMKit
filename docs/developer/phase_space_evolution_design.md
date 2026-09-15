# Phase-space evolution design

This page documents the current QMOMKit design for phase-space evolution and
time integration.

For raw-moment, tensor, phase-space, node, abscissa, and weight conventions,
see [numerical_conventions.md](numerical_conventions.md). For QMOM and CQMOM
reconstruction behavior, see [qmom_wheeler.md](qmom_wheeler.md) and
[cqmom_recursive.md](cqmom_recursive.md).

## Contents

- [Phase-Space Evolution Architecture Overview](#phase-space-evolution-architecture-overview)
- [Node Advection](#node-advection)
- [Moment Sources and Boundary Contributions](#moment-sources-and-boundary-contributions)
- [Time Integration](#time-integration)
- [Time Integration Details](#time-integration-details)
- [Phase-Space Evolution Composition](#phase-space-evolution-composition)
- [Failure And Status Handling](#failure-and-status-handling)
- [Public API Boundary](#public-api-boundary)
- [Batch Convenience Wrappers](#batch-convenience-wrappers)
- [Validation Expectations](#validation-expectations)

## Phase-space evolution architecture overview

QMOMKit phase-space evolution is built from three numerical operations:

- node advection: phase-space advection of quadrature abscissas in phase space;
- moment-source evolution: direct evolution of raw moments by a caller-supplied
  source term;
- operator-split phase-space evolution: composition of reconstruction, node advection,
  moment reconstruction, moment-source evolution, and post-source inversion.

The library provides the numerical utilities. Users provide physical models:

- phase-space velocity callbacks;
- moment-source callbacks;
- boundary contribution callbacks;
- outer timestep loops;
- physical interpretation of coordinates and moment components.

QMOMKit does not own meshes, cells, finite-volume or finite-element operators,
boundary conditions, MPI communication, particle models, or application
physics. Phase-space evolution operates on accepted QMOM/CQMOM state,
quadrature nodes, raw moments, and callback interfaces only.

The callback categories are independent runtime inputs. A workflow may use only
phase-space advection, only moment-source evolution, moment source plus
boundary exchange, velocity plus moment source, or all three categories. Case
and State objects store numerical configuration and evolving moment data; they
do not store physical models or callback objects.

The production composition has two numerical stages: advection
motion is represented by advancing quadrature nodes, while moment-source or
source-like effects are represented by advancing raw moments. Physical-space
transport of moments is not part of this design.

## Node advection

Node advection advances quadrature abscissas through a phase-space advection
velocity field. It preserves quadrature weights and terminal paths.

The public velocity callback is:

```cpp
using PhaseSpaceVelocity =
    std::function<void(double, const std::vector<double>&, std::vector<double>&, std::size_t)>;
```

The callback receives:

- current time;
- all node coordinates, packed node-major;
- the output velocity vector, sized to match and in the same layout;
- the phase-space dimension count.

It fills $\mathrm{d}\xi/\mathrm{d}t$ for every node. QMOMKit validates the
filled size and finiteness, but it does not assign physical meaning to any
coordinate. A non-finite entry is reported with the evaluation time, zero-based
phase-space dimension, offending value, the coordinate vector of the offending
node, and that node's complete velocity vector, so diagnostics stay per node
even though evaluation is packed.

Node-major packing is the same layout the Runge-Kutta integrator already uses
for its state vector. The callback therefore needs no per-node gather or
scatter and writes derivatives straight into the right-hand-side vector.

For callers whose physics is naturally per node, the convenience shape is:

```cpp
using PerNodePhaseSpaceVelocity =
    std::function<std::vector<double>(double, const std::vector<double>&)>;
```

`make_phase_space_velocity` adapts one into a `PhaseSpaceVelocity`, evaluating
it once per node in ascending node order and preserving its per-node
validation. This is the only per-node path; the library itself carries the
packed form everywhere, so there is exactly one integration path and no
per-node overload set to keep in sync.

The adapted form costs one call and one returned vector per node per
right-hand-side evaluation, and the terminal node count is the product of the
per-dimension node counts. Measured over `advance_nodes` with two nodes per
dimension, allocations per call scale as roughly 13 per terminal node for the
adapted form, while a directly written packed callback does not scale with the
node count at all:

| Case | Terminal nodes | Adapted per-node | Packed |
| --- | --- | --- | --- |
| 2D 2^2 | 4 | 80 | 28 |
| 3D 2^3 | 8 | 133 | 29 |
| 4D 2^4 | 16 | 238 | 30 |
| 5D 2^5 | 32 | 447 | 31 |

The residual packed growth is one allocation per phase-space dimension, from
the terminal conditioning-path table copied into the returned quadrature, not
from the velocity callback.

The evolution entry points carry the packed callback all the way
through, so an adapter, when one is used at all, is created once by the caller
rather than per stage or per state. Batch-index-aware wrappers bind the batch
index once per state around the same packed signature. Measured over a full
Strang phase-space step with a zero moment source:

| Case | Terminal nodes | Adapted per-node | Packed | Saved |
| --- | --- | --- | --- | --- |
| 2D 2^2 | 4 | 376 | 272 | 27.7% |
| 3D 2^3 | 8 | 621 | 413 | 33.5% |
| 4D 2^4 | 16 | 1096 | 680 | 38.0% |
| 5D 2^5 | 32 | 2029 | 1197 | 41.0% |

The saved fraction grows with phase-space dimension because the eliminated term
scales with the terminal node count while the rest of the step does not.

For one-dimensional QMOM, node advection advances the one-dimensional
abscissa array. For multidimensional CQMOM, node advection advances terminal
node coordinates stored by phase-space dimension:

```text
abscissas[dimension][node]
```

Both fixed rectangular CQMOM quadratures and adaptive ragged CQMOM quadratures
are represented through terminal-node arrays during advection. For adaptive
ragged quadratures, terminal paths are carried through unchanged.

Internally, node coordinates are packed into a vector state and integrated by
the selected `NodeAdvectionScheme`. The default is adaptive RK23; adaptive
RK45 and fixed-stage SSP-RK3 are also available. SSP-RK3 advances exactly over
the supplied node-advection timestep. After integration, the state is unpacked
back into the same terminal-node representation. Weights are not integrated by
node advection.

After a successful advection step, QMOMKit reconstructs raw moments
from the moved quadrature:
$$m_\alpha = \sum_j W_j \prod_r \xi_{r,j}^{\alpha_r}$$


This reconstruction is the handoff from phase-space advection node motion back to the
moment representation used by moment-source evolution and later inversions.

For adaptive CQMOM, the phase-space hierarchy node counts remain the
configured maximum branch capacity throughout. The current active
branch counts are separate runtime state. A newly accepted adaptive inversion
quadrature provides the terminal paths from which the next state's counts are
derived and stored on `CqmomState` beside the current moments. Node advection
preserves terminal paths and therefore propagates the state's
already accepted active-count metadata without re-deriving it. Moment
reconstruction updates values only; it must not replace configured maximum
counts with active counts or discard the active-count metadata.

## Moment sources and boundary contributions

Moment-source evolution advances raw moments directly:

```text
dm/dt = source(t, m)
```

The public source callback is:

```cpp
using MomentSource =
    std::function<std::vector<double>(double, const QmomSourceContext&)>;
```

The callback receives the current time and a non-owning view of the current
raw-moment SSP-RK stage. It must return a derivative vector with the same
size.

Returning the derivative rather than filling a supplied buffer is a deliberate
asymmetry: `MomentSource` is the one public callback that still returns, while
`PhaseSpaceVelocity` and `OdeRhs` both fill. The callbacks are invoked at
different rates. A moment source runs once per SSP-RK stage, three times per
phase-space step, independent of node count and phase-space dimension. The
velocity callback is evaluated for every quadrature node, and the terminal node
count is the product of the per-dimension node counts, so a per-call allocation
there grew with dimension and dominated the step; `OdeRhs` is the boundary the
velocity is carried across, and is evaluated at the same rate. Measured on a
QMOM step with three nodes, the three returned source vectors account for about
three percent of the step's allocations. The fill form is therefore worth its
extra awkwardness on the per-node path and is not worth it here, where
returning a vector keeps the most common callback the simplest to write.

The library absorbs that asymmetry in one place, where a `MomentSource` is
wrapped as an `OdeRhs`: the returned rates are copied into the supplied buffer
rather than moved, so the integrator keeps the buffer it reuses.

`context.moment(order)` returns a tracked stage value directly; during QMOM
phase-space source evolution, an untracked order is evaluated from a quadrature
reconstructed lazily from that same stage. The context-local cache is never
populated from an accepted quadrature from another stage. Generic source-only
operations without QMOM reconstruction metadata reject untracked requests.
`context.generalized_moment(order)` evaluates a finite, non-negative scalar
order from that same context-local stage quadrature through the shared
generalized-moment evaluator. QMOMKit validates the source result through the
generic ODE solver path, but it does not know what individual moment components
mean.

Phase-space boundary exchange uses a separate public concept:

```cpp
struct MomentBoundaryRate
{
    std::vector<double> values;
};

using MomentBoundaryContribution =
    std::function<MomentBoundaryRate(double, const std::vector<double>&)>;
```

`MomentBoundaryContribution` represents additive raw-moment rates from
nucleation, injection, removal, evaporation/condensation exchange, or similar
phase-space boundary physics. It is not a PDE boundary-condition framework and
does not introduce boundary geometry, node clipping, reflection, mesh faces, or
finite-volume concepts.

For QMOM, the vector is the one-dimensional raw moment sequence. CQMOM uses the
CQMOM-native source callback over a `CqmomSourceContext` view of the current
moment-source moments, so physical CQMOM source models do not need to depend on
dense tensor or lean required-moment storage. Equivalent dense-backed and
lean-backed `CqmomState` objects expose the same selected CQMOM moment layout;
dense-only tensor entries remain internal unless they are part of that layout.

The CQMOM source context exposes:

- the phase-space hierarchy;
- the SSP-RK stage time;
- the number of source moment entries;
- the raw moment multi-index for each entry;
- the raw moment value for each entry.

For source callbacks invoked through CQMOM evolution,
`context.moment(alpha)` returns a tracked entry directly or lazily reconstructs
the requested untracked raw moment from a quadrature built from the current
SSP-RK stage. The reconstructed quadrature is context-local and is never taken
from an accepted state at another RK stage. For quadrature-level source
closures, `context.for_each_terminal_node(...)` visits the global weight and
phase-space-order coordinates of each terminal node from that same cached
quadrature. The coordinate reference is read-only and valid only during the
visitor callback. Source callbacks do not reconstruct CQMOM quadrature
manually for either access pattern.
`CqmomSourceContext::generalized_moment(orders)` provides the corresponding
vector-order access, evaluating a finite, non-negative order vector from that
same current SSP-RK stage reconstruction.

Generic source-only CQMOM contexts have no source-stage reconstruction
metadata. They reject untracked `context.moment(alpha)` and terminal-node
visitor requests rather than selecting an independent reconstruction policy.

The source callback returns one derivative value per context entry in the same
order. CQMOM boundary exchange uses the same representation-neutral context:

```cpp
using CqmomMomentBoundaryContribution =
    std::function<CqmomMomentRate(double, const CqmomSourceContext&)>;
```

Neither CQMOM callback receives `CqmomMomentData`, `CqmomMomentLayout`,
`RawMomentTensor`, recursive branch-local storage, or solver internals.

Physics remains caller-defined. Internal or volumetric models enter through the
source callback. Phase-space boundary exchange enters through the boundary
contribution callback. The evolution code advances the numerical state and
reports failures; it does not define model equations.

Callbacks are combined only when the selected public routine receives them.
Lower-level workflows can call node advection or moment-source evolution
separately. The composed helpers combine the advection and source-like parts
requested by their overloads; absent source-like physics is represented
explicitly by a zero-rate source callback where the helper signature requires
one.

Moment-source updates use SSP-RK3. Optional acceptability predicates in
`SspRk3Options` allow callers to reject candidate moment states, including
states that violate externally defined realizability policies.

When both a source and a boundary contribution are supplied, QMOMKit builds a
single effective moment-source rate for each SSP-RK3 stage:

```text
effective_rate(t, m) = source(t, m) + boundary(t, m)
```

This preserves the existing source-update retry logic, acceptance predicates,
diagnostics, and post-source inversion behavior. Boundary callback exceptions
are wrapped before they enter the existing source failure pathway so
diagnostics identify the boundary contribution callback as the origin. Source
callback exceptions keep the existing source failure behavior.

## Time integration

QMOMKit exposes generic Runge-Kutta utilities, but evolution uses them in
different numerical roles.

### Node-advection integration

Node advection advances quadrature coordinates through:

```text
dxi/dt = velocity(t, xi)
```

The public `NodeAdvectionScheme` selects the trajectory integrator:

- `rk23`: the adaptive Bogacki-Shampine 3(2) scheme [[Bogacki & Shampine, 1989](https://www.doi.org/10.1016/0893-9659(89)90079-7)];
- `rk45`: the adaptive Dormand-Prince RK5(4)7M scheme [[Dormand & Prince, 1980](https://www.doi.org/10.1016/0771-050X(80)90013-3)];
- `ssp_rk3`: the fixed-stage SSP-RK3 scheme [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)].

The choice is a user-facing trajectory-integration policy, not a change to the
composition. SSP-RK3 is appropriate when the surrounding application
controls the timestep and needs fixed callback cost. The adaptive embedded
methods are appropriate when trajectory accuracy inside the requested
step should be controlled by internal error estimation. Neither scheme is
universally preferred; the useful choice depends on the application timestep
policy and velocity-field behavior.

The SSP-RK3 node path is implemented inside node advection, not through the
moment-source update wrapper. It does not use source acceptability predicates,
realizability checks, source retry semantics, or source-specific performance
counters. It performs one three-stage step over the supplied timestep and
reports scheme-neutral node-advection diagnostics.

RK23 and RK45 may take internal accepted and rejected trial steps to reach the
requested end time. Both use an initial step that defaults to `1.0e-6 s` and
may be overridden with a finite, strictly positive value; neither estimates an
initial step automatically.

Configuration-driven drivers load the full `NodeAdvectionOptions` object from
`[time_integration.node_advection]` and the nested
`[time_integration.node_advection.adaptive]` section. Startup reporting should
use the shared node-advection formatter so fixed-step runs report only the
active scheme, while embedded-RK runs also report effective tolerances, the
supplied initial-step value, and the maximum internal accepted-step count.

Node advection does not perform realizability checks. It moves already
reconstructed quadrature coordinates and preserves weights. Realizability
checks apply to raw moment states during moment-source SSP-RK3 evolution.
Quadrature reconstruction and inversion are separate numerical operations.

For adaptive CQMOM evolution through `CqmomState`, active branch counts are
part of accepted state metadata. `CqmomCase::initial_state()` seeds that
metadata from the initial trusted inversion and fails state construction if the
initial adaptive topology cannot be reconstructed. Later state inversions
validate supplied active-count metadata strictly: there must be exactly one
reachable entry per active branch, configured node-count metadata must match
the hierarchy, active counts must be within configured capacity, and stale or
unreachable descendant entries are rejected. During inversion after a
moment-source stage, if an ancestor branch changes active count, descendant
limits below the new branch are invalidated top-down and recomputed from
configured capacity so numeric path reuse does not apply stale branch limits
after a collapse, merge, or growth. Diagnostics report only true node-count
transitions and distinguish reduced, grown, created, and removed branch
records.

`CqmomPhaseSpaceEvolutionOptions` owns the stage policies:
node-advection/embedded-RK controls, moment-source SSP-RK3 controls, and stage
sequencing. Conditional moment controls, standalone QMOM controls, CQMOM
Wheeler controls, adaptive branch-count enablement, CQMOM diagnostic level, and
Vandermonde solver selection are case-owned policies stored in `CqmomState`.
Dense-backed and lean-backed states must pass these policies through to the
same recursive CQMOM, node-advection, moment-reconstruction, and source-update
algorithms; accepting a configured case policy and replacing it with `false`,
`detailed`, `automatic`, or another aggregate default is a bug.

### Moment-source SSP-RK3

The moment-source SSP-RK3 advances raw moment vectors through moment-source
updates. It is an explicit third-order strong-stability-preserving Runge-Kutta
step [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)]
over a generic vector state.

Phase-space evolution uses this path for the moment-source or source-like stage
because that stage acts directly on moments:

```text
m_new = SSP-RK3(m_old, source)
```

The solver supports:

- optional error estimation;
- final-state acceptance predicates;
- stage acceptance predicates;
- retry by timestep halving;
- minimum timestep enforcement.

The solver itself is generic. It does not know about QMOM, CQMOM, moment
realizability, particles, or phase-space coordinates. In phase-space evolution,
the moment-source wrapper supplies configured realizability checks through the
solver's acceptance interfaces [[Nguyen et al.,
2016](https://www.doi.org/10.1016/j.jcp.2016.08.017)]. Later QMOM or CQMOM
reconstruction is a separate numerical operation, not realizability checking.
Realizability checks do not validate inversion, node counts, CQMOM branch
topology, or reconstructed quadrature states.

## Time integration details

The public Runge-Kutta utilities are generic vector ODE solvers. Phase-space
evolution uses them through node-advection and moment-source wrappers where
appropriate, but the solvers themselves do not know about QMOM, CQMOM,
realizability, particles, phase-space coordinates, or physical model semantics.

Both solvers use the callback convention:

```cpp
using OdeRhs =
    std::function<void(double, const std::vector<double>&, std::vector<double>&)>;
```

The callback receives `(t, state)` and fills a derivative buffer the solver
supplies with the same size as `state`. Resizing that buffer, leaving an entry
unwritten, or writing a non-finite value is an error; the unwritten-entry check
is a debug-build check, as it is for the velocity callback. Callback exceptions
are not converted into solver-specific status values; they follow the low-level
solver behavior and propagate out of the generic RK utility.

Filling a supplied buffer rather than returning one makes the allocation cost
of an integration independent of the number of steps it takes: the stage
vectors are allocated once for the whole integration and every stage writes
through them.

Phase-space evolution-facing wrappers add domain-specific failure handling.
Node advection catches phase-space velocity failures and reports
`NodeAdvectionStatus::velocity_failed`; trajectory-integrator failures report
`NodeAdvectionStatus::integration_failed`, with optional embedded-RK status
retained when an adaptive method is selected. Moment-source updates catch
source callback or SSP-RK3 validation failures and report
`MomentSourceUpdateStatus::source_failed`. The composed helpers then map failed
node advection to `PhaseSpaceEvolutionStatus::advection_failed` and failed
moment-source updates to `PhaseSpaceEvolutionStatus::moment_source_failed`.

### SSP-RK3 details

SSP-RK3 advances one raw-moment state over a requested timestep:
$$y' = f(t, y)$$

Conceptually, the step evaluates three explicit stages:
```text
stage 1: forward Euler candidate from y_n
stage 2: convex combination of y_n and a second Euler-like update
final:   convex combination producing the third-order state y_(n+1)
```

The convex-combination structure is the reason SSP-RK3 is used for
moment-source evolution: source updates can be paired with caller-supplied
acceptability predicates while preserving a simple explicit update structure.
QMOMKit uses this path for moment-source or source-like moment evolution, where
the state is a raw-moment vector rather than node coordinates.

The fixed SSP-RK3 entry point takes explicit options:

```cpp
SspRk3Options options;
options.estimate_error = true;
const SspRk3Result result = ssp_rk3_step(state, time, time_step, rhs, options);
```

`SspRk3Options` controls optional solver behavior:

- `estimate_error` requests an embedded lower-order error estimate;
- `accept` checks the completed candidate state;
- `stage_accept` checks stage-1, stage-2, and final candidate states;
- `max_retries` limits retry attempts after predicate rejection;
- `min_time_step` prevents retrying below a configured timestep floor.

If error estimation is enabled, the solver compares the final third-order state
with an embedded lower-order state. The result is a relative Euclidean norm
when the final state norm is nonzero, and an absolute difference norm when the
final state norm is zero. This estimate is diagnostic and does not by itself
accept or reject a step.

Acceptance and rejection are explicit:

- if all supplied predicates accept the candidate, the candidate state is
  returned with `SspRk3Status::success`;
- if a predicate rejects the candidate and retries remain, the candidate is
  discarded, the timestep is halved, and the attempt restarts from the original
  input state;
- if retry capacity is exhausted, the solver returns
  `SspRk3Status::rejected_by_predicate`;
- if another retry is available but the next halved timestep would fall below
  `min_time_step`, the solver returns
  `SspRk3Status::minimum_timestep_reached`.

On SSP-RK3 rejection failure, the returned state is the original input state,
not the last rejected candidate. `time_step_used` records the accepted timestep
on success and the last rejected timestep on failure. `retries` records the
number of halving retries associated with the reported attempt.

Within generic SSP-RK3, the source driver evaluates configured
realizability checks through acceptance predicates. A `FinalOnly` setting
checks only the final raw-moment state; an `EveryStage` setting checks
intermediate SSP-RK3 raw-moment states and the final state. The solver remains
agnostic to the realizability mathematics.

The optional error estimate is diagnostic; any external timestep adaptation is
caller-owned. Phase-space evolution source updates use the predicate retry
behavior; they do not silently change the mathematical source model.

### Embedded Runge-Kutta schemes

The embedded Runge-Kutta schemes advances generic vector states from `t0` to
`t_bound`. Here the vector state is the packed collection of
quadrature node coordinates, and the RHS is built from `PhaseSpaceVelocity`.

The public API selects the built-in pair explicitly:

```cpp
const AdaptiveRkResult result =
    integrate_embedded_rk(initial_state, t0, t_bound, rhs, EmbeddedRkMethod::rk23, options);
```

The internal engine uses immutable built-in tableaus with stage abscissae,
strictly lower-triangular coefficients, higher and embedded weights, orders,
and FSAL metadata. It evaluates every tableau stage, forms both solutions, and
uses `y_high - y_low` as the error estimate. The accepted state is always the
higher-order solution. `rk23` is the four-stage Bogacki-Shampine 3(2) solver
[[Bogacki & Shampine, 1989](https://www.doi.org/10.1016/0893-9659(89)90079-7)];
`rk45` is the seven-stage Dormand-Prince RK5(4)7M solver [[Dormand & Prince,
1980](https://www.doi.org/10.1016/0771-050X(80)90013-3)]. Both are FSAL, so an
accepted final derivative is reused for the next trial in the same public
integration call. A rejected trial retains only the derivative at its unchanged
starting state.

Each trial step computes a scaled RMS error norm.

The error scale for component `i` is:

```text
scale_i = atol + max(abs(y_i), abs(y_new_i)) * rtol
```

The solver then forms an RMS norm of `error_i / scale_i`. A trial step is
accepted when the scaled error norm is less than `1`. Otherwise the trial step
is rejected and retried with a smaller step, unless the finite update was
already evaluated at the effective lower bound for the current interval.

`AdaptiveRkOptions` provides scalar tolerance and step-size controls for both
embedded methods:

- `relative_tolerance`;
- `absolute_tolerance`;
- `initial_step`;
- `max_step`;
- `min_step`;
- `max_steps`.

The default tolerances are:

```text
relative_tolerance = 1.0e-6
absolute_tolerance = 1.0e-9
```

`initial_step` defaults to `1.0e-6 s`, may be overridden, and must be finite
and strictly positive. QMOMKit does not estimate it from the state,
derivatives, or tolerances. `max_step == 0` means the full integration interval
is the maximum step size. `min_step == 0` means the solver uses a
floating-point spacing based floor near the current time. An explicit
`min_step` must not exceed the effective `max_step` (whether explicitly
configured or auto-derived from the integration interval); otherwise the
request is rejected as `AdaptiveRkStatus::invalid_input` before integration
begins. `max_steps` limits accepted steps, not rejected trial steps, and
remains an optional configuration override with an internal default. The actual
embedded-RK step is bounded by `min_step` and `max_step`, then clipped so
integration lands exactly on `t_bound`; the final remainder may therefore be
smaller than `min_step`. When the error controller would require a step below
the effective lower bound, the solver advances at that lower bound instead of
failing solely for tolerance reasons. Forward and backward integration are both
supported by applying the sign of `t_bound - t0` to the internal step
magnitude.

The controller exponent is derived from the embedded lower order: `1 / (p + 1)`
in magnitude, giving `1/3` for RK23 and `1/5` for RK45. Embedded-RK status
values are:

- `success`;
- `invalid_input`;
- `step_size_underflow`;
- `max_steps_exceeded`.

Invalid structural inputs and invalid tolerance or step-size options return
`AdaptiveRkStatus::invalid_input`. `AdaptiveRkStatus::step_size_underflow`
is retained for API compatibility and non-finite lower-bound updates, but a
finite update at the effective minimum step is accepted even when its error
estimate exceeds tolerance. If the number of accepted steps reaches `max_steps`
before `t_bound`, the solver returns `AdaptiveRkStatus::max_steps_exceeded`.

On embedded-RK failure, the result state and time are the last accepted state
and time. Diagnostics report accepted steps, rejected trial steps, RHS function
evaluations, final time, message, and the signed size of the last accepted
step. Node advection maps embedded-RK failure to
`NodeAdvectionStatus::integration_failed` and records the adaptive status
inside scheme-specific diagnostics. Phase-space evolution aggregates
scheme-neutral node-advection step, rejected-step, and RHS-evaluation counts.

## Phase-space evolution composition

The composed helpers combine the lower-level operations with a selected
operator-splitting scheme:

```cpp
advance_qmom_advection(...)
advance_qmom_moment_source(...)
advance_qmom_phase_space(...)
advance_cqmom_advection(...)
advance_cqmom_moment_source(...)
advance_cqmom_phase_space(...)
```

Each QMOM macro step starts from an accepted `QmomState`; the state's accepted
quadrature is used directly for phase-space advection, and its canonical raw
moments are used for moment-source evolution. `QmomState` exposes read-only
accessors for its configured maximum node count and reconstruction options, so
callers cannot mutate moments and quadrature independently or provide a second
reconstruction policy. It remains copyable and movable for application time
loops, but direct construction from independent moments and quadrature is
intentionally unavailable. Post-source QMOM reconstruction selects the next
accepted quadrature before moments are canonicalized from that quadrature.
CQMOM macro steps follow the same state-based pattern with `CqmomState`. If a
boundary contribution callback is supplied, it is added to the moment-source
rate inside the moment-source stage; it does not create another stage.

The public QMOM and CQMOM stage helpers expose the same advection and
moment-source-evolution stage implementations used by the composed helpers.
Standalone stage failure returns the stage input state unchanged. A complete
operator-split macro step remains transactional at the macro-step boundary: if
a later substage fails after an earlier substage succeeded, the returned result
state is still the state accepted before the complete split step began.

Successful CQMOM stage and macro-step states carry a private accepted
quadrature cache for the next CQMOM stage. A stage starting from a fresh
`CqmomCase::initial_state()` computes the initial CQMOM quadrature once before
applying the requested operation. A stage starting from a previously accepted
CQMOM phase-space evolution state reuses the cached accepted quadrature rather
than reinverting its stored moments. The cache is an implementation detail:
public callers pass `CqmomState` objects between stages and inspect the
`CqmomPhaseSpaceEvolutionResult::quadrature` value when they need the accepted
or diagnostic quadrature.

### Strang splitting

`PhaseSpaceEvolutionScheme::strang` is the default second-order Strang
splitting:

```text
advection half-step
moment reconstruction
moment-source full stage
post-source inversion
advection half-step
final moment reconstruction
```

The first advection half-step moves quadrature nodes and reconstructs raw
moments. The moment-source stage advances those moments with SSP-RK3, using the
combined source-plus-boundary rate when a boundary contribution is supplied.
The post-source inversion reconstructs quadrature from the updated moments.
The second advection half-step advances the new quadrature and reconstructs
final moments.

### Lie splitting

`PhaseSpaceEvolutionScheme::advsrc` applies first-order
advection-then-moment-source Lie splitting:

```text
advection full-step
moment reconstruction
moment-source full stage
post-source inversion
```

`PhaseSpaceEvolutionScheme::srcadv` applies first-order
moment-source-then-advection Lie splitting:

```text
moment-source full stage
post-source inversion
advection full-step
final moment reconstruction
```

The `advsrc` and `srcadv` schemes are useful for comparisons, diagnostics, or
applications that deliberately want a fixed first-order splitting order.

### Extending phase-space evolution schemes

Phase-space evolution split sequencing is centralized in the internal
phase-space evolution sequence executor. New phase-space evolution schemes or
time-integration patterns should extend that sequence code instead of
duplicating the existing `strang`, `advsrc`, and `srcadv` dispatch in each
representation workflow.

The shared phase-space evolution sequence should handle only:

- splitting and operation order;
- stage execution;
- common failure propagation;
- shared diagnostic records.

Representation-specific responsibilities remain in the individual phase-space
evolution paths:

- moment and quadrature reconstruction;
- QMOM and CQMOM inversion;
- dense CQMOM `RawMomentTensor` handling;
- lean CQMOM `CqmomMomentData` handling;
- source representation adaptation, including `CqmomSourceContext`;
- result and state construction, including `CqmomState` wrapping.

New schemes should not create a new generic phase-space evolution abstraction,
move QMOM or CQMOM representation details into the shared sequence code, or
require the dense and lean representation workflows to share storage-specific
code.

## Failure and status handling

Phase-space evolution APIs use explicit result and diagnostics objects. Callers
must check `ok()` or the status field before using returned moments or
quadrature.

Public phase-space evolution statuses are:

- `success`
- `invalid_input`
- `inversion_failed`
- `advection_failed`
- `moment_reconstruction_failed`
- `moment_source_failed`

Phase-space evolution diagnostics record:

- selected scheme;
- failed stage;
- post-source QMOM inversion status and initial/post-source CQMOM inversion statuses;
- CQMOM diagnostics when available;
- first and second phase-space advection statuses;
- node-advection step, rejected-step, and RHS-evaluation totals;
- moment reconstruction diagnostics;
- moment-source update status and SSP-RK3 retry count;
- message.

Failure handling is fail-fast. A failed standalone phase-space advection or
moment-source stage returns its stage input state. A failed complete
operator-split macro step returns the original macro-step input state. An
intermediate or failed-stage quadrature may be retained only as diagnostic
information and is not an accepted phase-space evolution state.

Status propagation follows subsystem ownership:

- input validation failures become `invalid_input`;
- failed initial or post-source reconstruction becomes `inversion_failed`;
- node-advection integration failures become `advection_failed`;
- raw-moment reconstruction exceptions become `moment_reconstruction_failed`;
- source callback, SSP-RK3 predicate rejection, or minimum-timestep failure
  becomes `moment_source_failed`.

Realizability interaction is explicit rather than hidden. Moment-source
candidate states may be rejected through SSP-RK3 acceptance predicates. After a
successful moment-source update, QMOM or CQMOM inversion is attempted again; if
the updated raw moments are not invertible under the configured reconstruction
options, phase-space evolution reports an inversion failure.

Adaptive CQMOM uses one context-dependent active-count policy. Initial
reconstruction starts each branch from its configured maximum and lets adaptive
Wheeler determine the active count. Phase-space advection reconstructs each
existing branch with its current active count as the ceiling, so counts may
remain or decrease but cannot grow. After an accepted moment-source update,
post-source inversion first reconstructs at the current active count; a
reduction is accepted immediately. If the current count is retained and the
configured maximum is larger, the branch tries exactly `Nactive + 1`, and
growth is accepted only when Wheeler returns that requested count. Branches
created by an invalidated parent topology start from configured maximum counts
and adapt downward.

After every successful inversion during CQMOM phase-space evolution through
`CqmomState`, the returned state is refreshed with the active counts extracted
from the returned quadrature. Structured phase-space evolution diagnostics
record branch path, configured maximum count, active count before and after the
stage, and whether a rank reduction or growth occurred.

Node advection preserves weights and moves an already accepted quadrature.
The moment reconstruction after node advection is a deterministic summation
from the moved quadrature; it is not a separate physical source model.

Optional QMOM and CQMOM abscissa bounds are retained as state reconstruction
metadata. Advection does not clip, project, or reject a node for crossing an
interval. Bounds are enforced only when a quadrature is reconstructed from
moments, including lazy source-stage reconstruction and post-source inversion;
a failed bounded inversion follows the ordinary stage or macro-step rollback
contract.

## Public API boundary

Public phase-space evolution and time-evolution entry points include:

- `advance_nodes(...)`
- `advance_moments_with_source(...)`
- `advance_qmom_advection(...)`
- `advance_qmom_moment_source(...)`
- `advance_qmom_phase_space(...)`
- `advance_cqmom_advection(...)`
- `advance_cqmom_moment_source(...)`
- `advance_cqmom_phase_space(...)`
- `ssp_rk3_step(...)`
- `integrate_embedded_rk(...)`

The main phase-space evolution helpers are:

```cpp
QmomPhaseSpaceEvolutionResult advance_qmom_advection(
    const QmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const QmomPhaseSpaceEvolutionOptions& options = {});

QmomPhaseSpaceEvolutionResult advance_qmom_moment_source(
    const QmomState& state,
    double t,
    double time_step,
    const MomentSource& source,
    const QmomPhaseSpaceEvolutionOptions& options = {});

QmomPhaseSpaceEvolutionResult advance_qmom_moment_source(
    const QmomState& state,
    double t,
    double time_step,
    const MomentSource& source,
    const MomentBoundaryContribution& boundary,
    const QmomPhaseSpaceEvolutionOptions& options);

QmomPhaseSpaceEvolutionResult advance_qmom_phase_space(
    const QmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const MomentSource& source,
    const QmomPhaseSpaceEvolutionOptions& options = {});

QmomPhaseSpaceEvolutionResult advance_qmom_phase_space(
    const QmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const MomentSource& source,
    const MomentBoundaryContribution& boundary,
    const QmomPhaseSpaceEvolutionOptions& options);

CqmomPhaseSpaceEvolutionResult advance_cqmom_advection(
    const CqmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const CqmomPhaseSpaceEvolutionOptions& options = {});

CqmomPhaseSpaceEvolutionResult advance_cqmom_moment_source(
    const CqmomState& state,
    double t,
    double time_step,
    const CqmomMomentSource& source,
    const CqmomPhaseSpaceEvolutionOptions& options = {});

CqmomPhaseSpaceEvolutionResult advance_cqmom_moment_source(
    const CqmomState& state,
    double t,
    double time_step,
    const CqmomMomentSource& source,
    const CqmomMomentBoundaryContribution& boundary,
    const CqmomPhaseSpaceEvolutionOptions& options);

CqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(
    const CqmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const CqmomMomentSource& source,
    const CqmomPhaseSpaceEvolutionOptions& options = {});

CqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(
    const CqmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const CqmomMomentSource& source,
    const CqmomMomentBoundaryContribution& boundary,
    const CqmomPhaseSpaceEvolutionOptions& options);
```

The CQMOM source callback used by phase-space evolution through `CqmomState`
has the shape:

```cpp
using CqmomMomentSource =
    std::function<CqmomMomentRate(double, const CqmomSourceContext&)>;
```

The CQMOM boundary contribution callback has the same context shape and return
ordering:

```cpp
using CqmomMomentBoundaryContribution =
    std::function<CqmomMomentRate(double, const CqmomSourceContext&)>;
```

The boundary-aware phase-space evolution overloads require the options argument
explicitly. This keeps overload resolution unambiguous while preserving default
options on the source-only overloads.

`CqmomSourceContext` exposes the CQMOM moment layout selected for phase-space
evolution. For equivalent states, dense-backed and lean-backed `CqmomState`
phase-space evolution must present identical context indices, values, ordering,
hierarchy, and lookup behavior to `CqmomMomentSource` and
`CqmomMomentBoundaryContribution`. Dense storage may retain additional raw
tensor moments internally, but those dense-only entries are not part of the
callback context unless they are also in the phase-space evolution CQMOM moment
layout.

CQMOM moment-source callbacks update the canonical phase-space evolution moment
layout. After a successful moment-source stage inversion, the accepted
quadrature becomes the canonical branch representation for the returned state.
Lean storage reconstructs the retained CQMOM layout from that quadrature for
every accepted post-source state, even when the active topology is unchanged.
Dense storage reconstructs the complete configured dense tensor from that
quadrature, including higher-order or mixed tensor entries that were not
exposed to the source callback. Stored post-source moments are therefore the
moments implied by the accepted quadrature, not independently retained source
variables. In adaptive mode, reduction or growth can replace source-candidate
higher-order information with the closure represented by the accepted
quadrature; the returned active branch counts are derived from that same
quadrature.

With coordinate support configured, the post-source invariant is:

```text
source candidate moments
    -> bounded recursive inversion
    -> accepted quadrature
    -> accepted active counts
    -> canonical stored moments
```

Bounded recursive inversion is phase-space-coordinate based, independent of the
conditioning level at which that coordinate is reconstructed. The accepted
quadrature is the only authority for the returned topology and stored moments.
This matters for lean storage because lean states carry only the phase-space
evolution CQMOM moment layout; after a reduction or rejected growth trial,
retaining the source candidate values would leave high-order entries
inconsistent with the accepted bounded quadrature.

Fixed-count bounded failures fail the stage and return the unchanged input
state. Adaptive bounded reconstruction may reduce a branch to the largest count
whose reconstructed abscissas satisfy the configured bounds. It may also try
one-node growth toward the configured maximum; if that growth candidate returns
an abscissa outside the configured bounds, the current accepted count is
retained. When a parent branch changes count, descendant branches and
active-count metadata are rebuilt from the accepted quadrature. Coordinate
support metadata is preserved across these state transitions. Failed
post-source inversions roll back both stored moments and active counts.

The public CQMOM phase-space evolution API is `CqmomState` based. The
high-level `advance_cqmom_phase_space(const CqmomState&, ...)` helper applies
the configured splitting scheme. The public `advance_cqmom_advection(...)` and
`advance_cqmom_moment_source(...)` stage helpers are for applications that own
stage orchestration or acceptance externally.

`CqmomState` stores the selected dense or lean moment representation,
reconstruction policy, current active branch-count metadata, and, after
successful evolution, a private accepted quadrature cache while remaining
representation neutral to callers. In adaptive mode, the hierarchy node counts
are configured maxima. `CqmomState::active_branch_counts()` reports the current
ragged topology using conditioning-level parent paths. The returned
`CqmomQuadrature` may contain adaptive ragged terminal storage, and successful
CQMOM evolution through `CqmomState` stores that topology for the next stage or
macro step.

Advection CQMOM stages preserve accepted terminal topology and active branch
counts, use only node-advection options from `CqmomPhaseSpaceEvolutionOptions`,
and reconstruct stored moments by summation from the moved quadrature. They do
not perform a post-stage CQMOM inversion. CQMOM moment-source stages use only
moment-source integration options, update the canonical CQMOM moment layout,
and apply the case-owned recursive reconstruction policy. Adaptive post-source
reconstruction may reduce active branch counts and may grow a retained branch
by at most one node toward the configured maximum. Standalone stage helpers
ignore the split-scheme option because the caller owns composition order.

`CqmomPhaseSpaceEvolutionOptions` contains only stage policies. CQMOM inversion
policies such as Wheeler options, conditional reconstruction options, adaptive
branch counts, diagnostic level, and Vandermonde solver are case-owned policies
preserved by `CqmomState`. `CqmomDenseCase` and `CqmomLeanCase` remain public
advanced/static inversion inputs, but public CQMOM evolution remains the
storage-independent `CqmomState` API.

The lower-level public helpers are:

- `advance_nodes(...)`: phase-space advection with unchanged
  weights;
- `advance_moments_with_source(...)`: moment-source evolution of raw moments through
  SSP-RK3;
- `ssp_rk3_step(...)`: generic vector SSP-RK3 utility;
- `integrate_embedded_rk(...)`: generic adaptive vector RK23/RK45 utility.

Implementation details include:

- node packing and unpacking;
- velocity callback wrappers;
- source callback wrappers;
- evolution stage lambdas;
- moment reconstruction handoff;
- diagnostic aggregation;
- timestep controller internals.

These details support the public APIs but should not become required concepts
for projects that use QMOMKit.

## Batch convenience wrappers

The accepted-state batch functions are serial convenience wrappers
only. Each wrapper reserves one result per input and delegates in increasing
input order to its corresponding public scalar accepted-state operation. It
must not call lower-level stages, reconstruction helpers, or time
integrators directly, and it must not introduce a second numerical path.

QMOM and CQMOM provide matching batch families for advection stages,
moment-source-evolution stages, and complete steps, including boundary-aware
overloads. Ordinary batch functions use the existing scalar callback types and
the `_batch` suffix. Batch-index-aware functions use the `_batch_with_index`
suffix and callbacks whose first argument is `batch_index`, the zero-based
position of the state in the supplied collection. The distinct names select the
callback family without callable-conversion ambiguity. Each scalar callback
adapter captures that position by value, forwards all remaining arguments
unchanged, and exists only for the duration of the scalar call. In particular,
`CqmomSourceContext` remains non-owning and stage-local.

The input/output contract is exactly a serial loop over scalar calls: an empty
batch performs no validation or callback invocation; ordinary unsuccessful
results are appended unchanged; and an escaping exception is not caught,
wrapped, or annotated. Results, adaptive node-advection substeps, source
retries, active node counts, active CQMOM branch topology, and diagnostics are
strictly entry-local.

Batch evolution does not add execution policies, threading, scheduling,
parallelism, mesh ownership, per-entry operation controls, batch diagnostics,
or persistent batch state. Hosts retain responsibility for spatial organization
and may instead call the scalar functions inside their own parallel loops.

## Validation expectations

Phase-space evolution changes should be validated through analytical or
manufactured cases that exercise the relevant representation path. For a
complete macro step this is:

```text
moments
    -> QMOM/CQMOM inversion
    -> node advection
    -> raw-moment reconstruction
    -> moment-source update
    -> post-source inversion
    -> final node advection/reconstruction
```

Existing verification covers:

- node advection against analytical advection trajectories;
- moment-source updates against analytical moment evolution;
- complete QMOM evolution for commuting and noncommuting operators;
- complete CQMOM evolution for a compact two-dimensional case;
- Strang and Lie scheme behavior;
- reconstruction error and zeroth-moment checks;
- diagnostics for successful runs.

Changes to splitting order, status mapping, callback validation,
realizability hooks, or diagnostic fields are design changes and
require matching verification coverage.
