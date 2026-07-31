# Transport Design

This page documents the current QMOMKit design for transport and time
integration. It is a developer reference for maintainers, contributors, and
AI-assisted development workflows.

For raw-moment, tensor, phase-space, node, abscissa, and weight conventions,
see `numerical_conventions.md`. For QMOM and CQMOM reconstruction behavior,
see `qmom_wheeler.md` and `cqmom_recursive.md`.

## Transport Architecture Overview

QMOMKit transport is built from three numerical operations:

- node advection: continuous evolution of quadrature abscissas in phase space;
- moment source evolution: direct evolution of raw moments by a caller-supplied
  source term;
- transport integration: composition of reconstruction, node advection, moment
  reconstruction, source evolution, and post-source inversion.

The library provides transport utilities. Users provide physical models:

- phase-space velocity callbacks;
- moment-source callbacks;
- boundary contribution callbacks;
- outer timestep loops;
- physical interpretation of coordinates and moment components.

QMOMKit does not own meshes, cells, finite-volume or finite-element operators,
boundary conditions, MPI communication, particle models, or application
physics. Transport operates on accepted QMOM/CQMOM state, quadrature nodes, raw
moments, and callback interfaces only.

The callback categories are independent runtime inputs. A workflow may use
only continuous phase-space motion, only moment-space source evolution, source
plus boundary exchange, velocity plus source, or all three categories. Case and
State objects store numerical configuration and evolving moment data; they do
not store physical models or callback objects.

The production transport composition has two numerical stages: continuous
motion is represented by advancing quadrature nodes, while discontinuous or
source-like effects are represented by advancing raw moments. Direct moment
transport is not part of the supported QMOMKit transport design.

## Node Advection Design

Node advection advances quadrature abscissas through a continuous phase-space
velocity field. It preserves quadrature weights and terminal paths.

The public velocity callback is:

```cpp
using PhaseSpaceVelocity =
    std::function<std::vector<double>(double, const std::vector<double>&)>;
```

The callback receives:

- current time;
- one node coordinate vector ordered by phase-space dimension.

It returns `dxi/dt` with the same number of components. QMOMKit validates the
returned size and finiteness, but it does not assign physical meaning to any
coordinate. A non-finite returned component is reported with the evaluation
time, zero-based phase-space dimension, offending value, input coordinate
vector, and complete returned velocity vector.

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
the selected `NodeAdvectionScheme`. The default is fixed-stage SSP-RK3, which
advances exactly over the supplied node-advection timestep. Adaptive RK23 is
available when trajectory accuracy should control internal substeps. After
integration, the state is unpacked back into the same terminal-node
representation. Weights are not integrated by node advection.

After a successful continuous node step, transport reconstructs raw moments
from the moved quadrature:

```text
m_alpha = sum_j W_j product_r xi_(r,j)^alpha_r
```

This reconstruction is the handoff from continuous node motion back to the
moment representation used by source evolution and later inversions.

For adaptive CQMOM, the phase-space hierarchy node counts remain the
configured maximum branch capacity throughout transport. The current active
branch counts are separate runtime state: they are extracted from the ragged
terminal paths of each successful adaptive CQMOM quadrature and stored on
`CqmomState` beside the current moments. Continuous node advection preserves
terminal paths and therefore preserves active branch topology. Moment
reconstruction updates values only; it must not replace configured maximum
counts with active counts or discard the active-count metadata.

## Moment Source And Boundary Contribution Design

Moment source evolution advances raw moments directly:

```text
dm/dt = source(t, m)
```

The public source callback is:

```cpp
using MomentSource =
    std::function<std::vector<double>(double, const std::vector<double>&)>;
```

The callback receives the current time and a caller-owned raw-moment vector.
It must return a derivative vector with the same size. QMOMKit validates the
source result through the generic ODE solver path, but it does not know what
individual moment components mean.

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

For QMOM, the vector is the one-dimensional raw moment sequence. CQMOM transport
uses the CQMOM-native source callback over a `CqmomSourceContext` view of the
current source moments, so physical CQMOM source models do not need to depend on
dense tensor or lean required-moment storage. Equivalent dense-backed and
lean-backed `CqmomState` objects expose the same selected CQMOM transport moment
layout; dense-only tensor entries remain internal unless they are part of that
layout.

The CQMOM source context exposes only:

- the phase-space hierarchy;
- the SSP-RK stage time;
- the number of source moment entries;
- the raw moment multi-index for each entry;
- the raw moment value for each entry.

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
contribution callback. The transport code advances the numerical state and
reports failures; it does not define model equations.

Callbacks are combined only when the selected public routine receives them.
Lower-level workflows can call node advection or source evolution separately.
Full transport helpers compose the continuous and source-like parts requested
by their overloads; absent source-like physics is represented explicitly by a
zero-rate source callback where the helper signature requires one.

Moment source updates use SSP-RK3. Optional acceptability predicates in
`SspRk3Options` allow callers to reject candidate moment states, including
states that violate externally defined realizability policies.

When both a source and a boundary contribution are supplied, transport builds a
single effective source rate for each SSP-RK3 stage:

```text
effective_rate(t, m) = source(t, m) + boundary(t, m)
```

This preserves the existing source-update retry logic, acceptance predicates,
diagnostics, and post-source inversion behavior. Boundary callback exceptions
are wrapped before they enter the existing source failure pathway so diagnostics
identify the boundary contribution callback as the origin. Source callback
exceptions keep the existing source failure behavior.

## Time Integration

QMOMKit exposes generic Runge-Kutta utilities, but transport uses them in
different numerical roles.

### Node-Advection Integration

Node advection advances quadrature coordinates through:

```text
dxi/dt = velocity(t, xi)
```

The public `NodeAdvectionScheme` selects the trajectory integrator:

- `ssp_rk3`: the default fixed-stage SSP-RK3 trajectory step;
- `adaptive_rk23`: adaptive Bogacki-Shampine RK23 with scalar tolerances.

The choice is a user-facing trajectory-integration policy, not a change to the
transport composition. SSP-RK3 is appropriate when the surrounding application
controls the timestep and needs fixed callback cost. Adaptive RK23 is
appropriate when trajectory accuracy inside the requested transport step should
be controlled by internal error estimation. Neither scheme is universally
preferred; the useful choice depends on the application timestep policy and
velocity-field behavior.

The SSP-RK3 node path is implemented inside node advection, not through the
moment-source update wrapper. It does not use source acceptability predicates,
realizability checks, source retry semantics, or source-specific performance
counters. It performs one three-stage step over the supplied timestep and
reports scheme-neutral node-advection diagnostics.

Adaptive RK23 may take internal accepted and rejected trial steps to reach the
requested end time. It is useful when QMOMKit should control trajectory
accuracy inside a larger transport macro step.

Configuration-driven drivers load the full `NodeAdvectionOptions` object from
`[time_integration.node_advection]` and the nested
`[time_integration.node_advection.rk23]` section. Startup reporting should use
the shared node-advection formatter so fixed-step runs report only the active
scheme, while RK23 runs also report effective tolerances, automatic step-size
controls, and the maximum internal accepted-step count.

Node advection does not perform realizability checks. It moves already
reconstructed quadrature coordinates and preserves weights; realizability is
checked through quadrature inversion or caller-supplied moment-state policies
after reconstruction or source updates.

For adaptive CQMOM transport through `CqmomState`, active branch counts are part
of accepted state metadata. `CqmomCase::initial_state()` seeds that metadata from
the initial trusted inversion and fails state construction if the initial
adaptive topology cannot be reconstructed. Later state inversions validate
supplied active-count metadata strictly: there must be exactly one reachable
entry per active branch, configured node-count metadata must match the hierarchy,
active counts must be within configured capacity, and stale or unreachable
descendant entries are rejected. During a source-stage inversion, if an ancestor
branch changes active count, descendant limits below the new branch are
invalidated top-down and recomputed from configured capacity so numeric path
reuse does not apply stale branch limits after a collapse, merge, or growth.
Diagnostics report only true node-count transitions and distinguish reduced,
grown, created, and removed branch records.

`CqmomTransportOptions` owns transport-stage policies: node-advection/RK23
controls, moment-source SSP-RK3 controls, and transport sequencing. Conditional
moment controls, standalone QMOM controls, CQMOM Wheeler controls,
adaptive branch-count enablement, CQMOM diagnostic level, and Vandermonde solver
selection are case-owned policies stored in `CqmomState`. Dense-backed and
lean-backed states must pass these policies through to the same recursive
CQMOM, node-advection, moment-reconstruction, and source-update algorithms;
accepting a configured case policy and replacing it with `false`, `detailed`,
`automatic`, or another aggregate default is a bug.

### Source SSP-RK3

Source SSP-RK3 advances raw moment vectors through moment-source updates. It is
an explicit third-order strong-stability-preserving Runge-Kutta step over a
generic vector state.

Transport uses this path for the discontinuous or source-like stage because
that stage acts directly on moments:

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
realizability, particles, or phase-space coordinates. Realizability checks
enter only through caller-supplied predicates or through subsequent
quadrature inversion.

## Time Integration Details

The public Runge-Kutta utilities are generic vector ODE solvers. Transport
uses them through node-advection and moment-source wrappers where appropriate,
but the solvers themselves do not know about QMOM, CQMOM, realizability,
particles, phase-space coordinates, or physical model semantics.

Both solvers use the callback convention:

```cpp
using OdeRhs =
    std::function<std::vector<double>(double, const std::vector<double>&)>;
```

The callback receives `(t, state)` and must return a derivative vector with the
same size as `state`. A wrong-sized derivative or non-finite derivative is an
error. Callback exceptions are not converted into solver-specific status values;
they follow the low-level solver behavior and propagate out of the generic RK
utility.

Transport-facing wrappers add domain-specific failure handling. Node advection
catches phase-space velocity failures and reports
`NodeAdvectionStatus::velocity_failed`; trajectory-integrator failures report
`NodeAdvectionStatus::integration_failed`, with optional RK23 status retained
when adaptive RK23 is selected. Moment-source updates catch source callback or
SSP-RK3 validation failures and report
`MomentSourceUpdateStatus::source_failed`. The full transport helpers then map
failed node advection to `TransportStatus::continuous_failed` and failed
moment-source updates to `TransportStatus::source_failed`.

### SSP-RK3 Details

SSP-RK3 advances one raw-moment state over a requested timestep:

```text
y' = f(t, y)
```

Conceptually, the step evaluates three explicit stages:

```text
stage 1: forward Euler candidate from y_n
stage 2: convex combination of y_n and a second Euler-like update
final:   convex combination producing the third-order state y_(n+1)
```

The convex-combination structure is the reason SSP-RK3 is used for moment-source
evolution: source updates can be paired with caller-supplied acceptability
predicates while preserving a simple explicit update structure. QMOMKit uses
this path for discontinuous or source-like moment evolution, where the state is
a raw-moment vector rather than node coordinates.

`SspRk3Options` controls optional solver behavior:

- `estimate_error` requests an embedded lower-order error estimate;
- `accept` checks the completed candidate state;
- `stage_accept` checks stage-1, stage-2, and final candidate states;
- `max_retries` limits retry attempts after predicate rejection;
- `min_time_step` prevents retrying below a configured timestep floor.

If error estimation is enabled, the solver compares the final third-order state
with an embedded lower-order state. The result is a relative Euclidean norm when
the final state norm is nonzero, and an absolute difference norm when the final
state norm is zero. This estimate is diagnostic and does not by itself accept or
reject a step.

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

Realizability policies interact with SSP-RK3 only through predicates. A
`FinalOnly` policy checks only the final candidate state; an `EveryStage` policy
checks intermediate SSP-RK3 stage states and the final state. The solver remains
agnostic to the meaning of those predicates.

The helper `adapt_ssp_rk3_timestep(...)` is a separate timestep-controller
utility for positive nonzero error estimates. It applies a safety factor and
clamps the suggested timestep between caller-supplied minimum and maximum
values. Transport source updates use the predicate retry behavior; they do not
silently change the mathematical source model.

### RK23 Details

Adaptive RK23 advances generic vector states from `t0` to `t_bound`. In
transport, the vector state is the packed collection of quadrature node
coordinates, and the RHS is built from `PhaseSpaceVelocity`.

The implementation uses a Bogacki-Shampine 3(2) embedded pair. Each trial step
computes:

- a third-order candidate state;
- an embedded lower-order error estimate;
- a scaled RMS error norm.

The error scale for component `i` is:

```text
scale_i = atol + max(abs(y_i), abs(y_new_i)) * rtol
```

The solver then forms an RMS norm of `error_i / scale_i`. A trial step is
accepted when the scaled error norm is less than `1`. Otherwise the trial step
is rejected and retried with a smaller step, unless the finite update was
already evaluated at the effective lower bound for the current interval.

`AdaptiveRk23Options` provides scalar tolerance and step-size controls:

- `relative_tolerance`;
- `absolute_tolerance`;
- optional `initial_step`;
- `max_step`;
- `min_step`;
- `max_steps`.

The default tolerances are:

```text
relative_tolerance = 1.0e-6
absolute_tolerance = 1.0e-9
```

`max_step == 0` means the full integration interval is the maximum step size.
`min_step == 0` means the solver uses a floating-point spacing based floor near
the current time. `max_steps` limits accepted steps, not rejected trial steps,
and remains an optional configuration override with an internal default. The
actual RK23 step is bounded by `min_step` and `max_step`, then clipped so
integration lands exactly on `t_bound`; the final remainder may therefore be
smaller than `min_step`. When the error controller would require a step below
the effective lower bound, the solver advances at that lower bound instead of
failing solely for tolerance reasons. Forward and backward integration are both
supported by applying the sign of `t_bound - t0` to the internal step magnitude.

RK23 status values are:

- `success`;
- `invalid_input`;
- `step_size_underflow`;
- `max_steps_exceeded`.

Invalid structural inputs and invalid tolerance or step-size options return
`AdaptiveRk23Status::invalid_input`. `AdaptiveRk23Status::step_size_underflow`
is retained for API compatibility and non-finite lower-bound updates, but a
finite update at the effective minimum step is accepted even when its error
estimate exceeds tolerance. If the number of accepted steps reaches `max_steps`
before `t_bound`, the solver returns `AdaptiveRk23Status::max_steps_exceeded`.

On RK23 failure, the result state and time are the last accepted state and time.
Diagnostics report accepted steps, rejected trial steps, RHS function
evaluations, final time, message, and the signed size of the last accepted step.
Node advection maps RK23 failure to
`NodeAdvectionStatus::integration_failed` and records the adaptive RK23 status
inside scheme-specific diagnostics. Transport aggregates scheme-neutral
node-advection step, rejected-step, and RHS-evaluation counts.

## Transport Composition

The full transport helpers compose the lower-level operations with a selected
operator-splitting scheme:

```cpp
advance_qmom_continuous(...)
advance_qmom_discontinuous(...)
advance_transport_qmom(...)
advance_cqmom_continuous(...)
advance_cqmom_discontinuous(...)
advance_transport_cqmom(...)
```

Each QMOM macro step starts from an accepted `QmomState`; the state's accepted
quadrature is used directly for continuous node motion, and its canonical raw
moments are used for discontinuous moment-source evolution. `QmomState` exposes
read-only accessors so callers cannot mutate moments and quadrature
independently. It remains copyable and movable for application time loops, but
direct construction from independent moments and quadrature is intentionally
unavailable. Post-source QMOM reconstruction selects the next accepted
quadrature before moments are canonicalized from that quadrature. CQMOM macro
steps follow the same state-based pattern with `CqmomState`. If a boundary
contribution callback is supplied, it is added to the source rate inside the
discontinuous source step; it does not create another transport stage.

The public QMOM and CQMOM stage helpers expose the same continuous and
discontinuous stage implementations used by the full transport helpers.
Standalone stage failure returns the stage input state unchanged. A complete
operator-split macro step remains transactional at the macro-step boundary: if
a later substage fails after an earlier substage succeeded, the returned result
state is still the state accepted before the complete split step began.

Successful CQMOM stage and macro-step states carry a private accepted
quadrature cache for the next CQMOM transport stage. A stage starting from a
fresh `CqmomCase::initial_state()` computes the initial CQMOM quadrature once
before applying the requested operation. A stage starting from a previously
accepted CQMOM transport state reuses the cached accepted quadrature rather
than reinverting its stored moments. The cache is an implementation detail:
public callers pass `CqmomState` objects between stages and inspect the
`CqmomTransportResult::quadrature` value when they need the accepted or
diagnostic quadrature.

### Strang Splitting

`TransportScheme::strang` is the default second-order Strang splitting
composition:

```text
continuous half-step
moment reconstruction
source full-step
post-source inversion
continuous half-step
final moment reconstruction
```

The first continuous half-step moves quadrature nodes and reconstructs raw
moments. The source step advances those moments with SSP-RK3, using the
combined source-plus-boundary rate when a boundary contribution is supplied.
The post-source inversion reconstructs quadrature from the updated moments.
The second continuous half-step advances the new quadrature and reconstructs
final moments.

### `cd`

`TransportScheme::cd` applies first-order continuous-then-discontinuous Lie
splitting:

```text
continuous full-step
moment reconstruction
source full-step
post-source inversion
```

### `dc`

`TransportScheme::dc` applies first-order discontinuous-then-continuous Lie
splitting:

```text
source full-step
post-source inversion
continuous full-step
final moment reconstruction
```

The `cd` and `dc` schemes are useful for comparisons, diagnostics, or
applications that deliberately want a fixed first-order splitting order. They
do not replace the default Strang splitting composition.

### Extending Transport Schemes

Transport split sequencing is centralized in the internal transport sequence
executor. New transport schemes or time-integration patterns should extend that
sequence code instead of duplicating the existing `strang`, `cd`, and `dc`
dispatch in each representation workflow.

The shared transport sequence should handle only:

- splitting and operation order;
- stage execution;
- common failure propagation;
- shared diagnostic records.

Representation-specific responsibilities remain in the individual transport
paths:

- moment and quadrature reconstruction;
- QMOM and CQMOM inversion;
- dense CQMOM `RawMomentTensor` handling;
- lean CQMOM `CqmomMomentData` handling;
- source representation adaptation, including `CqmomSourceContext`;
- result and state construction, including `CqmomState` wrapping.

New schemes should not create a new generic transport abstraction, move QMOM or
CQMOM representation details into the shared sequence code, or require the
dense and lean representation workflows to share storage-specific code.

## Failure And Status Handling

Transport APIs use explicit result and diagnostics objects. Callers must check
`ok()` or the status field before using returned moments or quadrature.

Public transport statuses are:

- `success`
- `invalid_input`
- `inversion_failed`
- `continuous_failed`
- `moment_reconstruction_failed`
- `source_failed`

Transport diagnostics record:

- selected scheme;
- failed stage;
- post-source QMOM inversion status and initial/post-source CQMOM inversion statuses;
- CQMOM diagnostics when available;
- first and second continuous node-advection statuses;
- node-advection step, rejected-step, and RHS-evaluation totals;
- moment reconstruction diagnostics;
- source update status and SSP-RK3 retry count;
- message.

Failure handling is fail-fast. A failed standalone continuous or discontinuous
stage returns its stage input state. A failed complete operator-split macro
step returns the original macro-step input state. An intermediate or
failed-stage quadrature may be retained only as diagnostic information and is
not an accepted transport state.

Status propagation follows subsystem ownership:

- input validation failures become `invalid_input`;
- failed initial or post-source reconstruction becomes `inversion_failed`;
- node-advection integration failures become `continuous_failed`;
- raw-moment reconstruction exceptions become `moment_reconstruction_failed`;
- source callback, SSP-RK3 predicate rejection, or minimum-timestep failure
  becomes `source_failed`.

Realizability interaction is explicit rather than hidden. Moment-source
candidate states may be rejected through SSP-RK3 acceptance predicates. After
a successful source update, QMOM or CQMOM inversion is attempted again; if the
updated raw moments are not invertible under the configured reconstruction
options, transport reports an inversion failure.

Adaptive CQMOM uses one context-dependent active-count policy. Initial
reconstruction starts each branch from its configured maximum and lets adaptive
Wheeler determine the active count. Continuous transport reconstructs each
existing branch with its current active count as the ceiling, so counts may
remain or decrease but cannot grow. After an accepted discontinuous source
update, post-source inversion first reconstructs at the current active count; a
reduction is accepted immediately. If the current count is retained and the
configured maximum is larger, the branch tries exactly `Nactive + 1`, and growth
is accepted only when Wheeler returns that requested count. Branches created by
an invalidated parent topology start from configured maximum counts and adapt
downward.

After every successful inversion during CQMOM transport through `CqmomState`,
the returned state is refreshed with the active counts extracted from the
returned quadrature. Structured transport diagnostics record branch path,
configured maximum count, active count before and after the stage, and whether
a rank reduction or growth occurred.

Node advection preserves weights and moves an already accepted quadrature.
The moment reconstruction after node advection is a deterministic summation
from the moved quadrature; it is not a separate physical source model.

## Public API Boundary

Public transport and time-evolution entry points include:

- `advance_nodes(...)`
- `advance_moments_with_source(...)`
- `advance_qmom_continuous(...)`
- `advance_qmom_discontinuous(...)`
- `advance_transport_qmom(...)`
- `advance_cqmom_continuous(...)`
- `advance_cqmom_discontinuous(...)`
- `advance_transport_cqmom(...)`
- `ssp_rk3_step(...)`
- `integrate_rk23(...)`

The main transport helpers are:

```cpp
QmomTransportResult advance_qmom_continuous(
    const QmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const QmomTransportOptions& options = {});

QmomTransportResult advance_qmom_discontinuous(
    const QmomState& state,
    std::size_t maximum_node_count,
    double t,
    double time_step,
    const MomentSource& source,
    const QmomTransportOptions& options = {});

QmomTransportResult advance_qmom_discontinuous(
    const QmomState& state,
    std::size_t maximum_node_count,
    double t,
    double time_step,
    const MomentSource& source,
    const MomentBoundaryContribution& boundary,
    const QmomTransportOptions& options);

QmomTransportResult advance_transport_qmom(
    const QmomState& state,
    std::size_t maximum_node_count,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const MomentSource& source,
    const QmomTransportOptions& options = {});

QmomTransportResult advance_transport_qmom(
    const QmomState& state,
    std::size_t maximum_node_count,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const MomentSource& source,
    const MomentBoundaryContribution& boundary,
    const QmomTransportOptions& options);

CqmomTransportResult advance_cqmom_continuous(
    const CqmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const CqmomTransportOptions& options = {});

CqmomTransportResult advance_cqmom_discontinuous(
    const CqmomState& state,
    double t,
    double time_step,
    const CqmomMomentSource& source,
    const CqmomTransportOptions& options = {});

CqmomTransportResult advance_cqmom_discontinuous(
    const CqmomState& state,
    double t,
    double time_step,
    const CqmomMomentSource& source,
    const CqmomMomentBoundaryContribution& boundary,
    const CqmomTransportOptions& options);

CqmomTransportResult advance_transport_cqmom(
    const CqmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const CqmomMomentSource& source,
    const CqmomTransportOptions& options = {});

CqmomTransportResult advance_transport_cqmom(
    const CqmomState& state,
    double t,
    double time_step,
    const PhaseSpaceVelocity& velocity,
    const CqmomMomentSource& source,
    const CqmomMomentBoundaryContribution& boundary,
    const CqmomTransportOptions& options);
```

The CQMOM source callback used by transport through `CqmomState` has the shape:

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

The boundary-aware transport overloads require the options argument explicitly.
This keeps overload resolution unambiguous while preserving default options on
the source-only overloads.

`CqmomSourceContext` exposes the CQMOM moment layout selected for transport.
For equivalent states, dense-backed and lean-backed `CqmomState` transport must
present identical context indices, values, ordering, hierarchy, and lookup
behavior to `CqmomMomentSource` and `CqmomMomentBoundaryContribution`. Dense
storage may retain additional raw tensor moments internally, but those dense-only
entries are not part of the callback context unless they are also in the
transport CQMOM moment layout.

CQMOM source callbacks update the canonical transport-layout moments. After a
successful source-stage inversion, the accepted quadrature becomes the canonical
branch representation for the returned state. Lean storage reconstructs the
retained CQMOM layout from that quadrature for every accepted post-source
state, even when the active topology is unchanged. Dense storage reconstructs
the complete configured dense tensor from that quadrature, including
higher-order or mixed tensor entries that were not exposed to the source
callback. Stored post-source moments are therefore the moments implied by the
accepted quadrature, not independently retained source variables. In adaptive
mode, reduction or growth can replace source-candidate higher-order information
with the closure represented by the accepted quadrature; the returned active
branch counts are derived from that same quadrature.

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
This matters for lean storage because lean states carry only the transport CQMOM
moment layout; after a reduction or rejected growth trial, retaining the source
candidate values would leave high-order entries inconsistent with the accepted
bounded quadrature.

Fixed-count bounded failures fail the stage and return the unchanged input
state. Adaptive bounded reconstruction may reduce a branch to the largest
accepted feasible count. It may also try one-node growth toward the configured
maximum; if that growth candidate is support-infeasible, the current feasible
count is retained. When a parent branch changes count, descendant branches and
active-count metadata are rebuilt from the accepted quadrature. Coordinate
support metadata is preserved across these state transitions. Failed
post-source inversions roll back both stored moments and active counts.

The public CQMOM transport API is `CqmomState` based. The high-level
`advance_transport_cqmom(const CqmomState&, ...)` helper owns the configured
C/D sequence. The public `advance_cqmom_continuous(...)` and
`advance_cqmom_discontinuous(...)` stage helpers are for applications that own
stage orchestration or acceptance externally.

`CqmomState` stores the selected dense or lean moment representation,
reconstruction policy, current active branch-count metadata, and, after
successful transport, a private accepted quadrature cache while remaining
representation neutral to callers. In adaptive mode, the hierarchy node counts
are configured maxima. `CqmomState::active_branch_counts()` reports the current
ragged topology using conditioning-level parent paths. The returned
`CqmomQuadrature` may contain adaptive ragged terminal storage, and successful
CQMOM transport through `CqmomState` stores that topology for the next stage or
macro step.

Continuous CQMOM stages preserve accepted terminal topology and active branch
counts, use only node-advection options from `CqmomTransportOptions`, and
reconstruct stored moments by summation from the moved quadrature. They do not
perform a post-stage CQMOM inversion. Discontinuous CQMOM stages use only source
integration options, update the canonical CQMOM transport moment layout, and
apply the case-owned recursive reconstruction policy. Adaptive post-source
reconstruction may reduce active branch counts and may grow a retained branch by
at most one node toward the configured maximum. Standalone stage helpers ignore
the split-scheme option because the caller owns composition order.

`CqmomTransportOptions` contains only transport-stage policies. CQMOM inversion
policies such as Wheeler options, conditional reconstruction options, adaptive
branch counts, diagnostic level, and Vandermonde solver are case-owned policies
preserved by `CqmomState`. `CqmomDenseCase` and `CqmomLeanCase` remain public
advanced/static inversion inputs, but public CQMOM transport remains the
storage-independent `CqmomState` API.

The lower-level public helpers are:

- `advance_nodes(...)`: continuous phase-space node advection with unchanged
  weights;
- `advance_moments_with_source(...)`: raw-moment source evolution through
  SSP-RK3;
- `ssp_rk3_step(...)`: generic vector SSP-RK3 utility;
- `integrate_rk23(...)`: generic adaptive vector RK23 utility.

Implementation details include:

- node packing and unpacking;
- velocity callback wrappers;
- source callback wrappers;
- transport stage lambdas;
- moment reconstruction handoff;
- diagnostic aggregation;
- timestep controller internals.

These details support the public APIs but should not become required concepts
for projects that use QMOMKit.

## Validation Expectations

Transport changes should be validated through analytical or manufactured
cases that exercise the relevant representation path. For a complete transport
macro step this is:

```text
moments
    -> QMOM/CQMOM inversion
    -> node advection
    -> raw-moment reconstruction
    -> moment source update
    -> post-source inversion
    -> final node advection/reconstruction
```

Existing verification covers:

- node advection against analytical continuous trajectories;
- moment-source updates against analytical moment evolution;
- complete QMOM transport for commuting and noncommuting operators;
- complete CQMOM transport for a compact two-dimensional case;
- Strang, `cd`, and `dc` scheme behavior;
- reconstruction error and zeroth-moment checks;
- diagnostics for successful transport runs.

Changes to splitting order, status mapping, callback validation,
realizability hooks, or diagnostic fields are transport-design changes and
require matching verification coverage.
