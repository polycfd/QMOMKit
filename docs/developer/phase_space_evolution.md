# Phase-space evolution

This page explains how a state advances by one step: the state itself, the
splitting of a step into node advection and the moment source, the
reconstruction after the source, how node counts change, and how failures are
reported. Public usage is in [Evolving in time](../user/5_evolving_in_time.md)
and [Coupling QMOMKit](../user/6_coupling_qmomkit.md); the reconstructions
called here are in [QMOM](qmom.md) and [CQMOM](cqmom.md).

## States

`QmomState` and `CqmomState` are immutable: they have no setters and no public
constructor. `initial_state()` and `state_candidate(...)` create them from
moments (`QmomCase`, `CqmomCase::initial_state()`, `CqmomPreparedCase`), and
the evolution stages create them from other states through the internal
access classes ([Architecture](architecture.md#public-and-internal-boundary)).

| | `QmomState` | `CqmomState` |
| --- | --- | --- |
| Moments | $m_0, \ldots, m_{K-1}$, as many as the case or candidate gave | the required moments of the hierarchy's layout |
| Reconstruction settings | maximum node count, `QmomOptions`, abscissa bounds, copied from the case | the settings of the `CqmomCase`, and the prepared branch structure, shared read-only along the lineage |
| Node counts | `active_node_count()`, the number of quadrature nodes | `active_branch_counts()`, one entry per branch; empty for a fixed case |
| Quadrature | `quadrature()`, the accepted `QmomResult` | `quadrature()`, the accepted `CqmomQuadrature` |

Creating a state reconstructs the quadrature once and then sums the stored
moments from it (`qmom_state_candidate` and `cqmom_state_candidate` in
[`src/public_api.cpp`](../../src/public_api.cpp)). QMOM reconstructs at the
case's node count; CQMOM starts every branch from its configured maximum. With
$n$ active QMOM nodes, the stored $m_0, \ldots, m_{2n-1}$ reproduce the given
moments to rounding and every higher one is that of the quadrature; $n$ is
below the case's node count when an adaptive reconstruction dropped nodes. For
CQMOM, an adaptive reduction likewise replaces the given required moments by
those of the reduced quadrature. Stored moments that overflow a double make
the candidate `reconstruction_failed`.

`initial_state()` runs the same construction as `state_candidate(...)` and
differs only in how it reports failure. It is called once, at setup, where a
failure means the case itself is wrong, so it throws; `state_candidate(...)`
runs for every location and step of a coupled host, where an outcome other than
`success`, such as an emptied cell, is an expected local event, so it returns a
status.

Moments and quadrature of a state therefore always agree, and a step starts
from the stored quadrature without reconstructing it. Every successful stage
builds its result state the same way, from the quadrature it accepted.

## Step structure

A step composes node advection and the moment source by operator splitting.
`execute_phase_space_evolution_sequence` in
[`src/phase_space_evolution.cpp`](../../src/phase_space_evolution.cpp) runs
the scheme:

| `PhaseSpaceEvolutionScheme` | Sequence over $[t, t+\Delta t]$ | Order |
| --- | --- | --- |
| `strang` (default) [[Strang, 1968](https://www.doi.org/10.1137/0705041)] | advection $\Delta t/2$ from $t$; source $\Delta t$ from $t$; advection $\Delta t/2$ from $t+\Delta t/2$ | second |
| `advsrc` | advection $\Delta t$, then source $\Delta t$, both from $t$ | first (Lie) |
| `srcadv` | source $\Delta t$, then advection $\Delta t$, both from $t$ | first (Lie) |

Strang splitting is second order for a smooth solution whose callbacks receive
time-dependent host data to second order at the stage times; how a host does
that is in
[Coupling QMOMKit](../user/6_coupling_qmomkit.md#host-data-that-changes-in-time).

The source always integrates over the whole interval from $t$, also in Strang
splitting, and evaluates its rates at the SSP-RK3 stage times $t$,
$t+\Delta t$ and $t+\Delta t/2$. A boundary contribution is added to the source
rate inside the source operator; it is not a third split operator.

A step is all-or-nothing. Each stage works on the state returned by the
previous one; if a stage fails, the result holds the input state and
`final_time` is $t$, otherwise $t+\Delta t$. The stage helpers
`advance_qmom_advection`, `advance_qmom_moment_source`,
`advance_cqmom_advection`, and `advance_cqmom_moment_source` call the same
stage functions, roll back to their own input, and ignore `options.scheme`.

`diagnostics.failed_stage` uses Strang's vocabulary:
`first_advection_half_step` and `first_moment_reconstruction`,
`moment_source_step`, `post_source_inversion`, `second_advection_half_step`
and `final_moment_reconstruction`. The advection of `advsrc` and of
`advance_*_advection` reports as first, that of `srcadv` as second.

A new splitting scheme adds a `PhaseSpaceEvolutionScheme` value, its case in
`validate_phase_space_evolution_scheme` and
`execute_phase_space_evolution_sequence`, and its `.qmom` value in
`parse_phase_space_evolution_scheme_value`
([`src/configuration.cpp`](../../src/configuration.cpp)). A scheme only
sequences the existing advection and source stages; reconstruction and state
construction stay in those stages, separately for QMOM and CQMOM.

## Node advection

Node advection integrates $\mathrm{d}\boldsymbol{\xi}_j/\mathrm{d}t =
\dot{\boldsymbol{\xi}}(t, \boldsymbol{\xi}_j)$ for every node $j$. Weights, terminal
paths, storage layout, and active counts are kept; the stage then sums the
stored moments from the moved nodes (`QmomStateAccess::from_quadrature`,
`internal::reconstruct_trusted_cqmom_moments`). Advection never reconstructs a
quadrature, in QMOM or CQMOM, so it cannot change a node count and runs no
realizability check. Abscissa bounds apply only where a quadrature is
reconstructed from moments; advection neither clips nor rejects a node that
leaves them.

`integrate_nodes` in [`src/node_advection.cpp`](../../src/node_advection.cpp)
packs the coordinates node-major, node $j$ at $[jD, (j+1)D)$ for $D$
dimensions. This is both the integrator's state vector and the layout of
`PhaseSpaceVelocity`, so one velocity call writes the derivative of all nodes
straight into the integrator's buffer. `make_phase_space_velocity` adapts a
`PerNodePhaseSpaceVelocity` by calling it once per node in ascending order and
validating each return, at the cost of one call and one returned vector per
node and evaluation.

| `NodeAdvectionScheme` | Method | Step control |
| --- | --- | --- |
| `rk23` (default) | Bogacki–Shampine 3(2), four stages [[Bogacki & Shampine, 1989](https://www.doi.org/10.1016/0893-9659(89)90079-7)] | embedded error estimate |
| `rk45` | Dormand–Prince 5(4)7M, seven stages [[Dormand & Prince, 1980](https://www.doi.org/10.1016/0771-050X(80)90013-3)] | embedded error estimate |
| `ssp_rk3` | SSP-RK3 [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)] | none: one step over the interval |

The embedded pairs run through `integrate_embedded_rk`, which also integrates
backwards when its end time lies before its start time. Phase-space evolution
always steps forward, since its step must be positive; the public
`advance_nodes` checks only that its step is finite, so a negative step moves
the nodes backwards. Node advection integrates in the local time
$\tau \in [0, \Delta t]$ and evaluates the velocity at $t+\tau$. The nodes then
move over exactly $\Delta t$, not over the rounded $(t+\Delta t)-t$, and the
step floor follows the spacing of $\tau$, not of $t$. `ssp_rk3` uses its own
three-stage formula, `integrate_nodes_ssp_rk3`, without the predicates and
retries of the source operator.

### Step-size control

Each trial step forms the higher-order solution $y$ and the embedded solution
$\hat{y}$ from the start value $y^0$, with the error norm

```math
E = \left[\frac{1}{n}\sum_{i=1}^{n}\left(\frac{y_i - \hat{y}_i}{a + r\max(|y^0_i|, |y_i|)}\right)^2\right]^{1/2},
```

$a$ and $r$ the absolute and relative tolerances. A step with $E < 1$ is
accepted and keeps $y$. The next step size is the last one times
$0.9\,E^{-1/(\hat{p}+1)}$, at most 10 after an accepted step (at most 1 if it
followed a rejection) and at least 0.2 after a rejected one; $\hat{p}$ is the
embedded order, so the exponent is $-1/3$ for `rk23` and $-1/5$ for `rk45`. Both
pairs are first-same-as-last: the last stage derivative of an accepted step is
the first of the next, and a rejected trial reuses the first derivative.

`AdaptiveRkOptions` in
[`include/qmomkit/runge_kutta.hpp`](../../include/qmomkit/runge_kutta.hpp)
holds the controls; the
[Reference](../user/7_reference.md#evolution-options) lists their defaults.
Step sizes are in the units of the time argument. The first trial step is never
estimated from the problem, rejected trials do not count toward `max_steps`,
and the step floor is the larger of `min_step` and ten spacings of the double
at the integration time.

The floor takes precedence over the tolerances. A step at the effective lower
bound, the smaller of the floor and the remaining interval, is accepted
whenever its solution and error are finite, so `success` means the interval
was covered, not that every step met the tolerance. `step_size_underflow`
reports a step at that bound whose solution or error estimate is not finite;
`max_steps_exceeded` reports that `max_steps` accepted steps did not cover the
interval. On either, node advection reports `integration_failed` and keeps
the input nodes. Invalid options, such as an explicit `min_step` above an
explicit `max_step`, are `NodeAdvectionStatus::invalid_input` before the
velocity is called.

## Moment source

The source operator advances the stored moment vector with SSP-RK3
[[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)],
`ssp_rk3_step` in [`src/runge_kutta.cpp`](../../src/runge_kutta.cpp):

```math
\begin{aligned}
y^{(1)} &= y^n + \Delta t\, f(t, y^n),\\
y^{(2)} &= \tfrac{3}{4} y^n + \tfrac{1}{4}\bigl(y^{(1)} + \Delta t\, f(t+\Delta t, y^{(1)})\bigr),\\
y^{n+1} &= \tfrac{1}{3} y^n + \tfrac{2}{3}\bigl(y^{(2)} + \Delta t\, f(t+\tfrac{1}{2}\Delta t, y^{(2)})\bigr).
\end{aligned}
```

It is explicit and has no error control; `SspRk3Options::estimate_error` is a
diagnostic only.

Acceptance predicates can reject a step
[[Nguyen et al., 2016](https://www.doi.org/10.1016/j.jcp.2016.08.017)]: a
realizability policy and the caller's `accept` and `stage_accept` in
`MomentSourceUpdateOptions::ssp_rk3_options`, combined by logical and.
`stage_accept` sees $y^{(1)}$, $y^{(2)}$, and $y^{n+1}$, and a rejected stage
ends the attempt before the next rate evaluation. After a rejection the step is
halved and restarted from $y^n$, at most `max_retries` times (default 0, so the
first rejection fails as `rejected_by_predicate`) and never below
`min_time_step` (default 0, no floor). A halving below `min_time_step`, or one
that underflows or no longer shrinks the step, fails as
`minimum_timestep_reached`. A requested step already shorter than
`min_time_step` is attempted, not refused.

A retried step covers only part of the interval. The driver
`advance_qmom_moments_with_source` in
[`src/moment_source.cpp`](../../src/moment_source.cpp) then calls
`ssp_rk3_step` again from the accepted state for the whole remaining interval,
with a fresh retry budget, until the interval is covered. An accepted substep
that does not shorten the remaining interval in floating point fails the update
as `source_failed`, so the loop always ends. A failed substep returns the input
moments, also after accepted substeps, and `diagnostics.retries` sums the
retries of all substeps.

`advance_qmom_moments_with_source` and `advance_cqmom_moments_with_source`
([`include/qmomkit/moment_source.hpp`](../../include/qmomkit/moment_source.hpp))
are this operator without a state: they advance a moment vector or
`CqmomMomentData` and reconstruct nothing. The stages call the
`internal::advance_trusted_qmom_moments_with_source` and
`internal::advance_trusted_cqmom_moments_with_source` overloads
([`src/internal/moment_source.hpp`](../../src/internal/moment_source.hpp)),
which give the source contexts the state's reconstruction settings.

## Source contexts

Each rate evaluation builds a `QmomSourceContext` or `CqmomSourceContext` on
the stage moments: a non-owning view, valid during the callback, whose
constructors from temporaries are deleted. Its accessors are listed in
[Evolving in time](../user/5_evolving_in_time.md#source).

A tracked moment comes from the stage values. An untracked moment, a
generalized moment, or a node visit uses a quadrature reconstructed from the
stage moments on first request and cached in the context, so at most once per
stage. It is never the state's quadrature, which belongs to the start of the
step. The reconstruction uses the state's settings and its active counts as
ceilings: QMOM calls `compute_qmom` with the state's active node count; CQMOM
passes `active_branch_counts()` to `internal::compute_cqmom_internal` with
growth disabled. A failed stage reconstruction throws inside the callback and
fails the source step.

Without a state, as in `advance_*_moments_with_source`, there are no
reconstruction settings: `moment(...)` works for tracked orders only, and
`moment(...)` for an untracked order, `generalized_moment(...)` and the node
visitors throw `std::out_of_range`, which fails the update.

## Node counts

Active node counts follow one policy:

- **Advection** keeps every count.
- **Source stages**: a stage reconstruction uses the state's counts as
  ceilings, so a branch may lose nodes but never grows. Only the branches
  below a branch that lost nodes, which condition on other nodes, start from
  their configured maximum.
- **QMOM post-source reconstruction**, `select_post_source_qmom`: a fixed case
  must return its configured count. An adaptive case reconstructs at the
  state's active count $n$ and accepts a reduction as is; if $n$ is kept and
  below the maximum, it tries $n+1$ and accepts that only when the
  reconstruction returns exactly $n+1$ nodes. With abscissa bounds,
  `internal::select_bounded_quadrature` applies the same rule and accepts the
  extra node only within the bounds.
- **CQMOM post-source reconstruction**, `internal::compute_cqmom_internal` with
  growth enabled: every branch follows the QMOM rule with its own active count;
  the mechanism, including the reset below a branch whose count changed, is in
  [CQMOM](cqmom.md#active-branch-counts).
- **Created states** start from the configured maximum ([States](#states)).

Counts therefore shrink readily and grow by at most one node per branch and
source step. `diagnostics.cqmom_node_count_transitions` records each branch
that was reduced, grew, was created, or was removed; a branch below a changed
ancestor is reported as removed and created even when its path and count are
unchanged (`record_cqmom_node_count_transitions`).

## Realizability policies during evolution

Policies act only as acceptance predicates of the source step; advection and
the post-source reconstruction run none. `FinalOnly` checks $y^{n+1}$,
`EveryStage` also $y^{(1)}$ and $y^{(2)}$, and a rejection follows the retry
rule above. Where policies are set and which moments they check is in
[Realizability](realizability.md#policies).

## Option validation

Every evolution entry point, `advance_*_phase_space`, `advance_*_advection`,
`advance_*_moment_source`, and their batches, checks the options of the stages
it runs before running any stage or calling any callback, and reports a
violation as `invalid_input` with `failed_stage == none`. A step failure
therefore always concerns the numerics or a callback. The checks are those the
stages apply themselves:

| Options | Validator | Checked |
| --- | --- | --- |
| `node_options.adaptive_options` | `internal::adaptive_rk_options_error` ([`src/internal/embedded_runge_kutta.hpp`](../../src/internal/embedded_runge_kutta.hpp)) | for `rk23` and `rk45` |
| QMOM `source_options` | `internal::validate_qmom_moment_source_options` ([`src/internal/moment_source.hpp`](../../src/internal/moment_source.hpp)) | when the source runs |
| QMOM `realizability_policy` | `internal::realizability_policy_error` ([`src/internal/realizability.hpp`](../../src/internal/realizability.hpp)) | when the source runs |
| CQMOM `source_options` | `internal::validate_cqmom_moment_source_options`, which applies `realizability_policy_error` per coordinate and rejects a policy beyond the dimension count | when the source runs |

Options a call does not read are not checked: advection helpers ignore source
options, source helpers ignore node options, and `ssp_rk3` ignores
`adaptive_options`. The stateless source operators run the same validators and
report `MomentSourceUpdateStatus::invalid_input`.

## Callbacks

| Callback | Shape | Contract |
| --- | --- | --- |
| `PhaseSpaceVelocity` | fills `velocities`, packed like `coordinates` | no resize; every entry written and finite |
| `QmomMomentSource`, `CqmomMomentSource`, boundary contributions | return a rate | one finite value per context entry, in context order |

The velocity fills a buffer because it runs at every stage of every
node-advection step, and filling keeps those calls free of allocation;
`OdeRhs`, the integrators' interface, fills one for the same reason. A source
runs once per SSP-RK3 stage, so returning a vector costs little and keeps the
most common callback simple.

Size and finiteness are checked after every call
(`internal::validate_packed_phase_space_velocity_result` in
[`src/internal/phase_space_velocity_validation.hpp`](../../src/internal/phase_space_velocity_validation.hpp),
`internal::rhs_checked` in
[`src/internal/ode_rhs_validation.hpp`](../../src/internal/ode_rhs_validation.hpp),
and the rate checks in [`src/moment_source.cpp`](../../src/moment_source.cpp))
and a violation throws `std::invalid_argument`. Node advection catches a
`std::exception` from the velocity or its checks as `velocity_failed`, the
source operator one from a source, boundary, acceptance predicate, or stage
reconstruction as `source_failed`, and the step maps them to `advection_failed`
and `moment_source_failed`, keeping the message. A non-finite or unwritten
velocity is reported with the time, dimension, and coordinates of its node.
Exceptions not derived from `std::exception` propagate. The generic
`ssp_rk3_step` and `integrate_embedded_rk` let callback exceptions propagate.
Both report invalid input, such as an empty or non-finite state, as
`SspRk3Status::invalid_input` or `AdaptiveRkStatus::invalid_input` with a
message before calling `rhs`, and throw `std::invalid_argument` when `rhs`
breaks its contract; `ssp_rk3_step` also rejects a step that is not positive and
throws when a candidate state is not finite.

Unwritten entries are detected in debug builds only. Without `NDEBUG`,
`internal::poison_output_buffer`
([`src/internal/numeric_utils.hpp`](../../src/internal/numeric_utils.hpp)) fills
the buffer with a NaN of a reserved bit pattern before the call, and an entry
holding that pattern afterwards is reported as unwritten, before the finiteness
check. Release builds cannot detect an unwritten entry; both node-advection
paths and `internal::rhs_checked` zero the buffer before each call, so such an
entry contributes a zero derivative.

## Statuses

| `PhaseSpaceEvolutionStatus` | Cause | Detail in `diagnostics` |
| --- | --- | --- |
| `invalid_input` | input or option check | `failed_stage == none` |
| `advection_failed` | velocity callback or node integration | `first_advection_status`, `second_advection_status`, `NodeAdvectionDiagnostics::adaptive_status` |
| `moment_reconstruction_failed` | moments of the moved or reconstructed quadrature not finite | `first_reconstruction_diagnostics`, `final_reconstruction_diagnostics` after advection; `failed_stage` and the message after the source |
| `moment_source_failed` | source or boundary callback, stage reconstruction, predicate rejection, minimum step | `moment_source_status`, `moment_source_diagnostics` |
| `inversion_failed` | post-source reconstruction | `post_source_qmom_status`, `post_source_cqmom_status`, `post_source_cqmom_diagnostics` |

`first_reconstruction_diagnostics` and `final_reconstruction_diagnostics`
describe the moment sums of the first and second advection stage, and are
written only after that stage's node integration succeeded; otherwise they keep
their defaults, `success == false` and `moment_count == 0`, also for a scheme or
helper without that stage. The moments summed after the post-source
reconstruction have no record of their own.

On failure, `CqmomPhaseSpaceEvolutionResult::diagnostic_quadrature` holds the
quadrature the failing stage produced or started from: the rejected
post-source candidate, the moved quadrature whose moments failed, or the stage
input. It is diagnostic only, and empty on success and for `invalid_input`.

## Batches

The `_batch` functions are serial loops (`apply_batch`) over the scalar
functions, with one time, step, set of options, and set of callbacks. They
return one result per state, in input order; a failed result is appended and
the loop continues. The `_batch_with_index` variants wrap each callback for
each state in a closure that passes `batch_index`, the position in the input
vector. An empty input returns an empty result without validation, and an
exception escaping a scalar call stops the loop.

## Verifying a change

A change to a stage, a scheme or an integrator needs the evidence of
[Testing](testing.md#verification-evidence), here from the exact solutions in
[`tests/verification/phase_space_evolution/`](../../tests/verification/phase_space_evolution/README.md):
commuting cases, whose split solution is exact, isolate the integrators, and
noncommuting cases measure the splitting order, about two for `strang` and one
for `advsrc` and `srcadv`. A new scheme adds its order there, and a new status
or diagnostic a component-regression test of the failing path.

## Where the code lives

| Function | File | Role |
| --- | --- | --- |
| `qmom_state_candidate`, `cqmom_state_candidate` | [`src/public_api.cpp`](../../src/public_api.cpp) | create states for `initial_state()` and `state_candidate(...)` |
| `execute_phase_space_evolution_sequence` | [`src/phase_space_evolution.cpp`](../../src/phase_space_evolution.cpp) | splitting schemes |
| `advance_qmom_advection_stage`, `advance_cqmom_advection_stage` | [`src/phase_space_evolution.cpp`](../../src/phase_space_evolution.cpp) | advection stage |
| `advance_qmom_moment_source_stage`, `advance_cqmom_moment_source_stage`, `select_post_source_qmom` | [`src/phase_space_evolution.cpp`](../../src/phase_space_evolution.cpp) | source stage and post-source reconstruction |
| `advance_nodes`, `integrate_nodes`, `make_phase_space_velocity` | [`src/node_advection.cpp`](../../src/node_advection.cpp) | node advection |
| `integrate_embedded_rk`, `ssp_rk3_step` | [`src/runge_kutta.cpp`](../../src/runge_kutta.cpp) | generic integrators |
| `take_embedded_rk_step`, `embedded_rk_tableau` | [`src/internal/embedded_runge_kutta.cpp`](../../src/internal/embedded_runge_kutta.cpp) | embedded stage arithmetic |
| `advance_qmom_moments_with_source`, `advance_cqmom_moments_with_source`, source contexts | [`src/moment_source.cpp`](../../src/moment_source.cpp) | source operator |
