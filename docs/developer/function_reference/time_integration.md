# Time integration and node motion

Public algorithms in [runge_kutta.cpp](../../../src/runge_kutta.cpp) are
declared in [runge_kutta.hpp](../../../include/qmomkit/runge_kutta.hpp).
Internal tableaus and stage arithmetic live in
[embedded_runge_kutta.cpp](../../../src/internal/embedded_runge_kutta.cpp) and
its [header](../../../src/internal/embedded_runge_kutta.hpp).
[node_advection.cpp](../../../src/node_advection.cpp) adapts them through
[node_advection.hpp](../../../include/qmomkit/node_advection.hpp). These
algorithms advance generic vectors; physical meanings reside in callbacks. See
[phase-space design](../phase_space_evolution_design.md).

## SSP-RK3 attempt and retries — Level A

Local `try_ssp_rk3_step` performs exactly one attempt from the supplied
original state `y` and positive step `h`:

1. Evaluate `k1=f(t,y)`, form `y1=y+h*k1`, then apply the Stage1 predicate.
2. Evaluate `k2=f(t+h,y1)`, form `y2=3*y/4+(y1+h*k2)/4`, then apply Stage2.
3. Evaluate `k3=f(t+h/2,y2)`, form `yn=(y+2*(y2+h*k3))/3`, require finite output,
   then apply the Final stage predicate and final acceptance predicate.

Each vector loop combines the same component across stages. An early rejected
stage prevents later RHS calls. One derivative buffer serves all three stages,
since each stage's derivative is consumed into the next stage state before the
following stage overwrites it. `rhs_checked` sizes that buffer, hands it to the
callback, and enforces that the callback neither resized it nor left an entry
unwritten, and that every derivative is finite. Intermediate state acceptance
depends on predicates and subsequent RHS checks; there is an explicit finite
check on the final candidate.

The second stage also constructs the second-order test state `(y+y1+h*k2)/2`.
With `estimate_error`, the reported estimate is the Euclidean norm of its
difference from `yn`, divided by `||yn||` unless that norm is zero. It is
**not** an error-controlled acceptance criterion. Norm computation can itself
overflow at extreme scale; the estimate is not a certified error bound.

Public `ssp_rk3_step` validates initial state/time/step and a finite
nonnegative minimum step, then loops over attempts. A step already shorter
than the minimum is attempted, just not halved. Predicate rejection halves
`h` and restarts at the **original** `t,y`; it does not continue a rejected
intermediate state. `max_retries` or the configured minimum stops retry with
original state and a rejection/minimum-step status. Arithmetic or callback
exceptions are not handled as predicate retries. Source wrappers translate
those exceptions into source failure. Before a retry is evaluated, its halved
step must be finite, positive, and strictly smaller than the rejected step;
otherwise the stepper returns `minimum_timestep_reached`. A zero configured
minimum remains valid and means no positive configured floor.

Evidence: [RK unit tests](../../../tests/unit/runge_kutta.cpp) inspect the
stage formula/times, predicate short-circuit, retries, and error-estimate
behavior; [source tests](../../../tests/unit/moment_source.cpp) establish that
the enclosing source driver finishes the whole interval. These are separate
acceptance layers.

## Embedded RK stage arithmetic — Level A

Internal `embedded_rk_tableau` selects the implemented Bogacki–Shampine RK23
(4 stages) or Dormand–Prince RK45 (7 stages); unknown methods throw.
`take_embedded_rk_step` loops over tableau stages, prior-stage derivatives, and
vector components to construct each stage state. It then combines the same
derivatives with higher-order and embedded weights, leaving both candidate
states and their componentwise difference in the caller's
`EmbeddedRkWorkspace`. The tableaux are the source of the coefficients; no
alternate backend is selected.

The workspace is where the stage vectors, the two candidate states, the error
vector, and the cached first derivative live. The controller owns one for the
whole integration and passes it to every attempt, so these buffers are
allocated on the first attempt instead of on each one; it is a controller
local, not shared state, so concurrent integrations do not interfere. Every
buffer is assigned its full length before it is read, which is what makes
reuse arithmetic-neutral.

A workspace whose `has_first_derivative` is set supplies the first derivative;
otherwise it is evaluated at the starting state. On return the flag is set,
because slot zero of the stage derivatives is then the slope at the starting
state whichever way it was obtained, so a rejected attempt's retry reuses it
without the controller moving it anywhere. An accepting controller must instead
install the outgoing slope of a first-same-as-last (FSAL) tableau — its last
stage, swapped into slot zero — or clear the flag. RHS shape and finite
derivative checks remain in this layer. Cache validity is an internal caller
invariant: the cached derivative must belong to the current starting
state/time. `embedded_rk_controller_exponent` returns `-1/(embedded_order+1)`
for the error controller.

Evidence: [RK tests](../../../tests/unit/runge_kutta.cpp) check method-specific
stage arithmetic and function-evaluation accounting, including FSAL reuse, and
that a reused workspace takes the same step as a fresh one while its cache
replaces exactly one evaluation. No independent arbitrary-tableau facility or
interpolation/dense-output contract is implemented by these helpers.

## Embedded controller — Level A

`integrate_embedded_rk` advances an initial finite vector from `t0` to
`t_bound`, including backward integration. It resolves the tableau, validates
dimensions, times, positive finite relative/absolute tolerances, and nonzero
`max_steps`. Equal endpoints return success before checking nonessential
step-size controls. For a nonzero interval, initial step must be positive; zero
maximum step means interval length, and zero minimum means ten floating-point
spacings at current time in the integration direction. An explicit minimum is
rejected only above an explicit maximum; an interval shorter than the minimum
is taken as one step, accepted if its error estimate is finite.

The outer loop tracks accepted state/time until the endpoint or accepted-step
limit. The inner loop retries one proposed step:

- Clamp the proposed magnitude to configured limits, then clip it to the endpoint.
  The effective lower bound is the smaller of minimum step and remaining interval.
- Compute an embedded pair. `error_norm` is the RMS of component errors scaled by
  `atol + rtol*max(abs(y_i),abs(yn_i))`.
- Accept when the finite norm is below one **or when the finite candidate and
  finite norm are at the effective lower bound**, even if the norm exceeds one.
- On acceptance, advance state/time and cache the endpoint derivative. Scale the
  next step using safety `0.9`, upper factor `10`, and the controller exponent;
  after a rejection, suppress immediate growth above factor one.
- On rejection, retain the starting derivative and reduce by at least the
  controller's factor floor `0.2`. A nonfinite candidate/error at the lower bound
  returns step-size-underflow status. The accepted-step budget counts accepted
  steps, not rejected attempts.

Failure statuses retain the last accepted state/time. Invalid RHS values may
throw instead of returning an integration status. Lower-bound acceptance is an
explicit implemented policy. It means success alone does not certify
satisfaction of the requested local-error norm. Counters and last step are
diagnostic, not a global error estimate.

Evidence: [RK tests](../../../tests/unit/runge_kutta.cpp), notably
`test_rk23_accepts_finite_update_at_minimum_step`, plus short-final-interval
and RK45 cases, test the policy directly. Autonomous analytical trajectories
and refinement tests provide accuracy evidence for specific cases; they do not
characterize floor-dominated integration or all time/scale combinations.

## Node trajectory adaptation — Level A

Local `pack_nodes` stores coordinates node-major in one ODE vector; weights are
excluded. `node_rhs_from_velocity` hands that vector to the velocity callback
unchanged at every stage — the packed callback layout is the ODE state layout,
so no gather or scatter is needed — and validates the filled rank/finiteness,
attributing any failure back to its node and phase-space dimension.
`unpack_nodes` restores per-coordinate arrays. This preserves terminal list
order and does not reconstruct a new CQMOM tree.

`integrate_nodes_embedded_rk` integrates over local time `[0, h]` with the
velocity evaluated at `t + tau`, so the nodes move over `h` itself rather than
the rounded interval `(t+h)-t`, and reports `t` plus the local final time. The
step floor then follows the spacing of `tau`, so the error control still
resolves a velocity that is stiff in the node positions in a step far shorter
than the spacing at `t`; its variation in time within that spacing is not
resolved, since `t + tau` rounds to it. It
translates controller failure into integration failure and returns original
coordinates on failure.
Callback exceptions become velocity failure. `integrate_nodes_ssp_rk3` uses the
same three-stage trajectory formula directly, with signed or zero step
supported; it has no source predicates, source retry loop, or embedded
controller. Thus selecting SSP-RK3 for advection is different from invoking the
source stepper. `integrate_nodes` dispatches the chosen scheme.

Public `advance_nodes` overloads validate finite shape, weights, coordinates,
and where applicable producing QMOM success status or CQMOM terminal paths.
Successful motion preserves weights, hierarchy, terminal paths, and
rectangular/ragged flags; only coordinates change. It does not constrain
intermediate trajectories to configured support intervals. Support-aware
reconstruction is a different stage. The caller's velocity need not conserve
kinetic energy or any physical quantity.

Evidence: [node-advection unit tests](../../../tests/unit/node_advection.cpp)
and [component tests](../../../tests/component_regression/node_advection.cpp)
cover constant/linear motion, schemes, callback failure, and retained
weights/topology. [Phase-space
tests](../../../tests/unit/phase_space_evolution.cpp) verify later moment
reconstruction and stage composition. Boundary crossings during a node
trajectory remain a separate question from bounded inversion acceptance.

## Supporting functions — Level C

| Family | Responsibility / important semantics | Evidence |
|---|---|---|
| `require_same_size`, `rhs_checked`, `add_scaled` | Componentwise shape-checked arithmetic and derivative checks; callback exceptions propagate to the owning wrapper. | RHS error tests |
| `euclidean_norm`, `error_norm` | Ordinary double norm reductions; different SSP diagnostic and embedded acceptance meanings. | RK formula/controller tests |
| `sign`, `validate_or_default_max_step`, `minimum_step` | Direction, zero-as-default maximum, representable-time lower bound. | Backward/step-boundary tests; `embedded_rk_resolves_a_stiff_velocity_far_from_the_origin` |
| `stage_accepted`, `invalid_embedded_rk_result` | Missing predicate means acceptance; invalid controller input preserves initial state/time. | Invalid/predicate tests |
| Node `all_finite`, `validate_weights`, `validate_abscissas` | Shape and finite input checks, independent of physical coordinate meaning. | Node invalid-input tests |
| `diagnostics_from_embedded_rk`, `scheme_name`, `format_node_advection_options` | Transfer integration counts/time and serialize relevant scheme options at explicit precision; do not assess accuracy. | Formatting/diagnostic tests |

Exact definitions and overload locations are in the [inventory](coverage_inventory.md).
