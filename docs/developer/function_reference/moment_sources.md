# Source integration and stage-local contexts

Definitions: [moment_source.cpp](../../../src/moment_source.cpp); public
contracts: [moment_source.hpp](../../../include/qmomkit/moment_source.hpp);
enriched internal entry points:
[internal/moment_source.hpp](../../../src/internal/moment_source.hpp). The
source operator advances raw moments using application callbacks. It does not
own a physical model. See [phase-space
design](../phase_space_evolution_design.md) and [SSP-RK3](time_integration.md).

## Completing a source interval — Level A

The public `advance_qmom_moments_with_source` overloads validate initial moments,
time, positive requested step, callbacks, and source controls, then adapt the
source to an `OdeRhs`. `QmomMomentSource` returns its rates and `OdeRhs` fills a
buffer, so the adapter copies the returned rates into the buffer the integrator
supplies; a source that returned the wrong number of rates resizes it there and
is reported as a wrong-sized right-hand side. The local overload with that same
name drives SSP-RK3. Its state is `state`, `current_time`, and
`remaining_time`.

On each iteration it asks `ssp_rk3_step` to advance the **whole remaining
interval**. A predicate rejection can make that step retry with a smaller
accepted step. Before committing an accepted step, the driver requires its
actual binary64 subtraction from the remaining interval to decrease that
interval. Thus a successful source operation completes the caller's interval
rather than silently returning after a halved step. Retry counts accumulate and
a rejection flag is retained; detailed error/step fields come from the most
recent substep.

Any unsuccessful substep returns the original full-interval moments, even if
previous substeps were accepted. Caught standard exceptions become source
failure. The loop has no separate cap on accepted substeps, but rejects an
accepted SSP step that is non-finite, non-positive, or fails to decrease the
remaining interval in actual floating-point arithmetic. Source controls must be
finite. The stepper's relative estimate is advisory and describes the final
accepted substep, not a global error bound for the full source interval; user
predicates control retry.

At large absolute physical times, distinct mathematical stage times can round
to the same binary64 value. Callbacks that depend explicitly on time must
therefore not assume that every nominal RK stage has a distinct representable
timestamp.

Evidence: [source unit tests](../../../tests/unit/moment_source.cpp) check
completion of the full requested interval after retries and rollback on later
failure; [component source
tests](../../../tests/component_regression/moment_source.cpp) cover integration
behavior. They do not define an envelope for extremely small accepted steps
relative to remaining time.

## CQMOM source vectors and marginal policy — Level A

`advance_cqmom_moments_with_source` validates canonical required data and
controls. Internal `advance_trusted_cqmom_moments_with_source` overloads reuse
metadata owned by a valid state. Local
`advance_validated_cqmom_moments_with_source` holds the layout fixed throughout
all RK stages and presents each stage's value vector in that layout to the
callback. Stage values and returned rates must retain matching size and finite
entries; a source cannot change hierarchy or moment indexing.

`extract_coordinate_marginal`, in `internal/cqmom_case.hpp` because
`CqmomPreparedCase::state_candidate` checks the same marginals, finds indices
whose other exponents are zero, sorts by the selected coordinate order, and
gathers their values. Upstream
canonical layout construction supplies the needed contiguous orders. It does
not marginalize by summing unrelated stored mixed moments.
`accepts_coordinate_realizability_policies` checks only coordinates with
configured policies, using their specified support and tolerance. These
marginal conditions are **not a mixed-moment CQMOM realizability theorem**.
User acceptance and configured realization predicates are combined;
`EveryStage` and `FinalOnly` determine when they are applied.

Failure preserves the original required data. Successful source-vector
integration alone does not imply successful multidimensional reconstruction;
the phase-space source stage subsequently reconstructs and canonicalizes an
accepted state.

Evidence: [source tests](../../../tests/unit/moment_source.cpp) check
coordinate selection, boundary policies, shape errors, trusted/ordinary
behavior, and callback failures. [Evolution
tests](../../../tests/unit/phase_space_evolution.cpp) cover post-source
inversion failures separately. Coordinate-wise acceptance does not establish
mixed-moment CQMOM realizability.

## Stage-local moment closure — Level A

`QmomSourceContext::moment` first returns a tracked raw moment directly. For an
untracked order, an enriched context reconstructs from **the current RK stage's
moments** using its active count, options, and support. It lazily caches that
stage quadrature, reconstructs orders through the requested order, and returns
the requested value. Without reconstruction metadata, untracked access throws.
A failed inversion throws; source wrappers translate the callback failure
rather than accepting a diagnostic quadrature.

`QmomSourceContext::generalized_moment` evaluates a finite non-negative scalar
order from that same cached current-stage quadrature through the shared
`evaluate_generalized_moment` implementation. It therefore preserves the
integer, zero-order, fractional-power, negative-abscissa domain, and overflow
semantics of the public quadrature evaluator. Generic source-only contexts
without reconstruction metadata reject the request.

`QmomSourceContext::for_each_node` rejects an empty callback, then passes each
weight and abscissa of that same cached current-stage quadrature to the
visitor, in quadrature order. Contexts without reconstruction metadata reject it
with the `std::out_of_range` that `generalized_moment` also throws. The weights
are the quadrature's absolute weights, so they carry the total number density,
and an adaptive reconstruction passes only its active nodes.

`CqmomSourceContext::moment` similarly uses exact required-index lookup first.
`reconstructed_stage_quadrature` constructs a stage `CqmomCase` and reuses the
prepared metadata, with the active counts of the state as ceilings: an adaptive
branch of the state may lose nodes but gets no growth trial, as the QMOM context
keeps its active node count, and grows only in the post-source inversion; a
branch that replaces others after its parent lost nodes starts from its
configured maximum, as in the post-source inversion. The cache belongs to this
context and therefore this RK stage; it is not a global or whole-time-step
cache. Missing stored moments are computed from that accepted quadrature, not
added to the stage's tracked ODE vector.

`CqmomSourceContext::generalized_moment` additionally requires finite
**nonnegative** exponents of matching dimension before generalized quadrature
evaluation. `for_each_terminal_node` traverses terminal weights/coordinates and
supplies each node to the visitor. Coordinates are a temporary vector;
callbacks must not retain references as persistent state. Global terminal
weights are already products of conditional weights; callers must not multiply
by another density implicitly.

These operations provide a finite quadrature closure. They do not recover
missing information uniquely from a truncated moment sequence, and an adaptive
reduction can change untracked higher moments. A stage may fail closure even if
its stored entries are finite. Fractional powers retain the negative-abscissa
domain restriction explained in [moment evaluation](moment_evaluation.md).

Evidence: [phase-space unit
tests](../../../tests/unit/phase_space_evolution.cpp) include
`test_qmom_source_context_uses_current_rk_stage_and_caches_quadrature`,
`test_qmom_source_context_generalized_moment_semantics`,
`test_qmom_source_context_generalized_moment_uses_current_rk_stage`,
`test_qmom_source_context_visits_nodes_from_each_rk_stage`,
`test_cqmom_source_context_reconstructs_untracked_moments_from_each_rk_stage`,
and `test_cqmom_source_context_generalized_moments_use_stage_quadrature`.
[Source unit tests](../../../tests/unit/moment_source.cpp) exercise standalone
metadata absence, a missing node-visitor callback, and indexed access.
`verification_qmom_node_source_matches_exact_translation` checks a growth
source written only through the node visitor against the exact translation of
the moments. These distinguish stage-local closure from reusing the pre-step
quadrature; they do not measure closure error for arbitrary physical kernels.

## Adding boundary contributions — Level A/B

`compose_source_and_boundary` and `compose_cqmom_source_and_boundary` return
one RHS that calls source and boundary at the same stage time/context,
validates each vector, then adds them componentwise. Invalid boundary
shape/value or callback failure aborts that RHS. The generic stepper checks the
combined derivative's finiteness, including overflow of two finite
contributions. Enriched QMOM/CQMOM trusted overloads preserve reconstruction
metadata for both callbacks.

This is addition **within** the source operator, not a third split operator.
Consequently source/boundary interaction uses the same RK stage states and
realizability predicates. No spatial mesh, inflow law, or physical units enter
these generic routines. Existing source and phase-space unit/component tests
exercise additive rates and failure propagation, not application-specific
balance laws.

## Context, validation, and result helpers — Levels B/C

| Functions/family | Responsibility and contract | Evidence |
|---|---|---|
| Context constructors; internal context `create` factories | Borrow stage values/layout and optionally attach reconstruction metadata. The caller owns referenced state for the callback lifetime. | Stage-cache tests |
| Context `time`, `moment_count`, `dimension_count`, `size`, `indices`, `values`, `index`, `value`, `position`, `hierarchy` | Expose stage metadata/tracked data; checked single-entry access throws for missing/out-of-range entries. `moment` is the distinct closure operation. | Source context tests |
| `validate_request`, `validate_moment_data_request`, `validate_trusted_moment_data_request`, `validate_cqmom_source_update_controls` | Separate public structural validation from trusted layout reuse; verify time input and callbacks before integration. | Invalid-input tests |
| `validate_qmom_moment_source_options`, `validate_cqmom_moment_source_options` | Check the source-step options: a finite non-negative `min_time_step`, no coordinate policy for QMOM, and for CQMOM no policy beyond the dimensions and every applied policy valid. Shared with phase-space evolution, which runs them before any stage. | Invalid-option tests |
| `validate_source_context`, `validate_fixed_cqmom_moment_data`, `validate_cqmom_stage_state`, `validate_rate`, `validate_source_rate`, `validate_boundary_rate`, `validate_cqmom_boundary_rate` | Validate the specific borrowed state or callback output boundary; no numerical repair. | Shape/nonfinite failure cases |
| `has_coordinate_realizability_policy` | Bypass conversion when no policy applies. | Marginal-policy tests |
| `status_from_ssp_rk3`, `message_from_status`, `diagnostics_from_ssp_rk3`, invalid/failed result factories, callback-error builders | Translate step failure and diagnostic fields while retaining original state on failure. | Retry/rollback/status tests |

All overload locations are in the [inventory](coverage_inventory.md).
