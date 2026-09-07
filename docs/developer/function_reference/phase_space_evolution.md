# Phase-space evolution, state transitions, and batches

Definitions are in [phase_space_evolution.cpp](../../../src/phase_space_evolution.cpp),
with public stage/macro/batch declarations in
[phase_space_evolution.hpp](../../../include/qmomkit/phase_space_evolution.hpp).
Local orchestration is in anonymous namespaces; public `advance_*` functions belong
to `qmomkit`. See [design](../phase_space_evolution_design.md),
[sources](moment_sources.md), and [node integration](time_integration.md).

## Operator order and clocks — Level A

Local template `execute_phase_space_evolution_sequence` dispatches three schemes.
`A(start,length)` denotes node advection followed by moment reconstruction;
`S(start,h)` denotes a full source interval followed by quadrature reconstruction.

| Scheme | Actual invocation sequence for a requested interval `(t,h)` |
|---|---|
| `strang` | `A(t,h/2)`, `S(t,h)`, `A(t+h/2,h/2)` |
| `advsrc` | `A(t,h)`, `S(t,h)` |
| `srcadv` | `S(t,h)`, `A(t,h)` |

Each stage must succeed before the next is called. The source lambda captures the
full macro step; it does not infer its length from adjacent advection times.
The final-advection flag selects diagnostic placement, not a new integrator.
Unknown scheme is invalid. Macro input validation requires finite positive `h`,
finite `t` and `t+h`, valid numerical state, and required callbacks.

For autonomous callbacks these are recognizable operator compositions. For
explicitly time-dependent callbacks, the source stage retains the full physical
interval `[t,t+h]`, despite appearing between the two advection half-steps.
SSP-RK3 samples that source at `t`, `t+h`, and `t+h/2`.
Do not infer nonautonomous order solely from autonomous convergence tests.

Evidence: [evolution unit tests](../../../tests/unit/phase_space_evolution.cpp)
inspect operator calls, failures, and stage-local contexts;
[analytical verification](../../../tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp)
checks autonomous commuting/noncommuting examples and refinement behavior. A
source-only linear-time unit case checks the nonautonomous source clock.

## QMOM accepted stages — Level A

`advance_qmom_advection_stage` starts from the state's accepted quadrature, advances
nodes, and reconstructs the stored raw sequence at its existing length. The accepted
node count and weights are retained. It does not first invert the same moments
again. Numerical failure returns the supplied rollback state and failed-stage
diagnostics.

`advance_qmom_moment_source_stage` supplies current stage moments plus
reconstruction metadata to the source integrator, combines configured realizability
and user predicates, and advances the full source interval. It then calls local
`select_post_source_qmom`:

- Fixed mode requests configured maximum exactly.
- Adaptive mode first requests the current active count. Reduction accepts the
  reduced result without an immediate growth attempt.
- If current count survives unchanged below capacity, try exactly one more node.
  Failed/reduced growth preserves the prior successful candidate.
- Bounded selection and automatic scaling wrap the same acceptance responsibilities;
  neither node clipping nor silent tolerance relaxation is performed.

The accepted post-source quadrature becomes authoritative: all stored moments are
reconstructed from it through `QmomStateAccess`. In adaptive reduction, this can
replace higher evolved moments with the lower-node closure. A finite source vector
is insufficient for success if post-source inversion fails. Configured maximum,
active count, and attempted growth count remain distinct.

Evidence: [evolution tests](../../../tests/unit/phase_space_evolution.cpp) cover
reduction, one-node growth, failed growth, source-stage quadrature caching,
bounded reconstruction, and original-state rollback. Existing tests intentionally
exercise canonicalization; it must not be mistaken for exact preservation of all
higher ODE moments after reduction.

## Preparing CQMOM working state — Level A/B

`prepare_cqmom_phase_space_evolution_working_state` reuses the state's accepted
quadrature and active-count vector when present. Normal `CqmomCase::initial_state`
already populates them and canonicalizes dense or lean stored moments from that
same accepted quadrature. If absent, it reconstructs once using existing active
constraints, with growth disabled, derives counts from the resulting quadrature,
and records initial transitions. Failure ends the stage with rollback state.
No fresh reconstruction is required merely because this is the first public
advection call.

The working structure borrows cached state data or owns newly reconstructed data;
it is local to the stage. Fixed hierarchy/layout metadata is immutable and shared
along accepted state lineage. Dense and lean storage are dispatch choices, not
different operator-splitting algorithms.

## CQMOM accepted advection and source — Level A

`advance_cqmom_advection_stage` moves the accepted terminal coordinates with fixed
weights, hierarchy, paths, and active topology. Dense states regenerate every
stored tensor moment; lean states regenerate the required layout. Successors keep
the accepted moved quadrature and active counts. No source growth/reduction policy
runs during node motion.

`advance_cqmom_moment_source_stage` first extracts the required layout from a dense
state or directly uses the lean required data. Both representations integrate the
**same required moment vector**. Dense-only higher moments are not independently
advanced by the callback's tracked vector. Context closure supplies untracked
quantities from each RK stage as needed.

After source integration, reconstruct with previous active counts and configured
adaptive growth permission. Failure is an inversion failure, distinct from a
callback/predicate failure. On success, derive active counts and rebuild stored
moments from the accepted quadrature: all dense extents for dense storage, required
indices for lean storage. `dense_state_from_accepted_quadrature` and
`lean_state_from_accepted_quadrature` retain policies, hierarchy, and compatible
prepared metadata while replacing evolving values/cache.

`record_cqmom_node_count_transitions` compares before/after branch maps keyed by
level, conditioning parent path, and physical dimension. Its loops match records,
identify changed ancestors, and classify growth, reduction, creation, or removal.
When a strict ancestor changed, descendants are removed/created even if their
textual path/count happens to match. This records conditional-distribution identity,
not just numerical equality of integers. Recursive descendant constraints are
invalidated as described in [CQMOM reconstruction](cqmom_reconstruction.md).

Evidence: [evolution tests](../../../tests/unit/phase_space_evolution.cpp) include
`test_public_cqmom_advection_stage_preserves_topology`,
`test_public_cqmom_moment_source_stage_growth_is_limited_to_one_node`,
`test_lean_cqmom_moment_source_stage_canonicalizes_adaptive_reduction`, and
`test_dense_state_reconstructs_dense_only_moments_after_moment_source_source`.
Ancestor-collapse/growth tests check transition identity. These finite examples do
not characterize topology chatter or accuracy across many collapse/reappearance
cycles.

## Transaction and result materialization — Level B

`advance_qmom_phase_space_impl` and `advance_cqmom_phase_space_state_impl` compose
stages against a macro-step rollback state. Public standalone advection/source
calls instead roll back to their own stage input. A failure may include a partial
quadrature and a failed-stage time for diagnosis while the returned state is the
original state. Do not interpret that diagnostic time as proof the returned state
has advanced there.

`materialize_cqmom_phase_space_evolution_result` obtains the accepted quadrature
from a successful state's cache; absence is an internal logic error. Failure can
materialize a diagnostic quadrature separately. `neutral_cqmom_message`,
`neutral_cqmom_diagnostics`, and `neutral_cqmom_result` translate representation-
specific wording for neutral public APIs; they do not change numerical status.
`qmom_from_quadrature`, `qmom_failure`, `cqmom_state_failure`, and `cqmom_success`
package already decided results. Diagnostic record/finalization helpers preserve
per-stage counters and aggregate work, not a rigorous global error bound.

Evidence: failure and partial-quadrature assertions in
[evolution unit tests](../../../tests/unit/phase_space_evolution.cpp),
[component tests](../../../tests/component_regression/phase_space_evolution.cpp),
and analytical verification distinguish accepted state from diagnostic data.

## Batch and public forwarding — Level B/C

Local template `apply_batch` loops serially over input states in index order and
appends one result per invocation. Public advection, moment-source, and macro-step
batch overloads share callbacks/options or wrap callbacks with the current batch
index (`*_batch_with_index`). Boundary-enabled overloads retain source-plus-boundary
semantics. The helper does not schedule threads, combine states, or stop merely
because one result has a failure status. An uncaught exception can still interrupt
the loop. Empty batches perform no scalar calls.

Scalar public wrappers select the stage/macro helper and boundary composition;
they do not introduce additional numerical solvers. Validation helpers separate
common time/callback checks from dense/lean state requirements. `active_branch_key`,
key `operator<`, `path_has_prefix`, and `has_changed_strict_ancestor` supply the
transition map operations described above.

Evidence: [batch tests](../../../tests/unit/phase_space_evolution_batch.cpp)
compare scalar and batch results, indexed callbacks, and failure isolation at the
returned-status level. The [inventory](coverage_inventory.md) identifies every
forwarding overload and local validation/diagnostic definition.
