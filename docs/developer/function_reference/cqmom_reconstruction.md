# Recursive CQMOM reconstruction

The implementation performs successive one-dimensional inversions and
conditional linear solves in **conditioning order**, while public moment
indices and coordinates remain in **phase-space dimension order**. See
[recursive design](../cqmom_recursive.md) for the moment equations and
[conditional systems](conditional_systems.md) for the weighted solve.

Definitions are in
[recursive_quadrature.cpp](../../../src/recursive_quadrature.cpp), with
declarations in [the internal
header](../../../src/internal/recursive_quadrature.hpp) and
[cqmom.hpp](../../../include/qmomkit/cqmom.hpp). Internal symbols below belong
to `qmomkit::internal`; local helpers are in anonymous namespaces. Public
layout functions belong to `qmomkit`. Definition locations and overloads are in
the [inventory](coverage_inventory.md).

## Required moments and branch layouts — Level A

`create_cqmom_moment_layout` maps a `PhaseSpaceHierarchy` to the fixed required
set of raw moment indices. Its outer loop visits conditioning levels. At level
`l`, earlier conditioning axes range over their configured node counts, the
current axis over twice its configured count, and later axes only over exponent
zero. The inner enumeration converts ordered indices back to phase-space order,
retaining first occurrences and suppressing duplicates. These are
capacity-based requirements even when adaptive reconstruction ultimately uses
fewer nodes. `extract_cqmom_moments` gathers that set from a sufficiently large
full tensor; `create_cqmom_moment_data` checks values against the layout.
Missing required moments are errors, not implicitly zero values.

`create_ordered_root_layout` changes index interpretation to conditioning
order. `create_branch_moment_layout` constructs immutable lookup metadata for
one conditional slice: offsets for the pure leading-coordinate sequence and
layouts for its possible children. A one-dimensional slice terminates
recursion. Fixed metadata needs only the configured leading count; adaptive
metadata enumerates each possible count from one through that maximum. This
trades metadata work for reuse during numerical branch traversal; it is not
numerical quadrature caching.

`create_child_layout_case` removes the leading coordinate, retains unique
trailing indices whose leading order is below the candidate parent count, and
maps every child index back to parent orders `0..n-1`. Its nested loops
describe the RHS columns and parent-moment rows of the weighted Vandermonde
equation. It then recurses over the remaining extents.
`extract_wheeler_moments` gathers the leading raw sequence using cached
offsets; `make_branch_rhs_matrix` gathers the matrix with parent order as row
and trailing multi-index as column. Size disagreement or an unavailable child
layout throws rather than fabricating a conditional moment.

`PreparedCqmomMetadata` and `create_prepared_cqmom_metadata` pair
the canonical layout with its branch layout. Compatibility checks compare
hierarchy counts, conditioning order, required index order, and the adaptive
branch-count flag. Normal case construction prepares this once and shares
immutable metadata along state lineage.

Adaptive cases are prepared too. `create_branch_moment_layout` already
materializes a child case for every parent node count Wheeler could realize, so
an adaptive tree covers every runtime outcome and is a pure function of the
hierarchy, exactly as a fixed one is. It is wider, not count-dependent. The
flag is part of the compatibility check because the two trees differ: a fixed
tree carries only the configured count, so reusing one for an adaptive request
would leave the realized count without a child case. A case-derived
reconstruction that fails to find compatible metadata throws rather than
rebuilding, so a cache miss is loud. The trusted numerical path assumes this
ownership boundary has already established compatibility; it is not a public
bypass for arbitrary layouts.

Evidence: [recursive unit tests](../../../tests/unit/recursive_quadrature.cpp)
exercise required offsets, non-identity order, unequal node counts, and reduced
child layouts; [configuration tests](../../../tests/unit/configuration.cpp) and
[public API tests](../../../tests/unit/public_api.cpp) inspect metadata
ownership/reuse and tensor-to-required extraction. Preparing all adaptive child
layouts costs 2^d - 1 layouts at two nodes per axis, 63 at six dimensions, and
121 for five dimensions at three nodes; peak resident size of the five-
dimensional adaptive example moves from 7.0 MB to 7.2 MB. Growth beyond those
ranks and counts is not characterized.

## Conditional slices — Level A

Local `reconstruct_child_slices` receives the current branch's mixed moments
and the accepted parent quadrature. It builds the child layout for the
**returned** parent count, then solves all trailing moment columns together
through `solve_conditional_vandermonde_system`. Each solution row becomes the
raw moment sequence of one conditional child distribution; global branch
probability is carried separately in the tree's accumulated weight.

The function validates shapes, parent values, and the presence of another
coordinate. A failed solve stops the branch. Its solver request disables
residual acceptance checking: detailed diagnostics may calculate
residual/condition information, but enabling diagnostics does not introduce a
numerical acceptance test. The outer reconstruction decides whether to retain
successful branch records. This differs from the public standalone
conditional-moment API, whose solver route is described in [conditional
systems](conditional_systems.md).

Evidence: the solver-policy and diagnostic-level cases in [recursive unit
tests](../../../tests/unit/recursive_quadrature.cpp) compare moments and status
across dense/BP/QR/SVD/Moore-Penrose paths. Agreement on these cases is not an
independent conditioning limit or evidence that every finite solution is
accurate.

## Active topology and branch limits — Level A

The optional prior active-count records are used sorted by parent path, the
order in which a reconstruction records them; records in another order are
sorted into a copy first. In that order every branch follows its parent, so
`validate_active_branch_counts` checks them in one pass without building a
lookup structure: level, parent-path length, physical dimension, configured
maximum, and positive active count bounded by capacity for each record; unique
parent paths between neighbours; and, keeping the latest record of each level,
that every record's parent precedes it with a child slot for it. The number of
such reachable records must equal the number of reachable branches, so every
reachable branch has exactly one record and stale or unreachable records are
rejected. This runs on every constrained reconstruction, so a state whose
active-count metadata is malformed fails its next inversion. An
empty/unconstrained request has a distinct meaning from a partially populated
topology.

`recursive_branch_node_count_limit` finds a branch's record by binary search on
its parent path, which names the branch because its length is the conditioning
level. It returns the previous active count for an unchanged, constrained
branch, or configured capacity if no limit applies.
`branch_limit_invalidated` checks whether a changed ancestor path prefixes the
branch. Ancestor changes invalidate descendant constraints, since those
children belong to a newly reconstructed conditional distribution. They are not
necessarily limited to a one-node increase from their previous records.

`branch_growth_trial_allowed` permits one additional node only after an
adaptive branch reconstructs successfully **without reducing its current
count**, when growth is requested and configured capacity remains. A successful
exact-count trial replaces the current candidate; failed or reduced growth
trials leave the current accepted candidate intact. A reduction does not
immediately trigger growth. Bounded branches use `select_bounded_quadrature`
with the same distinction between current count and growth ceiling. Bounds
apply to reconstructed nodes; they do not clip nodes or establish support of
the underlying distribution.

Evidence: [recursive unit tests](../../../tests/unit/recursive_quadrature.cpp)
check that supplied counts in any order are accepted and that stale, missing,
duplicate, and over-capacity records are rejected; [phase-space unit
tests](../../../tests/unit/phase_space_evolution.cpp), including
`test_adaptive_state_parent_growth_creates_unconstrained_descendants`,
`test_adaptive_state_failed_growth_trial_keeps_current_count`, and
`test_adaptive_state_reduction_prevents_source_growth_trial`, explicitly
exercise these policies. Root and second-level collapse tests check descendant
invalidation. [Adaptive
verification](../../../tests/verification/cqmom/adaptive_branch_growth.cpp)
provides additional workflow coverage. Repeated threshold crossings can still
produce topology changes across successive reconstruction steps.

## Recursive driver — Level A

`compute_recursive_quadrature_from_validated_moment_data` is the file-local
numerical backbone called by the validated and trusted moment-data entry points.
Inputs are raw global moments, hierarchy, Wheeler and conditional-solve options,
support intervals, diagnostic policy, and optional active topology/growth
permission. Outputs include quadrature storage, the constructed tree, status,
failure location, selected branch diagnostics and, when the request asks for
them in adaptive mode, the active branch counts. Configured counts, the current
trial's requested count, and returned active counts are separate values.

1. The entry points have already rejected invalid options. The driver
   validates coordinate support and any supplied active branch counts; invalid
   support, topology, or metadata becomes an invalid request.
2. Construct or reuse the ordered root layout and gather the root slice. The root
   starts with global weight one and a coordinate vector of the hierarchy rank.
3. The outer loop visits conditioning levels; the inner loop visits current
   branches. Gather the pure leading moments and invert them with Wheeler, with
   bounds/active ceiling/growth as applicable. A failed branch ends reconstruction
   immediately and records its location even when successful diagnostics are off.
4. For a nonterminal level, recover conditional slices for all accepted parent
   nodes. For every local node, multiply the accumulated parent weight by its
   conditional weight, place its abscissa in the corresponding physical dimension,
   append its conditioning index to the path, and schedule its child slice.
5. The final level emits terminal weights and coordinates. Adaptive storage is a
   ragged terminal list. Fixed storage is indexed using the phase-space rectangular
   multi-index. Termination follows the finite hierarchy rank and each accepted
   finite node list; there is no convergence loop over CQMOM levels.
6. When the request records active branch counts, each accepted adaptive branch
   adds its level, parent path, dimension, configured count, and accepted node
   count. The traversal visits branches level by level; on success the counts
   are sorted by parent path, which lists each branch before its descendants,
   the order in which the terminal paths visit them. A `CqmomState` takes its
   counts from here.

Failure leaves partial tree/storage as diagnostic data. Public result
conversion maps status and failed dimension; consumers must check success
before treating that data as a quadrature. Exceptions caught by the driver map
to failure status; this is not a promise that every lower-level API is
nonthrowing.

Numerical sensitivity compounds across levels: conditional weights enter the
linear solve, conditional moments enter another nonlinear inversion, and
terminal weights multiply along paths. Conditioning order changes these
intermediate problems even for an unchanged global moment set. Small weights,
clustered nodes, near-degenerate marginals, and extreme coordinate scales
therefore deserve joint experiments.

Evidence: [recursive tests](../../../tests/unit/recursive_quadrature.cpp)
inspect failure diagnostics, partial branches, support tolerances, and policy
parity; [analytical high-dimensional
verification](../../../tests/verification/cqmom/cqmom_high_dimensional_analytical.cpp)
and [recursive
verification](../../../tests/verification/cqmom/recursive_cqmom.cpp) recover
independently specified distributions/moments. Fixed terminal paths enumerate
physical row-major indices before conversion to conditioning paths, so they
align with fixed weights and abscissas under non-identity order.

## Entry points and helpers — Levels B/C

| Functions | Responsibility and contract | Evidence |
|---|---|---|
| `compute_recursive_quadrature_from_moment_data` | Validate canonical required moment data against the hierarchy, then route to the recursive driver. | Recursive/public API unit tests |
| `compute_recursive_quadrature_from_trusted_moment_data` | Reuse established structural metadata; retain required numerical-value checks. Trusted does not mean arbitrary user data is valid. | Fixed metadata and evolution tests |
| `create_ordered_root_slice`, `one_dimensional_moments`, `child_layout_for_parent_node_count` | Take the root values in their canonical stored order, gather leading values, and select precomputed child metadata; no new solver policy. | Recursive layout tests |
| `record_failure`, `create_recursive_result`, `record_bounded_candidate_diagnostics`, `records_success_branch_diagnostics` | Keep failure location, candidate rejection, and requested diagnostic detail distinct from success acceptance. | Failure/diagnostic tests |
| `valid_conditional_options`, `valid_diagnostic_level`, `validate_recursive_request_preflight`, `validate_coordinate_support`, `support_interval_for_dimension` | Resolve valid enum/options and per-coordinate bounds before numerical work. No physical meaning inferred from coordinate names. | Invalid-request and support tests |
| `path_has_prefix`, `parent_path_less`, `append_unique_index`, `local_branch_offset` | Deterministic record ordering, ancestor membership, first-occurrence layout deduplication, and checked offset discovery. | Topology/layout tests indirectly |

Index arithmetic, hierarchy methods, and tree mutation are covered in
[state and representations](state_and_representations.md); remaining exact
helper/overload locations are routed by the inventory.
