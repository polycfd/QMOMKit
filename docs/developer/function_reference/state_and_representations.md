# State, hierarchy, and representations

Public ownership lives in [qmom.hpp](../../../include/qmomkit/qmom.hpp),
[cqmom.hpp](../../../include/qmomkit/cqmom.hpp),
[conditioning_path.hpp](../../../include/qmomkit/conditioning_path.hpp), and
related moment/quadrature headers listed in the
[inventory](coverage_inventory.md). Definitions are primarily
[public_api.cpp](../../../src/public_api.cpp),
[recursive_quadrature.cpp](../../../src/recursive_quadrature.cpp), and
[conditioning_path.cpp](../../../src/conditioning_path.cpp). Private
construction is declared in
[cqmom_case.hpp](../../../src/internal/cqmom_case.hpp) and
[qmom_state_access.hpp](../../../src/internal/qmom_state_access.hpp). See
[numerical conventions](../numerical_conventions.md).

## Indices and capacity — Level A

`product_extents` validates positive extents and checks multiplication
overflow. `flatten_multi_index` validates rank and coordinate bounds, then
repeatedly multiplies by the next extent and adds the next index. The **last
physical dimension varies fastest**. `unflatten_multi_index` checks the flat
range and reverses that operation by remainder/division.
`enumerate_multi_indices` iterates flat indices in this order; it does not use
conditioning order.

`PhaseSpaceHierarchy` constructors validate positive configured counts and a
permutation of physical dimensions (identity when omitted). `moment_extents`
returns twice each configured count; capacities are not active branch counts.
`flatten_node_index` and `unflatten_node_index` use phase-space order.
`ConditioningPath` instead stores local node indices in conditioning order.
Path extension/prefix access describes tree ancestry, not a dense offset.
Required extents are checked before dependent storage is constructed.

`CqmomMomentData` associates canonical required indices with values. Its lookup
returns the stored raw moment at an exact required index; missing entries are
not interpolated. Dense `RawMomentTensor` access uses extents/flat offsets;
full tensor storage contains moments not needed for inversion. Extraction into
required storage preserves raw zeroth density and the canonical enumeration
order.

Evidence: [recursive tests](../../../tests/unit/recursive_quadrature.cpp),
[conditioning-path cases](../../../tests/unit/recursive_quadrature.cpp), and
[public API tests](../../../tests/unit/public_api.cpp) cover index round trips,
non-identity order, invalid paths/layouts, and missing required moments. These
small-grid tests do not establish allocation limits for large capacities.

## Quadrature storage and tree mutation — Level A

Internal `RecursiveQuadrature` fixed construction allocates rectangular weights
and per-dimension abscissas. `set_weight`/`set_abscissa` write by physical
multi-index. `add_tree_node` assigns an id and attaches the child id to its
parent; node records separately carry level, conditioning path, accumulated
weight, and coordinate information. Parent ids are an internal construction
invariant.

`set_adaptive_branch_counts_enabled` switches to ragged operation and clears
rectangular terminal arrays. `set_terminal_quadrature` checks terminal shapes
and path count before installing the supplied terminal vectors. Adaptive output
uses those aligned traversal-order vectors. Fixed reconstruction fills
rectangular values and creates paths from the same physical row-major indices,
then converts them to conditioning paths.

Public `CqmomQuadrature` owns hierarchy, terminal values, paths, and flags. Its
internal construction factory transfers these fields; it is not an independent
numerical validator. `weight`, `abscissa`, and `node` by multi-index require
rectangular storage and throw for ragged storage. Iterate terminal arrays for
an adaptive result. Neither the container nor a nonempty array alone proves
that the producing reconstruction succeeded.

## Deriving active counts from paths — Level A

Internal `active_branch_counts_from_quadrature`, defined in
[cqmom_active_counts.hpp](../../../src/internal/cqmom_active_counts.hpp),
returns an empty constraint vector for fixed quadratures. For adaptive
quadratures it checks dimension and terminal-path count, then loops over every
terminal path and conditioning level. Each encountered child index is
accumulated once under its branch key `(level,parent path,physical dimension)`.
`parent_path_at_level` constructs that prefix; `find_active_branch_accumulator`
and `contains_node_index` perform linear searches preserving first-encounter
order.

`validate_contiguous_child_indices` sorts each branch's unique child indices
and requires `0,1,...,n-1`; the resulting count is therefore more than an
unchecked maximum-index-plus-one. Out-of-capacity, incomplete, or noncontiguous
paths throw. The resulting records describe reachable topology, not inferred
numerical rank. They do not check path/value alignment; fixed mode bypasses
this derivation. Evidence: topology and malformed-path cases in [phase-space
tests](../../../tests/unit/phase_space_evolution.cpp) and [recursive
tests](../../../tests/unit/recursive_quadrature.cpp). Large ragged hierarchies
could make the linear accumulator searches material; no scaling measurement was
performed.

## Case and accepted state — Levels A/B

`QmomCase::initial_state` calls `compute_qmom` and throws if reconstruction
fails. `QmomStateAccess::from_quadrature` rebuilds the stored raw sequence, at
its original length, from the accepted quadrature. Consequently adaptive
acceptance may replace higher input moments with the reduced quadrature's
closure. `QmomState` exposes moments, quadrature, bounds, active count,
configured maximum node count, and read-only reconstruction options. The
maximum and reconstruction options are copied from `QmomCase` and remain
authoritative for every subsequent evolution transition.

`CqmomCase` and `CqmomState` hold a `CqmomLeanCase` and do not expose it.
Public accessors read the hierarchy and policies from it; numerical algorithms
use the held case directly. A `CqmomLeanCase` is promoted through the public
`CqmomCase` constructor.
`CqmomCase::initial_state` computes a successful quadrature, derives active
branch counts, reconstructs the stored representation from that accepted
quadrature, and stores the same quadrature in the state.
`CqmomState::quadrature()` exposes that stored accepted quadrature by const
reference. A state reconstructs exactly its required layout in its existing
order. Thus every successfully constructed initial state has stored moments and
cached accepted quadrature describing the same current state. Adaptive
reduction may replace higher supplied moments with the reduced quadrature's
closure, as it already does for QMOM and later CQMOM accepted stages. Raw and
generalized moment views are therefore views of the same accepted state.

`CqmomConstructionAccess::{create_case,create_state}` enforces fixed metadata
ownership/compatibility and creates the variant payload. Fixed cases eagerly
prepare immutable hierarchy/layout metadata; copies deep-copy case/state data
while sharing that immutable metadata. Successors retain the same fixed
metadata when the layout is unchanged. Adaptive variants reject inappropriate
fixed metadata. `CqmomConstructionAccess` supplies internal variant/cache
references; it is a trusted access boundary, not an ordinary public
representation API.

`CqmomState::moment_indices`, `moment_values`, and `moment` expose the
required moment data. A moment outside that set returns no value;
source-context closure is a separate operation. A mutable `CqmomLeanCase` must
pass public validation before normal reconstruction. Internally trusted calls rely on their creator's structural
checks.

Evidence: [configuration tests](../../../tests/unit/configuration.cpp) inspect
case construction/replacement, copying, metadata reuse, and validation;
[phase-space tests](../../../tests/unit/phase_space_evolution.cpp) cover fixed
and adaptive cache ownership, stored-moment versus
accepted-quadrature consistency, source-context agreement, and zero-source
invariance.

## Public dispatch and conversion — Level B

`compute_qmom` overloads resolve case/options/bounds and map Wheeler status;
[QMOM reconstruction](qmom_reconstruction.md) explains scaling and selection.
`compute_cqmom` validates a public `CqmomLeanCase` or dispatches an opaque case
to `compute_cqmom_internal`. Fixed normal workflows supply prepared metadata;
adaptive workflows construct count-dependent layouts. Scaling wraps the
internal reconstruction with reversible coordinate transformations.

`to_qmom_status`, `to_cqmom_status`, `to_cqmom_diagnostics`, `to_cqmom_result`,
and `failed_phase_space_dimension` translate internal status/location and
requested diagnostic fields. They do not retry a failed solver or independently
validate a numerical result. `terminal_paths_for_public_quadrature` preserves
existing paths when present and synthesizes them otherwise. The separate
diagnostics and value conversion steps explain why diagnostic data can
accompany a failed result.

## Access and ownership helpers — Level C

| Family | Purpose and significant semantics | Evidence |
|---|---|---|
| Case/state constructors, assignments, destructors | Own variant data; ordinary copying preserves independent evolving values. Moves follow C++ ownership rules; access to a moved-from opaque payload is not established here as a supported workflow. | Public API/configuration tests |
| Hierarchy, option, support, storage, weight, abscissa, path, and count getters | Return the corresponding stored metadata/value/reference; no recomputation or realizability checks. Reference lifetime follows the owner. | Public-header and API tests |
| `CqmomConstructionAccess` case and metadata accessors | Internal access and cache sharing; they hand out the held case by reference rather than copying it. | State/metadata tests |
| `default_conditioning_order`, `ordered_extents`, `from_ordered_index`, `phase_space_node_index_to_path`, `phase_space_index_to_conditioning_path`, `rectangular_terminal_paths` | Translate the two index conventions; fixed paths are generated from physical row-major indices to align with rectangular value offsets. | Non-identity hierarchy tests |
| `validate_positive_extents`, `checked_multiply`, `validate_conditioning_order`, `validate_abscissas` | Structural preconditions for indexing/allocation and finite quadrature storage. | Invalid-input tests |
| Result `ok`/boolean helpers; `format_cqmom_settings` | Read status or serialize selected settings. They do not establish numerical acceptance beyond the stored status. | API/configuration formatting tests |

Every concrete definition/overload is locatable through the inventory, including
inline moment-container and result helpers that need no separate algorithm section.
