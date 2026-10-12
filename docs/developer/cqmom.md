# CQMOM

CQMOM [[Yuan & Fox, 2011](https://www.doi.org/10.1016/j.jcp.2011.07.020)]
reconstructs a multidimensional quadrature one coordinate at a time, each
conditioned on the coordinates before it. QMOMKit applies this recursively
through all coordinates
[[Bouguezzoul, 2026](https://publications.polymtl.ca/76505/)]. The user-level
picture is in [Concepts](../user/3_concepts.md#qmom-and-cqmom).

The hierarchy stores node counts in phase-space order, and the conditioning
order is a separate permutation that only sets the traversal; see
[Numerical conventions](numerical_conventions.md). Below, $c_l$ is the
coordinate at conditioning level $l$ and $N_l$ its configured node count.

## The recursion

A *branch* at level $l$ is identified by its parent path, the node indices
chosen at levels $0,\ldots,l-1$. Each branch carries a *slice*: the moments it
needs, indexed in conditioning order with $c_l$ first.

1. **Root.** The root slice is the required moment data, reordered into
   conditioning order. Its leading sequence $M_{k\,e_{c_0}}$, $k<2N_0$, is the
   marginal of $c_0$, which the Wheeler algorithm
   [[Wheeler, 1974](https://www.doi.org/10.1216/RMJ-1974-4-2-287)] inverts as
   in [QMOM](qmom.md). The root weights add up to $M_{0,\ldots,0}$.
2. **Conditional moments.** For a nonterminal branch whose inversion accepted
   $n$ nodes $\xi_i$ with weights $w_i$, the slices of its $n$ children solve one
   weighted Vandermonde system with one column per trailing multi-index
   $q$:

   ```math
   \sum_{i=0}^{n-1} w_i\, \xi_i^{p}\, C_{i,q} = M_{p,q},
   \qquad p = 0,\ldots,n-1,
   \qquad\text{that is}\qquad
   V^{\mathsf T}\operatorname{diag}(w)\, C = M,\quad V_{ip} = \xi_i^{p}.
   ```

   $M_{p,q}$ is the slice moment of order $p$ in $c_l$ and $q$ in the later
   coordinates; row $i$ of $C$ is the slice of child $i$. All columns
   are solved together by `solve_conditional_vandermonde_system`;
   [Linear algebra](linear_algebra.md) describes how. Only parent orders
   $p<n$ enter, so a branch with fewer nodes than configured solves a square
   system of its own size. Column $q=0$ reproduces moments that Wheeler
   matched, so $C_{i,0}=1$ up to rounding: child slices are normalized
   conditional moments, and weights below the root are conditional
   probabilities.
3. **Next level.** Wheeler inverts the leading sequence of every child slice.
   Conditional moments carry the rounding of the solve that produced them,
   which is relative to the scale of their coordinate, not to their own
   magnitude. In adaptive mode each inversion therefore receives the scale $s$
   of the coordinate it inverts, which sets the floor of its rounding bound
   ([QMOM](qmom.md#rounding-bound-on-the-recurrence)); how $s$ is computed is
   in [QMOM](qmom.md#coordinate-scales). The root inverts given moments and
   uses $s=0$, as does every level in fixed mode.
4. **Terminal nodes.** At the last level, every accepted node becomes a
   terminal node. Its weight is the product of the local weights along its
   path, and its abscissas are collected in phase-space order.

The traversal visits levels in order and branches within a level in path
order. It stops at the first failing branch. It also records a traversal tree
in `internal::RecursiveQuadrature`, which only tests read. The driver is
`compute_recursive_quadrature_from_validated_moment_data` in
[`src/recursive_quadrature.cpp`](../../src/recursive_quadrature.cpp); the
conditional step is the file-local `reconstruct_child_slices`.

## Required moments

`create_cqmom_moment_layout` lists the moments the recursion reads. At level
$l$ it enumerates, in conditioning order:

| Coordinates | Orders | Reason |
| --- | --- | --- |
| $c_j$, $j<l$ | $0,\ldots,N_j-1$ | rows of the conditional systems at earlier levels |
| $c_l$ | $0,\ldots,2N_l-1$ | the Wheeler sequence at level $l$ |
| $c_j$, $j>l$ | $0$ | not yet reached |

The union over levels, deduplicated in first-occurrence order, is stored in
phase-space order in `CqmomMomentLayout::required_indices`. Each level adds the
indices whose order in $c_l$ is positive, so the layout holds
$2N_0 + \sum_{l\ge1}\bigl(\prod_{j<l}N_j\bigr)(2N_l-1)$ moments; both the set
and its size depend on the conditioning order. The layout is the same in
adaptive mode, where a reduced branch reads only part of it.

`CqmomMomentData` aligns one value with each required index. It is built with
`extract_cqmom_moments` from a full tensor, which throws `std::out_of_range`
when a required moment lies outside the tensor, or with
`create_cqmom_moment_data` from a successful `SelectedRawMomentResult` with
exactly the layout's indices in layout order. A missing moment is never
treated as zero. Reconstruction requires exactly the required layout of the
hierarchy and compares against it once: when a `CqmomPreparedCase` is
constructed, in `CqmomCase::initial_state()`, `compute_cqmom(cqmom_case)`,
`state_candidate(...)` and `config::with_initial_moments(...)`. An aggregate
`CqmomCase` checks nothing when it is built. The root slice then takes the
values in stored order without lookups.

### Branch layouts

`create_branch_moment_layout` turns a slice's indices into lookup metadata: the
offsets of the Wheeler sequence, a table from each multi-index to its position, and,
for each parent node count, a *child case*: the slice offsets that form each
right-hand-side column and the child slice's own layout, built recursively.
Child cases depend only on the parent node count, so every branch with the
same accepted count uses the same one. A fixed layout has a child case only
for the configured count; an adaptive layout has one for every count from one
to the configured maximum.

## Fixed and adaptive reconstruction

`QmomOptions::adaptive` decides both whether each Wheeler inversion may
reduce its node count and whether branches may differ in size. In fixed mode, a
branch that returns any other count than its configured one fails with
`qmom_failed`. In adaptive mode, each branch accepts at most its ceiling, and
the children of a reduced branch use the child case of the accepted count. The
two storage layouts are in
[Numerical conventions](numerical_conventions.md#quadrature-storage).

### Active branch counts

A reconstruction that creates an adaptive `CqmomState` records one
`CqmomActiveBranchCount` per accepted branch at every level: the level, the
parent path, the phase-space dimension, the configured count and the accepted
count.

A later reconstruction can receive these records. `validate_active_branch_counts`
then checks, in one pass over the sorted records, that each matches the
hierarchy, that every reachable branch has exactly one record and that no
record belongs to an unreachable branch; records in another order are sorted
into a copy first, and a malformed set is `invalid_input`. A fixed
reconstruction ignores the records. `recursive_branch_node_count_limit` makes
each branch's recorded count its ceiling. If growth is allowed, which only the
post-source reconstruction does, a branch that returns exactly its recorded
count, below its configured maximum, also tries one more node. It
keeps that rule only if Wheeler returns all of its nodes.

A branch whose accepted count differs from its record has different children:
they condition on other nodes, so the records below it do not describe them.
The traversal therefore drops the records of its whole subtree, and those
branches are reconstructed up to their configured counts. Which
reconstructions of phase-space evolution pass the records, and whether they
allow growth, is in
[Phase-space evolution](phase_space_evolution.md#node-counts).

### Bounds

`CqmomAbscissaBounds` holds one optional `AbscissaBounds` per phase-space
dimension, or nothing at all. A non-empty vector of the wrong length, or an
interval that is not finite with lower below upper, is `invalid_input`. Every
branch of a bounded coordinate goes through the same `select_bounded_quadrature`
as QMOM, with that branch's ceiling, recorded count and growth permission; see
[QMOM](qmom.md#abscissa-bounds).

## Prepared cases

`create_prepared_cqmom_metadata` builds a `PreparedCqmomMetadata`: the
required moment layout, the root branch layout with every child case below
it, and the adaptive flag. It holds no moment values, quadratures, or solver
workspaces, and depends only on the hierarchy and the flag. The adaptive
structure holds a child case for every count Wheeler can return, so it too is
built once and never depends on the moments.

`CqmomPreparedCase` builds it once in its constructor, which throws
`std::invalid_argument` for a layout other than the required one, and shares it through a
`std::shared_ptr<const PreparedCqmomMetadata>` with:

- every `compute_cqmom(prepared_case)`;
- every state from `initial_state()` or `state_candidate(...)`, their copies
  and every evolution successor;
- every case from `config::with_initial_moments(...)`, which keeps the
  hierarchy and the flag and so stays compatible.

`prepared_cqmom_metadata_is_compatible` compares node counts, conditioning
order, required indices and the flag; a prepared case or state given
incompatible metadata throws `std::invalid_argument`. Every reconstruction of
a prepared case or state, including those inside evolution steps, goes
through `compute_cqmom_internal`, which takes the metadata of the lineage by
reference and never builds branch layouts; metadata prepared for the other
adaptivity is `invalid_input`.

`CqmomCase` holds no metadata. `compute_cqmom(cqmom_case)` checks the layout
and builds the branch layouts on every call; `CqmomCase::initial_state()`
builds metadata once for the state it returns, which then shares it with its
successors. Both routes give the same quadrature.

## Diagnostics

`CqmomResult::diagnostics` holds `CqmomDiagnostics` records chosen by
`CqmomDiagnosticLevel`:

| Level | Records |
| --- | --- |
| `none` (default) | Only the failing branch, if a branch failed. |
| `basic` | One record per reconstructed branch: level, parent path, dimension, `requested_node_count` and accepted count, `adaptive_node_count_reduced`, Wheeler status as `QmomStatus`, conditional-moment status, requested and selected Vandermonde solver, and bounds details. |
| `detailed` | As `basic`, plus the relative residual and condition estimate of each conditional system and `conditioning_warning` when the estimate exceeds `ConditionalMomentOptions::conditioning_warning_threshold`. |

The detailed quantities describe the scaled system the solvers decide on
([Linear algebra](linear_algebra.md#diagnostics)), so they do not depend on
the number density, and a change of unit moves the condition estimate by a
bounded factor only. They need the explicit matrix and its full inverse, which
is why they have their own level, and they never fail a solve. The level never
changes a quadrature, a status or a node count.

`requested_node_count` is the branch's ceiling, raised by one when a growth
trial is accepted; for a bounded coordinate it is the candidate count the
bounds accepted. The selected solver is the one `automatic` resolves to.
Records of the last level have no conditional solve, so their solver fields
keep their defaults, `automatic` and `dense`, which mean nothing there.

## Statuses

| `CqmomStatus` | Cause |
| --- | --- |
| `invalid_input` | Invalid conditional options or diagnostic level, a layout or values that do not fit the hierarchy, non-finite moments, invalid bounds or active counts, node counts whose $2N$ overflows. |
| `qmom_failed` | A branch's Wheeler inversion failed, or a fixed branch returned another count. |
| `conditional_moment_failed` | A branch's conditional system could not be solved or had non-finite inputs. |

`failed_phase_space_dimension` is the coordinate of the failing level and has
no meaning on success or for input rejected before the traversal. Check `ok()`
first: after a failure, rectangular storage holds zeros where no terminal node
was assembled and ragged storage is empty; either is diagnostic only.

## Moments from a quadrature

[`src/moment_reconstruction.cpp`](../../src/moment_reconstruction.cpp) evaluates
$\sum_j W_j\prod_d \xi_{d,j}^{k_d}$ for both QMOM and CQMOM:

- `reconstruct_raw_moments` returns the first `moment_count` moments of a
  successful `Quadrature1D`, a full CQMOM tensor of given extents, or the
  moments of a
  `RawMomentRequest` in request order with duplicates kept;
- `reconstruct_cqmom_moments` returns the required moments of a layout, from
  rectangular or ragged storage;
- `evaluate_generalized_moment` takes finite real orders; a non-integer order
  needs non-negative abscissas.

These functions throw on malformed input or overflow
([Architecture](architecture.md#errors-and-statuses)). States store the moments of their
accepted quadrature: QMOM states through `reconstruct_raw_moments`, CQMOM
states through `internal::reconstruct_trusted_cqmom_moments`, which skips the
layout checks. Source contexts and output writers use the public functions.

## Where the code lives

| Function | File | Role |
| --- | --- | --- |
| `compute_cqmom`, `compute_cqmom_internal` | [`src/public_api.cpp`](../../src/public_api.cpp) | Public entries; status and diagnostics mapping. |
| `create_cqmom_moment_layout`, `extract_cqmom_moments`, `create_cqmom_moment_data` | [`src/recursive_quadrature.cpp`](../../src/recursive_quadrature.cpp) | Required moments. |
| `create_branch_moment_layout`, `create_prepared_cqmom_metadata` | [`src/recursive_quadrature.cpp`](../../src/recursive_quadrature.cpp) | Branch layouts and prepared metadata. |
| `compute_recursive_quadrature_from_moment_data`, `..._from_trusted_moment_data` | [`src/recursive_quadrature.cpp`](../../src/recursive_quadrature.cpp) | Checked and trusted entries into the recursion. |
| `validate_active_branch_counts`, `recursive_branch_node_count_limit` | [`src/recursive_quadrature.cpp`](../../src/recursive_quadrature.cpp) | Active counts and ceilings. |
| `solve_conditional_vandermonde_system` | [`src/internal/vandermonde.cpp`](../../src/internal/vandermonde.cpp) | Conditional systems. |
| `reconstruct_raw_moments`, `reconstruct_cqmom_moments`, `evaluate_generalized_moment` | [`src/moment_reconstruction.cpp`](../../src/moment_reconstruction.cpp) | Moments from a quadrature. |

Coverage is in `tests/unit/recursive_quadrature.cpp`,
`tests/unit/public_api.cpp`, `tests/unit/moment_reconstruction.cpp`,
`tests/component_regression/` and `tests/verification/cqmom/`.
