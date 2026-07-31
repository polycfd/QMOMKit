# Recursive CQMOM

This page documents the current QMOMKit design for multidimensional CQMOM
construction. It is a developer reference for maintainers, contributors, and
AI-assisted development workflows.

For raw-moment definitions, tensor storage, phase-space dimension order,
conditioning order, and node/weight terminology, see
`numerical_conventions.md`. For the one-dimensional Wheeler inversion used at
each recursive branch, see `qmom_wheeler.md`.

## Role Of Recursive CQMOM

QMOM reconstructs a one-dimensional quadrature rule from a one-dimensional raw
moment sequence. CQMOM extends this idea to multidimensional raw moment
tensors by reconstructing dimensions conditionally. A multidimensional moment
has the form:

```text
m_alpha = sum_j W_j product_r xi_(r,j)^alpha_r
```

where `alpha` is ordered by phase-space dimension. Directly inverting a full
multidimensional tensor is not the QMOMKit design. Instead, CQMOM processes
one phase-space dimension at a time. The first processed dimension is inverted
with Wheeler. Later dimensions are reconstructed as conditional quadrature
rules under the already selected parent nodes.

Recursive CQMOM is the implementation strategy that applies this conditional
construction through all dimensions in the configured conditioning order. The
public API remains CQMOM:

```text
compute_cqmom(...)
```

The recursion is how QMOMKit computes CQMOM internally; users still call the
CQMOM API.

## Conditional Moment Representation

CQMOM separates two orderings:

- phase-space dimension order, used by `RawMomentTensor`, terminal
  coordinates, node counts, and public coordinate vectors;
- conditioning order, used by the recursive construction sequence.

At a conditioning level:

```text
phase_space_dimension = conditioning_order[level]
```

Node counts are owned by phase-space dimensions:

```text
node_counts[r] == requested node count for phase-space dimension r
```

For a parent dimension with known nodes and local weights:

```text
x_i, w_i, i = 0, ..., n_parent - 1
```

CQMOM reconstructs conditional raw moments for the next dimension by solving a
weighted Vandermonde system:

```text
A * C = M
```

where:

```text
A[p, i] = w_i * x_i^p
C[i, q] = c_i,q
M[p, q] = mixed raw moment data
```

Here `c_i,q` is the conditional raw target moment of order `q` at parent node
`i`. Rows of `A` and `M` are parent moment orders, columns of `A` are parent
nodes, and columns of `C` and `M` are target moment orders or flattened
remaining branch-local moment indices.

The public one-level function `reconstruct_conditional_moments(...)` exposes
the weighted Vandermonde convention for global raw moment tensors. Recursive
CQMOM uses the same mathematical solve on branch-local moment slices while
walking the full conditioning tree.

## Root Moment Storage

QMOMKit's scalable CQMOM representation is the required moment set: the minimal
set of raw moments required by the configured hierarchy and conditioning order.
The normal user workflow loads a `CqmomCase` and computes its quadrature with
`compute_cqmom(CqmomCase)`. If the distribution evolves in time, the same case
can create a `CqmomState` for transport. Users on this path do not manage
internal branch-local storage, dense tensor extents, or required moment indices.
Representation choices are hidden by default. Callers that need direct control
can opt into an explicit dense representation workflow or an explicit lean
representation workflow.

### Explicit Dense Representation Workflow

The public dense CQMOM case accepts:

```text
RawMomentTensor
```

as an explicit full rectangular tensor moment representation. A dense tensor is
convenient for callers that need full tensor access, already own rectangular
moment data, exchange data with other tools, debug tensor values, or validate
against full manufactured tensors. Dense is a first-class supported explicit
representation workflow. The dense public representation remains
phase-space ordered:

```text
values[flatten_multi_index(alpha, extents)] == m_alpha
```

Recursive CQMOM does not require every entry of that rectangle. The public
`compute_cqmom(CqmomDenseCase)` path validates the dense tensor, constructs the
fixed-node CQMOM requirement layout, selects the required entries, and then
delegates using the lean CQMOM moment representation. This keeps dense tensor
control at the API boundary while avoiding dense tensor lookup as the recursive
working representation.

### Explicit Lean Representation Workflow

The lean root representation is described by:

```text
CqmomMomentLayout
CqmomMomentData
```

`CqmomMomentLayout` stores the phase-space hierarchy and the deterministic list
of phase-space-order raw moment multi-indices required by fixed-node recursive
CQMOM. It stores no values and is not a general sparse tensor description.

`CqmomMomentData` stores one value per layout index:

```text
values[i] == m_(layout.required_indices[i])
```

The values are aligned with the layout order. Public helper functions construct
this data either by extracting from a dense `RawMomentTensor` or by validating a
selected raw-moment result produced elsewhere. The public case type named:

```text
CqmomLeanCase
```

accepts `CqmomMomentData` directly for callers that can produce the lean moment
set directly.

Inside recursive reconstruction, lean moment data is adapted at the root
boundary to:

```text
IndexedMomentSet
CqmomMomentView
```

`IndexedMomentSet` is internal lookup storage. It preserves deterministic
construction order for indices and values, rejects duplicate keys, and builds a
sorted flat-key offset table for moment lookup. It is intentionally scoped to
CQMOM-required moment storage, not arbitrary sparse tensor manipulation.

`CqmomMomentView` is an internal non-owning adapter over either a dense
`RawMomentTensor` or an `IndexedMomentSet`. Recursive code reads moments through
the view so the traversal is independent of the owning root representation.
The lean public path stores indexed moments for the duration of reconstruction,
then delegates through this same view-based recursive algorithm.

## Validated And Reused Setup Paths

Recursive CQMOM has two validation boundaries.

Explicit representation workflows perform full metadata and value validation
before reconstruction:

- user-provided `CqmomMomentData`;
- explicit `CqmomLeanCase`;
- explicit `CqmomDenseCase`.

These paths validate hierarchy consistency, conditioning order, the required
fixed CQMOM moment layout, value count, and finite moment values. The dense path
also validates rectangular tensor extents before extracting the required lean
moment set.

`CqmomCase` and `CqmomState` are normal public workflows. For these objects,
QMOMKit can reuse setup information that was already checked when the case or
state was created or loaded. Reconstruction therefore reuses the existing
hierarchy, conditioning order, and moment layout instead of regenerating and
comparing the global required layout on every reconstruction.

Reusing setup information does not mean unchecked numerical data. These paths
still validate evolving moment values, vector sizes, finite values, runtime
numerical options, and reconstruction failure conditions. Branch layouts,
adaptive branch counts, Wheeler inversions, conditional Vandermonde solves,
terminal assembly, statuses, and diagnostics follow the same recursive
algorithm as the explicit representation paths.

## Recursive Construction

The recursive construction starts from a root CQMOM moment view whose indices
are still phase-space-order raw moment indices. Internally, the root moment
data is reordered into a branch-local slice whose axes follow conditioning
order. This reordering is internal; public tensor storage and lean layout
indices remain phase-space ordered.

At each conditioning level, the implementation:

1. Determines the phase-space dimension for the level.
2. Extracts the leading one-dimensional raw moment sequence from the current
   branch-local slice.
3. Calls Wheeler inversion with the node count for that phase-space dimension.
4. Uses the returned Wheeler nodes and local weights to extend the recursive
   branch path.
5. If more dimensions remain, solves the branch-local weighted Vandermonde
   problem to build one child moment slice per Wheeler node.
6. If this is the final dimension, stores the accumulated terminal quadrature
   node.

The leading-axis moment sequence passed to Wheeler is:

```text
slice[order, 0, 0, ...]
```

where the leading axis corresponds to the current conditioning level, not
necessarily phase-space dimension `0`.

The branch-local Vandermonde solve treats the current Wheeler nodes as parent
nodes. Its right-hand side columns are the flattened trailing axes of the
branch-local slice. Each solved parent-node row becomes the child
branch-local moment slice for one recursive child.

Wheeler weight interpretation depends on the conditioning level. At the first
conditioning level, Wheeler produces the marginal quadrature weights for the
first processed phase-space dimension. At subsequent levels, Wheeler produces
branch-local conditional weights associated with the current parent path. The
implementation forms accumulated terminal weights by multiplying the current
level's weight into the incoming branch weight:

```text
child_accumulated_weight =
    parent_accumulated_weight * local_weight
```

The root accumulated weight is `1`. At a terminal branch, the accumulated
weight is the global multidimensional quadrature weight.

## Vandermonde Solver

Conditional moment reconstruction requires solving weighted Vandermonde
systems. At each non-terminal recursive branch, Wheeler provides parent
abscissas and local weights:

```text
x_i, w_i, i = 0, ..., n_parent - 1
```

The branch-local reconstruction solves:

```text
A * C = M
```

with:

```text
A[p, i] = w_i * x_i^p
C[i, q] = c_i,q
M[p, q] = mixed branch-local raw moment data
```

The row index `p` is the parent moment order. The column index `i` is the
parent node. The right-hand-side column `q` is a flattened trailing
branch-local moment index.

QMOMKit calls the numerical component that solves this system the Vandermonde
solver. The user-facing solver choices are:

- Dense solver: a generic dense linear solve used as the robustness reference.
- Björck-Pereyra Vandermonde solver: the Vandermonde-specific
  Björck-Pereyra algorithm without dense fallback.
- Automatic solver: the default. It selects the fallback-capable
  Björck-Pereyra Vandermonde solver.
- Björck-Pereyra with dense fallback: tries the Vandermonde-specific
  Björck-Pereyra algorithm, validates the residual when requested by
  diagnostics or automatic mode, and uses the dense solver if the fast solve
  fails validation or returns a failure.

The distinction between requested solver and actual algorithm is deliberate.
Diagnostics record the requested Vandermonde solver, the selected solver after
automatic resolution, the algorithm actually used (`dense` or
`bjorck_pereyra`), and whether dense fallback occurred. For example, the
default request is `automatic`; recursive CQMOM resolves that to
`bjorck_pereyra_with_dense_fallback`; a difficult branch may then report that
the dense algorithm was actually used and that dense fallback occurred.

The default behavior remains automatic Björck-Pereyra with residual validation
and dense fallback. Solver thresholds and performance measurements belong in
source code and benchmark documentation.

## Terminal Quadrature Representation

QMOMKit terminal quadrature data is exposed through `CqmomQuadrature`.
Coordinates are always stored by phase-space dimension:

```text
abscissas[dimension][node]
```

Weights are global terminal weights, not branch-local conditional weights.

For fixed rectangular results:

```text
has_rectangular_storage() == true
weights()[flat_node_index]
abscissas()[dimension][flat_node_index]
```

The flat node index is the phase-space-order row-major index from
`PhaseSpaceHierarchy`. Rectangular accessors are valid:

```cpp
quadrature.weight(node_index)
quadrature.abscissa(dimension, node_index)
quadrature.node(node_index)
```

For adaptive branch-count results:

```text
has_rectangular_storage() == false
weights()[terminal_node]
abscissas()[dimension][terminal_node]
terminal_paths()[terminal_node]
```

The terminal arrays are deterministic traversal-order arrays. Rectangular
multi-index access is not valid for ragged adaptive output.

`ConditioningPath` records node indices by conditioning level:

```text
path[level] corresponds to conditioning_order[level]
```

For fixed rectangular storage, complete paths can be mapped to phase-space
node multi-indices. For adaptive ragged output, paths are exposed so callers
and diagnostics can identify the conditioning branch that produced each
terminal node.

## Adaptive Recursive CQMOM

Adaptive recursive CQMOM separates standalone QMOM controls from CQMOM branch
controls:

- `adaptive_branch_counts` is the CQMOM reconstruction policy that controls
  whether recursive CQMOM accepts branch-local node counts smaller than the
  configured hierarchy;
- when `adaptive_branch_counts` is enabled, recursive CQMOM automatically uses
  adaptive Wheeler for every marginal and conditional inversion.

The adaptive branch-count policy is available through the normal
`CqmomCase`/`CqmomState` workflow and through explicit `CqmomDenseCase` and
`CqmomLeanCase` workflows.

When adaptive branch counts are disabled, recursive CQMOM expects fixed
rectangular output storage. In that mode, adaptive Wheeler is rejected because
a branch with fewer nodes than requested would not match the rectangular node
hierarchy.

When adaptive branch counts are enabled:

- `node_counts[dimension]` is the requested maximum for that dimension;
- a successful branch-local Wheeler result may have fewer nodes than the
  configured maximum;
- different branches may have different child counts;
- terminal output uses ragged terminal arrays and conditioning paths;
- diagnostics record requested and actual node counts per branch.

Adaptive reduction is not a failure when adaptive branch counts are enabled
and Wheeler returns success. A branch whose actual node count is smaller than
the requested count records `adaptive_node_count_reduced = true` in the
public diagnostics.

This policy is useful for deterministic or nearly degenerate conditional
relations, near-boundary moment sequences, and high-dimensional cases where
some branches are effectively lower order than the configured maximum
hierarchy.

## Numerical Assumptions And Failure Handling

Recursive CQMOM assumes:

- multidimensional moments are raw moments of the number density function;
- `RawMomentTensor` extents and values are phase-space ordered;
- conditioning order is a permutation of phase-space dimensions;
- node counts are positive and owned by phase-space dimensions;
- all raw moments and numerical options are finite;
- dense raw moment tensor extents contain at least
  `2 * node_counts[dimension]` orders, even when adaptive branch counts are
  enabled;
- branch-local one-dimensional slices follow the same raw moment convention
  expected by Wheeler;
- weighted Vandermonde systems must be nonsingular under the configured
  singularity tolerance.

The implementation reports failures with explicit status values. Public
`CqmomStatus` values are:

- `success`: recursive construction completed and the returned quadrature was
  accepted.
- `invalid_input`: hierarchy, tensor, option, storage, or adaptive-mode
  requirements were violated.
- `qmom_failed`: a branch-local Wheeler inversion failed or returned an
  unacceptable node count.
- `conditional_moment_failed`: a branch-local Vandermonde solver
  reconstruction failed.

Internal recursive statuses map to public CQMOM statuses. In particular,
internal Wheeler failures are reported publicly as `qmom_failed`.

Failure handling is fail-fast. The recursive implementation stops at the first
failing branch and records diagnostics for the attempted branch. Diagnostic
entries include:

- conditioning level;
- parent path;
- phase-space dimension processed at that level;
- requested node count;
- actual node count;
- adaptive branch-count flags;
- mapped Wheeler/QMOM status;
- conditional-moment status;
- requested Vandermonde solver, selected solver, algorithm used, and dense
  fallback flag when available;
- Vandermonde residual, condition estimate, and warning flag when available;
- message.

Callers must check `CqmomResult::ok()` or `status` before using quadrature
data. A failed result can contain partially assembled diagnostic data, but
that data is not a valid reconstructed multidimensional quadrature.

### Bounded Coordinate Support

Recursive CQMOM accepts optional bounded intervals in phase-space dimension
order. An empty support vector means every coordinate is unconstrained, which
preserves the ordinary Wheeler reconstruction path. A nonempty vector must have
one entry per phase-space dimension; only dimensions with an interval use the
support-aware candidate path.

For a bounded coordinate, each branch-local `N`-node candidate is checked before
acceptance with:

```cpp
check_realizability(MomentSet(m_0, ..., m_(2*N-1)),
                    BoundedInterval{lower, upper},
                    RealizabilityOptions{});
```

Only support-realizable candidates proceed to Wheeler inversion. The inversion
is run for the candidate count, and the returned nodes must lie inside the
configured interval within the realizability tolerance. This is an acceptance
criterion, not a node clipping or projection step.

Adaptive branch-count reconstruction retries bounded candidates from the
current ceiling down to one node. When active branch-count metadata permits a
one-node growth trial, failed bounded realizability, failed inversion, or
out-of-interval nodes reject only that growth attempt and preserve the current
accepted branch count. Fixed-count reconstruction does not silently reduce a
bounded requested count; support infeasibility is reported as a reconstruction
failure.

Bounded-support diagnostics are attached to the affected CQMOM branch
diagnostic. They identify the phase-space coordinate, conditioning path and
level, interval, previous active count when available, attempted candidate
count, attempt type, candidate conditional moments, bounded realizability
status, effective tolerance, limiting matrix family and signed margin, accepted
fallback count, and any returned nodes and weights from attempted inversion.
The acceptance decision uses `check_realizability`; determinant or
Hankel-Hadamard/localizing determinant analysis remains diagnostic-only.

### Diagnostic Levels

CQMOM diagnostics balance numerical transparency with reconstruction cost. Some
diagnostic quantities, especially branch-local Vandermonde residual norms and
condition estimates, are useful for analysis and debugging but are not required
to construct the quadrature. The CQMOM diagnostic setting on the loaded case or
explicit representation case selects how much diagnostic work is performed:

`CqmomDiagnosticLevel` controls diagnostic retention and reporting only. It does
not alter Vandermonde solver selection, residual acceptance, dense fallback
behavior, inversion results, active node counts, quadrature values, transport
results, or success/failure status.

- `CqmomDiagnosticLevel::detailed` preserves full diagnostic reporting. It
  provides residual norms, condition estimates, and conditioning warnings for
  branch-local conditional-moment solves when those solves are performed. These
  values are informational diagnostics; reconstruction correctness is still
  determined by algorithm status and validated quadrature output.
- `CqmomDiagnosticLevel::basic` records structural branch diagnostics, statuses,
  path and level information, phase-space dimensions, requested and actual node
  counts, and adaptive branch-count information. Expensive Vandermonde residual
  and condition-estimate diagnostics are omitted.
- `CqmomDiagnosticLevel::none` minimizes success-path informational diagnostics.
  Failure statuses, failure messages, and available context for unsuccessful
  reconstruction are still preserved so callers can distinguish successful
  output from failed reconstruction.

Expensive diagnostic fields use optional values. Callers should treat absence as
"not computed" rather than as a numerical value. The public API defaults to
`detailed` to preserve diagnostic transparency unless the caller explicitly
selects a lower diagnostic level.

Configuration files control the same request field with:

```ini
[cqmom]
diagnostic_level = detailed|basic|none
```

## API Boundary

The public CQMOM API is declared in `<qmomkit/cqmom.hpp>`:

```cpp
qmomkit::CqmomResult compute_cqmom(const qmomkit::CqmomCase& case_);
qmomkit::CqmomResult compute_cqmom(const qmomkit::CqmomDenseCase& dense_case);
qmomkit::CqmomResult compute_cqmom(const qmomkit::CqmomLeanCase& lean_case);
```

Public CQMOM concepts include:

- `compute_cqmom`
- `CqmomCase`
- `CqmomState`
- `CqmomDenseCase`
- `CqmomLeanCase`
- `CqmomResult`
- `CqmomQuadrature`
- `CqmomDiagnostics`
- `CqmomStatus`
- `RawMomentTensor`
- `CqmomMomentLayout`
- `CqmomMomentData`
- `PhaseSpaceHierarchy`
- `ConditioningPath`

`CqmomCase` and `CqmomState` are the neutral public objects for the normal
CQMOM workflow. The case stores initialized moment data and reconstruction
policies; a state is created from the case when users choose CQMOM transport.
This workflow avoids exposing dense or lean moment storage to users.

`CqmomDenseCase` bundles the phase-space hierarchy, raw moment tensor,
conditional-moment numerical options, QMOM/Wheeler options, and the adaptive
branch-count policy. It is the advanced explicit dense representation workflow
for callers that intentionally use full rectangular tensor moment storage.

`CqmomLeanCase` bundles lean CQMOM moment data, the same
conditional-moment and QMOM/Wheeler options, and the same adaptive branch-count
and diagnostic controls. It is the public case type for callers that can
produce required CQMOM moments directly. It is the advanced explicit lean
representation workflow and provides direct control over the required CQMOM
moment representation. The dense `CqmomDenseCase` path internally selects the
same required moment data before recursive reconstruction.

`CqmomMomentLayout` and `CqmomMomentData` are public data-transfer types for
the lean CQMOM moment representation. They are not intended to expose
recursive tree state, branch-local slices, solver workspaces, or arbitrary
sparse tensor operations.

`CqmomResult` exposes the quadrature, public status, message, branch
diagnostics, and failed phase-space dimension. `CqmomQuadrature` exposes fixed
rectangular or adaptive ragged terminal quadrature storage through a common
public container.

The public `reconstruct_conditional_moments(...)` API is a supported
lower-level mathematical utility for one conditional-moment reconstruction
from a global raw moment tensor. It documents and exposes the weighted
Vandermonde convention, but it is not the public multidimensional CQMOM entry
point and it does not perform recursive quadrature assembly.

Implementation details include:

- `src/internal/recursive_quadrature.hpp`
- `compute_recursive_quadrature`
- `RecursiveQuadrature`
- `RecursiveQuadratureRequest`
- `RecursiveQuadratureResult`
- `IndexedMomentSet`
- `CqmomMomentView`
- branch-local moment slices;
- recursive tree nodes;
- branch-state traversal;

## Validation Expectations

Successful CQMOM reconstruction is validated primarily by reconstructing raw
moments from the recovered terminal quadrature:

```text
m_alpha ~= sum_j W_j product_r xi_(r,j)^alpha_r
```

Deterministic verification also checks terminal coordinates, global weights,
local branch weights, conditioning paths, branch counts, diagnostics, and
conditioning-order behavior where the analytical structure is known.

Current verification material covers:

- two-dimensional independent-product cases;
- deterministic affine conditional relations;
- branch-dependent conditional rules;
- non-identity conditioning orders;
- fixed rectangular CQMOM;
- adaptive ragged CQMOM;
- deeper high-dimensional recursive traversals;
- raw-moment reconstruction from terminal quadrature.

Changes to conditioning order, path mapping, branch-local slice layout,
adaptive acceptance, status propagation, or terminal storage are convention
changes and require matching unit and verification coverage.
