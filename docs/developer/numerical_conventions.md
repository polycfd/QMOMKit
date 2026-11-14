# Numerical conventions

This page fixes the conventions every algorithm and test shares: moments,
indices, orderings, quadrature storage, tolerances, supports, bounds and
units. Changing one changes a numerical contract
([AGENTS.md](../../AGENTS.md#numerical-contracts)).

## Moments

QMOMKit works with raw moments of the number density, defined in
[Concepts](../user/3_concepts.md#raw-moments). $m_0$ and $M_{0,\ldots,0}$
are the total number density, and no algorithm divides by them implicitly.

Moments generated from a distribution are those of a normalized shape scaled
by a number density $n_0$, so $m_0 = n_0$; see [Distributions](distributions.md).

Moments have integer orders. Generalized moments with real exponents are
derived from a quadrature, never an input to reconstruction or evolution.

## Indices and storage

| Container | Indexed by | Storage |
| --- | --- | --- |
| `MomentSet` | order $k$ | `values()[k]` is $m_k$; $N$ nodes of the classical family need $m_0, \ldots, m_{2N-1}$, and each family states its count ([QMOM families](qmom_families.md)). |
| `MultiIndex` | phase-space dimension | `std::vector<std::size_t>`, one order or node index per dimension. |
| `RawMomentTensor` | multi-index $\alpha$ | `values[flatten_multi_index(alpha, extents)]`. |
| `CqmomMomentData` | position in its layout | `values[i]` is the moment of `layout.required_indices` row `i`. |

`flatten_multi_index` is row-major in phase-space order, with the last dimension
varying fastest. For extents $n_0, \ldots, n_{D-1}$:

```math
\mathrm{flat}(\alpha) = \bigl(\cdots(\alpha_0 n_1 + \alpha_1) n_2 + \cdots\bigr) n_{D-1} + \alpha_{D-1} .
```

It throws for a rank mismatch, a component outside its extent, or an overflow.

Maximum orders are inclusive: a maximum order $K$ gives an extent $K + 1$.
`PhaseSpaceHierarchy::moment_extents()` is $2N_d$ in dimension $d$. The
required layout from `create_cqmom_moment_layout` is a subset of that box.
Its entries are listed level by level of the conditioning order, and each is a
multi-index in phase-space order. The layout depends only on the hierarchy;
[CQMOM](cqmom.md) derives it.

`MultiIndexTable` and `ConditioningPathTable` store equal-length rows
contiguously, row-major, in one allocation. A row is read as a
`MultiIndexView` or `ConditioningPathView`, which must not outlive its table.

## Dimension order and conditioning order

The phase-space dimension order is the order in which the coordinates are
declared. It indexes moments, node counts, node multi-indices, abscissa arrays,
abscissa bounds, and per-coordinate realizability policies. In the packed
velocity vectors, the coordinates within each node follow it too.

`PhaseSpaceHierarchy` holds the node counts in that order, together with
`conditioning_order`, a permutation in which `conditioning_order[l]` is the
dimension reconstructed at level `l`. The default is the identity. The
constructor rejects no node counts, a zero node count, a total node count
beyond `std::size_t`, and an order that is not a permutation.

The conditioning order changes the reconstruction, never the meaning of an
index or a node count: every per-coordinate index, from multi-index components
to abscissa arrays, bounds and policies, is in phase-space order. The
conditioning order sets only sequences: conditioning paths, whose entry `l` is
a node index in dimension `conditioning_order[l]`; the internal per-branch data
of the recursion ([CQMOM](cqmom.md)); the order of the required layout's
entries; and the terminal order of an adaptive quadrature.

## Quadrature storage

**QMOM.** Abscissas are ascending, and the weights follow the same order.
The same holds for every inversion inside CQMOM: node index $i$ of a branch is
its $i$-th smallest abscissa. Conditioning paths, the terminal order, the
active-count records and the output rows all rely on this order.

**CQMOM.** A quadrature consists of terminal nodes. `weights()[j]` is the global weight of
terminal node `j`, the product of the local weights along its path.
`abscissas()[d][j]` is its coordinate in dimension `d`, and
`terminal_paths()` row `j` is its conditioning path. Two layouts exist:

| | Fixed (`adaptive = false`) | Adaptive (default) |
| --- | --- | --- |
| `has_rectangular_storage()` | `true` | `false`, even when no branch reduced |
| Terminal order | Row-major over `node_counts()` in phase-space order, as `flatten_node_index` | Lexicographic in the conditioning path |
| `weight(node_index)`, `abscissa(d, node_index)`, `node(...)` | Valid | Throw `std::logic_error` |
| `CqmomState::active_branch_counts()` | Empty | One record per branch |

Active branch counts (`CqmomActiveBranchCount`) are sorted by parent path, so
every branch follows its parent. Code for both layouts iterates over
`weights()`, `abscissas()` and `terminal_paths()`.

## Adaptive and fixed reconstruction

One flag, `QmomOptions::adaptive` (default `true`), governs every inversion; in
CQMOM, every marginal and conditional inversion and every branch count. An
adaptive reconstruction returns at most the family's maximum node count,
`maximum_node_count()`, which is $N + 1$ for `hyperbolic` and $N_G$ for
`generalized`. A fixed one returns exactly that count in every inversion, or
fails, except that a fixed `generalized` reconstruction returns the classical
$N$ nodes of moments on the boundary of moment space. Criteria:
[QMOM](qmom.md); branches: [CQMOM](cqmom.md).

## Tolerances

Numerical decisions use relative tolerances, so they depend neither on the
unit of a coordinate nor on $m_0$: the Wheeler rounding bounds, `rmin` and
`eabs` (relative despite its name), the realizability test on standardized
moment matrices, the singularity tests of conditional solves, and the margin of
abscissa bounds. The conditional solves scale their systems by powers of two,
so their decisions are exactly unit-free only for a change of unit by a power of
two ([Linear algebra](linear_algebra.md)).

A few thresholds are absolute: `QmomOptions::cutoff` on $m_0$, in units of
number density; the absolute tolerance of adaptive node advection, in units of
the coordinates; and the step sizes of the integrators (`initial_step`,
`min_step`, `max_step`, `min_time_step`), in units of the time argument.

Each value is given once, on the page that owns it:
[QMOM](qmom.md), [Linear algebra](linear_algebra.md),
[Realizability](realizability.md), and
[Phase-space evolution](phase_space_evolution.md).

## Supports and abscissa bounds

The two concepts are distinct and never derived from each other.

A **support** (`MomentSupport::real_line`, `nonnegative`, or
`bounded_interval` with an explicit `BoundedInterval`) is an argument of a
realizability check or policy and does not constrain reconstruction
([Realizability](realizability.md)).

**Abscissa bounds** (`AbscissaBounds`, or `CqmomAbscissaBounds` per coordinate)
are acceptance criteria for reconstructed abscissas. A candidate quadrature with
a node outside the bounds, beyond a rounding margin, is rejected. An adaptive
reconstruction then tries fewer nodes, and a fixed one fails. Nodes are never
clipped or projected. Bounds exist only when the caller gives them. In a `.qmom`
file, they are the keys `phase_space.support_lower` and
`phase_space.support_upper`, which set abscissa bounds despite their names. They
are never inferred from a realizability support, a coordinate name, or the
moments. The margin is in [QMOM](qmom.md).

## Units and scales

Coordinate names are metadata: they carry no unit, support, or meaning.
QMOMKit never rescales moments or coordinates for the caller; results are in
the caller's units. Internal kernels may scale the systems they solve, invisibly
in the results ([Linear algebra](linear_algebra.md)).

The only internal coordinate scales bound the rounding of CQMOM conditional
moments in adaptive reduction; see [QMOM](qmom.md#coordinate-scales).

## Terms

Terms the [user guide](../user/3_concepts.md#terms) does not define:

| Term | Meaning |
| --- | --- |
| active node count | Nodes that a branch or a QMOM state uses; at most the configured count, or for a QMOM state its family's maximum. |
| base node count | `QmomResult::base_node_count` and `QmomState::base_node_count()`, the number of nodes $k\le N$ of the classical quadrature on which a QMOM quadrature rests. It is the number of nodes for `classical`, $N$ for the $N+1$ nodes of `hyperbolic`, $k$ for the $N_G$ nodes of `generalized` continued from $k$ nodes, and 0 for a failed result. A value below $N$ shows an adaptive reduction in every family. It is also the node count the evolution passes when it reconstructs a state at its own count ([Phase-space evolution](phase_space_evolution.md#node-counts)). |
| canonical moments | Moments recomputed from an accepted quadrature; a state stores these, not the moments it was given. |
| ceiling | The largest node count a branch may return in one reconstruction, apart from the growth trial: its recorded active count when a state's counts are passed in, otherwise its configured count. |
| child case | Internal layout data for the children of a branch with a given node count ([CQMOM](cqmom.md#branch-layouts)). |
| conditional moments | Moments of later coordinates conditioned on a node of earlier ones, from a conditional Vandermonde solve. |
| conditioning level | Position `l` in the conditioning order. |
| conditioning path | `ConditioningPath`: the node index chosen at each level from the root to a node. |
| configured node count | Node count of a dimension in the hierarchy: the maximum of every branch. For QMOM, the case's `node_count` $N$; `maximum_node_count()` is the family's maximum. |
| full tensor | `RawMomentTensor`: every moment in a rectangular box; an input or output format that cases and states never hold. |
| growth | A branch trying one more node than its active count, allowed only in the post-source reconstruction. |
| local weight | Weight from one branch inversion; terminal weights are products of local weights. |
| marginal moments | Moments with every order zero except one. |
| maximum node count | `QmomState::maximum_node_count()`, the largest number of nodes a QMOM state's quadrature can have. Its case fixes it at $N$ for `classical`, $N+1$ for `hyperbolic` and $N_G$ for `generalized`. |
| parent path | Conditioning path of the node on which a branch is conditioned; its length is the level of the branch. |
| prepared metadata | `internal::PreparedCqmomMetadata`: the immutable branch structure that a prepared case shares with its states and derived cases. |
| required layout | The user guide's *layout*, `CqmomMomentLayout`: the moments a hierarchy's reconstruction reads, in a fixed order. |
| slice | The moments one CQMOM branch reconstructs from, indexed in conditioning order. |
