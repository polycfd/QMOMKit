# Numerical conventions

All code and tests share these conventions; changing one changes a numerical
contract (see [AGENTS.md](../../AGENTS.md)).

## Moments

QMOMKit works with raw moments of the number density, defined in
[Concepts](../user_new/3_concepts.md#raw-moments). $m_0$ and $M_{0,\ldots,0}$
are the total number density, and no algorithm divides by them implicitly.

A moment generated from a distribution scales a normalized shape $p$ by a
number density $n_0$:

```math
n(\xi) = n_0\, p(\xi), \qquad m_k = n_0 \int \xi^k p(\xi)\,\mathrm{d}\xi, \qquad m_0 = n_0 .
```

Mixtures and independent products follow from this rule; see
[Distributions](7_distributions.md).

Moments have integer orders. Generalized moments with real exponents are
derived from a quadrature, never an input to reconstruction or evolution.

## Indices and storage

| Container | Indexed by | Storage |
| --- | --- | --- |
| `MomentSet` | order $k$ | `values()[k]` is $m_k$; $N$ nodes need $m_0, \ldots, m_{2N-1}$. |
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
required layout from `create_cqmom_moment_layout` is a subset of that box,
enumerated one conditioning level at a time and stored in phase-space order.
Its order depends only on the hierarchy, so a host creates a layout once and
replaces the values. [CQMOM](4_cqmom.md) derives the layout.

`MultiIndexTable` and `ConditioningPathTable` store equal-length rows
contiguously, row-major, in one allocation. A row is read as a
`MultiIndexView` or `ConditioningPathView`, which must not outlive its table.

## Dimension order and conditioning order

The phase-space dimension order is the order in which the coordinates are
declared. It indexes moments, node counts, node multi-indices, abscissa arrays,
abscissa bounds, and per-coordinate realizability policies, and it orders the
coordinates within each node of the packed, node-major velocity-callback
vectors.

`PhaseSpaceHierarchy` holds the node counts in that order, together with
`conditioning_order`, a permutation in which `conditioning_order[l]` is the
dimension reconstructed at level `l`. The default is the identity. The
constructor rejects a sequence that is not a permutation and a zero node count.

The conditioning order changes the reconstruction, never a storage order or the
meaning of a node count. Only two things are ordered by conditioning level:
conditioning paths, whose entry `l` is a node index in dimension
`conditioning_order[l]`, and the internal branch layouts and slices of the
recursive traversal. Everything public is ordered by dimension.

## Quadrature storage

**QMOM.** Abscissas are ascending, and the weights follow the same order.

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

Reconstruction is adaptive by default (`QmomOptions::adaptive = true`): node
counts are maxima, and Wheeler returns fewer nodes where the moments do not
support more. In CQMOM the one flag governs every marginal and conditional
inversion and every branch count. A fixed reconstruction returns exactly the
configured count in every inversion, or fails. Criteria:
[QMOM](3_qmom.md); branches: [CQMOM](4_cqmom.md).

## Tolerances

Tolerances are relative, so a decision depends neither on the unit of a
coordinate nor on $m_0$. Examples are the Wheeler roundoff bounds, the
weight ratio `rmin`, the realizability test on standardized moment matrices,
the relative pivot floor of conditional solves, and the margin of abscissa
bounds. Despite its name, `eabs` is relative: it compares the smallest node
separation with the largest, or for two nodes with their larger magnitude.

Few thresholds are absolute. `QmomOptions::cutoff` applies to $m_0$, in units of
number density. Below it, reconstruction reports `not_realizable`, and
`state_candidate` reports `empty`. The error control of the embedded Runge-Kutta
schemes of node advection has an absolute tolerance in the units of the
coordinates.

Each value is given once, on the page that owns it:
[QMOM](3_qmom.md), [Linear algebra](5_linear_algebra.md),
[Realizability](6_realizability.md), and
[Phase-space evolution](8_phase_space_evolution.md).

## Supports and abscissa bounds

The two concepts are distinct and never derived from each other.

A **support** (`MomentSupport::real_line`, `nonnegative`, or
`bounded_interval` with an explicit `BoundedInterval`) is an argument of a
realizability check or policy and does not constrain reconstruction
([Realizability](6_realizability.md)).

**Abscissa bounds** (`AbscissaBounds`, or `CqmomAbscissaBounds` per coordinate)
are acceptance criteria for reconstructed abscissas. A candidate rule with a
node outside the bounds, beyond a rounding margin, is rejected. An adaptive
reconstruction then tries fewer nodes, and a fixed one fails. Nodes are never
clipped or projected. Bounds exist only when the caller gives them. In a
`.qmom` file, they are the keys `phase_space.support_lower` and
`phase_space.support_upper`, which set abscissa bounds despite their names.
They are never inferred from a realizability support, a coordinate name, or the
moments. The margin is in [QMOM](3_qmom.md).

## Units and scales

Coordinate names are metadata: they carry no unit, support, or meaning.
QMOMKit never rescales moments or coordinates for the caller; results are in
the caller's units. Internal kernels may scale the systems they solve, invisibly
in the results ([Linear algebra](5_linear_algebra.md)).

The only internal coordinate scales bound the rounding of CQMOM conditional
moments in adaptive reduction; see [QMOM](3_qmom.md#coordinate-scales).

## Terms

Terms the [user guide](../user_new/3_concepts.md#terms) does not define:

| Term | Meaning |
| --- | --- |
| active node count | Nodes that a branch or a QMOM state uses; at most the configured count. |
| canonical moments | Moments recomputed from an accepted quadrature; a state stores these, not the moments it was given. |
| conditional moments | Moments of later coordinates conditioned on a node of earlier ones, from a conditional Vandermonde solve. |
| conditioning level | Position `l` in the conditioning order. |
| conditioning path | `ConditioningPath`: the node index chosen at each level from the root to a node. |
| configured node count | Node count of a dimension in the hierarchy: the maximum of every branch. |
| full tensor | `RawMomentTensor`: every moment in a rectangular box; an input or output format that cases and states never hold. |
| local weight | Weight from one branch inversion; terminal weights are products of local weights. |
| marginal moments | Moments with every order zero except one. |
| parent path | Conditioning path of the node on which a branch is conditioned; its length is the level of the branch. |
| prepared metadata | `internal::PreparedCqmomMetadata`: the immutable branch structure that a prepared case shares with its states and derived cases. |
