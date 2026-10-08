# QMOM Wheeler inversion

This page documents the current QMOMKit design for one-dimensional QMOM
inversion with the Wheeler algorithm [[Wheeler,
1974](https://www.doi.org/10.1216/RMJ-1974-4-2-287)].

For raw-moment, node, abscissa, weight, and ordering conventions, see
[numerical_conventions.md](numerical_conventions.md). This page describes how
those conventions are used by the Wheeler inversion implementation.

## Contents

- [Role of Wheeler Inversion](#role-of-wheeler-inversion)
- [QMOM Representation](#qmom-representation)
- [Algorithm Overview](#algorithm-overview)
- [Adaptive Wheeler](#adaptive-wheeler)
- [Optional Abscissa-Bound Acceptance](#optional-abscissa-bound-acceptance)
- [Numerical Assumptions And Failure Handling](#numerical-assumptions-and-failure-handling)
- [API Relationship](#api-relationship)
- [Validation Expectations](#validation-expectations)

## Role of Wheeler inversion

QMOM represents a one-dimensional number density function through a finite set
of raw moments:
$$m_k = \int \xi^k n(\xi) \ \mathrm{d}\xi$$

Many algorithms need an equivalent discrete quadrature:
$$m_k = \sum_i w_i x_i^k$$


Wheeler inversion maps a supplied moment sequence to quadrature abscissas
$x_i$ and weights $w_i$. For an $N$-node quadrature, the inversion uses the first
$2N$ raw moments: $m_0, m_1, ..., m_{2N-1}$.


The resulting quadrature rule is the local one-dimensional representation used
by public QMOM reconstruction and by CQMOM branch-local reconstructions. CQMOM
conditioning, branch paths, and terminal multidimensional weights are separate
topics for [cqmom_recursive.md](cqmom_recursive.md).

## QMOM representation

The one-dimensional input is a `qmomkit::MomentSet` whose entries are stored in
ascending moment order:

```text
moments[k] == m_k
```

The output quadrature contains:

- abscissas: one coordinate value per node;
- weights: one quadrature contribution per node;
- status and message fields that describe success or numerical failure.

Returned abscissas are sorted in ascending order. Weights are permuted with the
same ordering:

```text
weights[i] corresponds to abscissas[i]
```

The source order of any exact nodes used to generate a moment sequence is not
part of the input convention. Raw moments are sums over nodes, and Wheeler
returns the deterministic sorted convention.

The zeroth moment is preserved as the total quadrature weight on successful
rules:
$$\sum_i w_i = m_0$$

With `automatic_scaling = false`, QMOMKit passes physical raw moments to
Wheeler unchanged. With optional automatic reconstruction scaling enabled, the
reconstruction wrapper supplies coordinate-scaled raw moments as numerical
preconditioning; Wheeler still operates on raw moments and does not perform
centering or $m_0$ normalization. A moment sequence with $m_0 = 1$ is treated
as a unit-number-density sequence, not as a separate probability convention.

## Algorithm overview

QMOMKit uses the Wheeler recurrence to construct recurrence coefficients for
the orthogonal polynomial sequence associated with the supplied raw moments.
Conceptually, the inversion has four main steps:

1. Validate the case/input, options, moment count, moment finiteness, and zeroth
   moment.
2. Construct recurrence coefficients from the moment sequence with the
   sigma-table recurrence.
3. Build the symmetric tridiagonal Jacobi matrix defined by those recurrence
   coefficients.
4. Diagonalize the Jacobi matrix, then form abscissas and weights from the
   eigenpairs.

For a Jacobi matrix $J$, the quadrature abscissas are the eigenvalues:
$$x_i = \lambda_i(J)$$


The weights are computed from the first component of the normalized
eigenvectors:
$$w_i = m_0 \, v_{0,i}^2$$


The multiplication by $m_0$ preserves the total number density represented by
the supplied moment sequence.

The implementation uses QMOMKit's internal small dense linear-algebra
utilities for the symmetric eigenproblem. External numerical backends such as
Eigen, BLAS, or LAPACK are not part of the current core dependency model. The
recurrence tables, Jacobi matrix construction, and eigensolver plumbing live
in implementation files and are not exposed as public API objects.

## Adaptive Wheeler

The adaptive Wheeler alogrithm [[Marchisio & Fox,
2013](https://www.doi.org/10.1017/CBO9781139016599)] allows a case for up to
$N$ nodes to return a smaller accepted quadrature rule when the supplied
moments are effectively supported by fewer distinct nodes. This is important
near the boundary of realizable moment space, where a fixed-node inversion can
encounter nearly duplicate nodes, tiny weights, or degenerate recurrence data.

Adaptive behavior is controlled by `QmomOptions` or `WheelerOptions`:

- `adaptive`: enables adaptive node reduction; on by default.
- `rmin`: minimum accepted ratio between the smallest and largest weight.
- `eabs`: dimensionless minimum accepted relative node separation.
- `cutoff`: small positive threshold for near-zero mass.

The adaptive algorithm starts from the requested maximum node count. It checks
the Wheeler recurrence data for degeneracy, may reduce the active node count,
then attempts candidate quadratures from the active count down to one node.

The recurrence is degenerate at count $k$ when $\sigma_{k,k}$, which equals
$m_0 \prod_{i=1}^{k-1} b_i$, is zero within its own rounding. Adaptive Wheeler bounds
that rounding by running the $\sigma$ recurrence again with every term in
absolute value,
$M_{k,l} = M_{k-1,l+1} + |a_{k-2}|\,M_{k-1,l} + |b_{k-2}|\,M_{k-2,l}$ from
$M_{1,l} = |m_{l-1}|$, and treats $\sigma_{k,k}$ as zero when
$|\sigma_{k,k}| \le 100\,\varepsilon\,M_{k,k}$. Both sides scale as $m_0$
times the coordinate to the power $2(k-1)$, so the node count depends neither
on the unit of the coordinate nor on the number density, and the bound grows
with the cancellation of a distribution far from the origin. In recursive
CQMOM, conditional moments come from a solve whose rounding is relative to the
scale $s$ of their coordinate, $\sqrt{m_2/m_0}$ of its marginal, rather than to
their own magnitudes; for them the bound has the floor $m_0\,s^{2(k-1)}$.
A multi-node adaptive rule is accepted only when both the `rmin` weight-ratio
and the `eabs` node-separation criteria above pass. If no multi-node candidate
passes them, the current implementation can still return a successful reduced
one-node quadrature when that one-node rule is mathematically defined.
Therefore, `success` means the returned rule was accepted; it does not
guarantee that the originally requested node count was retained.

## Optional abscissa-bound acceptance

`AbscissaBounds` adds an optional acceptance condition around ordinary Wheeler
inversion. It is not a bounded Wheeler variant and does not change recurrence,
Jacobi, or eigensolver mathematics. A fixed candidate succeeds only when its
returned nodes lie in the configured interval within the shared `1e-12`
abscissa-bound tolerance. Adaptive bounded reconstruction evaluates explicit
candidate counts from the requested count down to one and accepts the highest
admissible exact-count candidate. An out-of-interval node is never clipped or
projected.

The same internal one-dimensional selection is used by direct QMOM and by
recursive CQMOM branches. `AbscissaBounds` remains separate from
`BoundedInterval`, which is the input type for bounded-interval moment
realizability analysis.

If a mathematically defined one-node fallback is available, the result may
contain that fallback with either success or a non-success status, depending
on the condition that produced it. Callers must use `ok()` or the explicit
status value, not vector non-emptiness, to determine whether an inversion
succeeded and whether the returned node count is acceptable for the caller's
workflow.

For a candidate with more than two nodes, the separation ratio is the minimum
pairwise abscissa distance divided by the maximum pairwise distance. For a
two-node candidate, it is the node distance divided by the larger abscissa
magnitude. Thus `eabs` is dimensionless for every multi-node candidate. A
two-node candidate whose abscissas both equal zero is degenerate and does not
pass the separation criterion.

Adaptive Wheeler is a one-dimensional inversion policy. In recursive CQMOM, the
`adaptive` option sets it for every marginal and conditional inversion and also
decides whether branch node counts adapt.

## Numerical assumptions and failure handling

The Wheeler implementation assumes:

- moments are raw moments of the number density function;
- one-dimensional moments are supplied in ascending order;
- an $N$-node case provides at least $2N$ moments;
- all moments consumed by the requested reconstruction and all numerical
  options are finite;
- $m_0$ is the total number density and must not be negative;
- recurrence coefficients must remain finite;
- recurrence coefficients that imply negative Jacobi off-diagonal squares
  indicate a non-realizable moment sequence;
- abscissas and weights must remain finite after diagonalization.

Cases with zero nodes, too few moments, non-finite moments, invalid options,
negative $m_0$, near-zero $m_0$, non-finite recurrence data, negative
recurrence data, or non-finite eigen output are reported through status
values. Adaptive rejection of a multi-node candidate is not by itself a
failure if the implementation reaches and accepts a reduced one-node rule.

Current status values are:

- `success`: the returned quadrature rule is accepted, possibly at a lower node
  count after adaptive retries.
- `invalid_input`: the case/input is structurally invalid or contains invalid
  numerical options.
- `not_realizable`: the supplied moments fail a realizability or recurrence
  check.

Non-success statuses may still carry diagnostic fallback data. This behavior is
intentional and supports downstream diagnostics or controlled failure handling.
Fallback values may not preserve the original moment sequence and must not be
treated as a successful quadrature rule. Callers should check the status before
using abscissas or weights in reconstruction, phase-space evolution,
validation, or output.

Numerical conditioning is most sensitive when moments are high order, nodes
are close together, weights differ by many orders of magnitude, or the
sequence lies near the boundary of realizable moment space. Changes to
thresholds, recurrence checks, eigenvalue sorting, or fallback behavior are
convention changes and require matching unit and verification coverage.

## API relationship

The preferred public QMOM entry points are declared in `<qmomkit/qmom.hpp>`:

```cpp
qmomkit::QmomResult compute_qmom(
    const qmomkit::MomentSet& moments,
    std::size_t node_count,
    const qmomkit::QmomOptions& options = {});

qmomkit::QmomResult compute_qmom(const qmomkit::QmomCase& qmom_case);
```

`QmomCase` bundles:

- the raw moment sequence;
- the requested node count;
- the QMOM numerical options.

`QmomResult` exposes:

- sorted abscissas;
- matching weights;
- `QmomStatus`;
- a diagnostic message;
- `ok()` for success checking.

The lower-level Wheeler entry point is declared in `<qmomkit/wheeler.hpp>`:

```cpp
qmomkit::Quadrature1D compute_quadrature(
    const qmomkit::MomentSet& moments,
    std::size_t max_nodes,
    const qmomkit::WheelerOptions& options = {});
```

`compute_quadrature` is a public numerical utility for direct
one-dimensional Wheeler inversion. `compute_qmom` wraps this routine and maps
`QuadratureStatus` to `QmomStatus`. Recursive CQMOM also uses Wheeler
inversion internally for branch-local one-dimensional moment sequences.

The following are public API concepts:

- `MomentSet`
- `Quadrature1D`
- `QuadratureStatus`
- `WheelerOptions`
- `compute_quadrature`
- `QmomOptions`
- `QmomCase`
- `QmomResult`
- `QmomStatus`
- `compute_qmom`

The following are implementation details:

- sigma-table storage;
- recurrence-coefficient work arrays;
- Jacobi matrix assembly;
- eigenvalue/eigenvector workspace;
- adaptive scan and fallback mechanics inside `src/wheeler.cpp`;
- internal dense linear-algebra helpers under `src/internal/`.

Public headers should not expose the internal recurrence or eigensolver data
structures unless QMOMKit deliberately adds a supported diagnostic API.

## Validation expectations

Wheeler changes should preserve the core invariant that successful quadrature
reconstructs the supplied raw moments to the documented tolerance:
$$m_k \simeq \sum_i w_i x_i^k$$


Existing coverage includes:

- public API wrapper behavior;
- one-node, two-node, and three-node analytical cases;
- scaled weights and scaled abscissas;
- source-node order invariance;
- adaptive reduction for degenerate, near-duplicate, and tiny-weight cases;
- invalid input and invalid option handling;
- non-realizable moment sequences;
- deterministic verification against known quadrature rules.

Any Wheeler behavior change must update the relevant unit and verification
tests before QMOM, CQMOM, phase-space evolution, configuration parsing, or
validation workflows rely on it.