# Internal Linear Algebra

This document defines the architecture, mathematical conventions, numerical
behavior, and current QMOMKit roles of the internal dense linear-algebra kernels
and the conditional Vandermonde adapter.

The APIs described here live under `src/internal/` and are implementation
details, not public APIs. This page is for QMOMKit developers who need to
maintain the numerical kernels or their CQMOM integration. User-facing
configuration syntax belongs in [configuration.md](../user/configuration.md);
recursive CQMOM workflow and branch traversal details belong in
[cqmom_recursive.md](cqmom_recursive.md).

## Scope And Design Principles

The internal linear-algebra layer provides dependency-minimal dense numerical
building blocks for small systems whose mathematical inputs are owned by higher
level QMOMKit algorithms. It prioritizes mathematical correctness, explicit
matrix conventions, clear separation between numerical kernels and adapter
policy, and algorithms proportional to the small dense systems used by CQMOM.

Solver failures remain visible. QMOMKit does not silently change the
mathematical problem by substituting pseudoinverse, truncated, regularized, or
minimum-norm recovery paths when a configured solve is not accepted.

This layer does not own phase-space hierarchy, configuration parsing, recursive
CQMOM traversal, application physics, or public workflow design. It also does
not provide a general-purpose matrix package.

## Internal Component Map

### Numerical Kernels

The numerical kernels own their algorithms, supported matrix shapes,
singularity or numerical-rank detection, and direct solve semantics:

```text
solve_dense(...)
solve_bjorck_pereyra(...)
factorize_qr(...) / solve_qr(...)
factorize_svd(...) / solve_svd(...)
```

`solve_dense(...)`, `factorize_qr(...)`, `solve_qr(...)`, `factorize_svd(...)`,
and `solve_svd(...)` are generic dense kernels. `solve_bjorck_pereyra(...)` is a
structured, unweighted transposed Vandermonde kernel.

### Conditional Vandermonde Adapter

The conditional Vandermonde adapter owns the weighted production system used by
CQMOM conditional reconstruction:

```text
make_conditional_vandermonde_matrix(...)
solve_conditional_vandermonde_system(...)
```

The adapter constructs the weighted matrix, resolves the selected
`VandermondeSolver` policy, dispatches to the numerical kernels, performs
Björck-Pereyra (BP) weight recovery, validates residuals when required, applies
BP-to-dense fallback, records performance counters, and returns solver
diagnostics.

### Recursive CQMOM Integration

Recursive CQMOM constructs branch-local conditional systems, propagates the
configured `VandermondeSolver` policy to the adapter, maps failed adapter solves
to CQMOM status values, and records branch diagnostics. It does not implement
the dense, BP, QR, or SVD algorithms and does not perform numerical-rank-based
recovery. See [cqmom_recursive.md](cqmom_recursive.md) for the recursive
construction.

## Common Matrix And Solve Conventions

`DenseMatrix` is the shared internal dense matrix representation. It stores
dimensions and row-major values for local numerical kernels; the storage
convention is useful for maintainers but is not part of a public matrix API.

Solves may contain multiple right-hand-side columns. For square systems, the
returned solution has one row per unknown and one column per right-hand side.
Failure results must be checked before using returned solution data; data
attached to failed results is diagnostic.

The dense Gaussian-elimination and BP kernels solve square systems only. The QR
and SVD kernels support square and tall systems with `m >= n`; tall accepted
systems produce full-column-rank least-squares solutions. Conditional
Vandermonde systems used in production CQMOM reconstruction are square.

## Solver Comparison

| Solver | Mathematical system | Type | Supported shapes | Multiple RHS | Rank diagnostics | Current Vandermonde policy | Fallback role |
|---|---|---|---|---|---|---|---|
| Dense Gaussian elimination | \(AX=B\) | Generic dense | Square | Yes | No; reports singular-pivot failure | `dense` | Dense fallback target for fallback-capable BP |
| Björck-Pereyra | \(V^T Y=B\) | Structured, unweighted Vandermonde | Square Vandermonde | Yes | No; reports node-separation failure | `bjorck_pereyra`; first backend for `automatic` | May trigger dense fallback through adapter policy |
| Pivoted Householder QR | \(AP=QR\), solve \(AX\approx B\) | Generic dense | Square and tall, \(m\ge n\) | Yes | Numerical rank | `pivoted_qr` | None |
| One-sided Jacobi SVD | \(A=U\Sigma V^T\), solve \(AX\approx B\) | Generic dense | Square and tall, \(m\ge n\) | Yes | Singular values and numerical rank | `svd` | None |

The table summarizes capabilities and ownership. It is not a benchmark or an
accuracy ranking.

## Dense Gaussian-Elimination Solver

### Purpose

`solve_dense(...)` is the simple generic square dense solver. It provides a
baseline for dense square systems, the explicit `dense` Vandermonde policy, and
the dense fallback target for fallback-capable BP policies.

### Mathematical Problem

The dense solver computes `X` in

\[
A X = B,
\qquad
A \in \mathbb{R}^{n \times n}.
\]

`B` may contain one or more right-hand-side columns.

### Algorithm

The implementation uses Gaussian elimination with row pivoting, followed by
triangular back-substitution. All right-hand-side columns are advanced through
the same elimination, so multiple RHS columns are solved simultaneously.

### Supported Shapes And Right-Hand Sides

The dense solver supports square matrices only. It accepts multiple
right-hand-side columns. It does not support rectangular least-squares problems.

### Numerical Characteristics

The solver is general-purpose and does not exploit Vandermonde structure. It
does not produce a numerical-rank decomposition or singular values. Its main
role is a straightforward independent baseline and a robustness reference for
the conditional Vandermonde adapter.

### Failure Semantics

At elimination step `k`, the pivot is the largest absolute entry in column `k`
over rows `k..n-1`. The solve fails when the selected pivot satisfies

\[
\max_{i \ge k} |a_{ik}|
\le
\tau_{\mathrm{sing}}\max(1,\|A\|_\infty),
\]

where the infinity norm is computed from the original matrix passed into the
solve. Failure returns the caller-provided singular message and no accepted
solution.

### Current QMOMKit Role

The dense solver is used by the explicit `dense` Vandermonde policy and as the
dense fallback target for `bjorck_pereyra_with_dense_fallback`. The public
one-level conditional-moment utility also solves its weighted system directly
with the dense solver.

### Limitations And Non-Goals

The dense solver is square-only. It does not report numerical rank, compute
least-squares solutions, regularize singular systems, or expose a public matrix
API.

## Björck-Pereyra Solver

### Purpose

`solve_bjorck_pereyra(...)` exploits Vandermonde structure for the conditional
Vandermonde pathway. It is the structured kernel used before optional dense
fallback in the default CQMOM policy.

### Mathematical Problem

For nodes `x_i`, define

\[
V_{i,p}=x_i^p.
\]

The BP kernel solves the unweighted transposed Vandermonde system

\[
V^T Y = B.
\]

Rows of `B` correspond to powers, rows of `Y` correspond to nodes, and columns
of `B` are independent right-hand sides. The kernel is unweighted; quadrature
weights belong to the conditional Vandermonde adapter.

### Algorithm

The implementation applies structured \(\Theta(n^2 r)\) Vandermonde transformations
equivalent to divided-difference and Newton-basis operations. It avoids forming
or factoring a generic dense matrix and applies the same structured operations
to every right-hand-side column, where \(n\) is the Vandermonde size and \(r\)
is the right-hand-side width. For fixed \(r\), the kernel work is
\(\Theta(n^2)\). Optional internal instrumentation counts the executed
transform updates exactly as \(n(n-1)r\); this is a deterministic work-structure
check, not a wall-clock performance guarantee.

### Supported Shapes And Right-Hand Sides

The BP kernel supports square Vandermonde systems with one node per unknown and
one RHS row per power. It supports multiple right-hand-side columns. It is not a
generic rectangular solver.

### Numerical Characteristics

BP exploits exact Vandermonde structure and avoids the overhead of a generic
dense factorization. Its behavior depends on the structured problem, the node
values, and node ordering. It is sensitive to node separation and does not
provide a general numerical-rank estimate.

### Failure Semantics

The kernel rejects empty systems, inconsistent dimensions, invalid tolerances,
non-finite inputs, non-finite outputs, and insufficiently separated nodes. A
node-separation denominator is treated as singular when

\[
|x_i-x_j|
\le
\tau_{\mathrm{sing}}\max(1,|x_i|,|x_j|).
\]

The kernel reports this as a transposed Vandermonde singular failure. It does
not perform weighted residual validation or fallback.

### Current QMOMKit Role

BP is used by the explicit `bjorck_pereyra` policy, as the first backend for
`automatic`, and as the first backend for
`bjorck_pereyra_with_dense_fallback`.

### Limitations And Non-Goals

The BP kernel is Vandermonde-specific, unweighted, and square-only. It does not
construct conditional systems, divide by quadrature weights, provide rank
diagnostics, or own fallback policy.

## Pivoted Householder QR Solver

### Purpose

`factorize_qr(...)` and `solve_qr(...)` provide a generic square/tall dense
solver with numerical-rank detection. In the conditional Vandermonde adapter,
QR is an explicit direct backend for the complete weighted system.

### Mathematical Problem

The factorization is

\[
A P = Q R,
\]

for `A` with `m >= n`. The stored permutation satisfies

```text
column_permutation[pivoted_column] = original_column.
```

### Algorithm

The implementation uses Householder reflectors with column pivoting. Solves
apply \(Q^T\) implicitly to the right-hand side and then use triangular
back-substitution in the pivoted coordinates before unpermuting the solution.
The factorization stores reflector data compactly in the factor matrix and a
coefficient vector.

### Supported Shapes And Right-Hand Sides

QR supports square and tall matrices with `m >= n`. Square full-column-rank
systems are ordinary linear solves. Tall full-column-rank systems produce
least-squares solutions. Factorizations can be reused for multiple RHS solves,
and each solve may contain multiple RHS columns.

### Numerical Characteristics

QR uses orthogonal transformations and column pivoting. It reports a numerical
rank estimate and does not exploit Vandermonde structure.

### Failure Semantics

Numerical rank is determined from the diagonal of `R` with the threshold

\[
\tau_{\mathrm{rank}} = 10^{-14}|R_{00}|.
\]

The current relative rank tolerance is the internal default
`relative_rank_tolerance = 1.0e-14`. Rank-deficient factorizations remain
inspectable, but `solve_qr(...)` refuses them. The QR rank tolerance is not
configurable through `.qmom` files.

### Current QMOMKit Role

QR is available to recursive CQMOM only through the explicit `pivoted_qr`
Vandermonde policy. It solves the complete weighted conditional matrix
directly. `automatic` does not select QR, and QR failures do not trigger dense
fallback.

### Limitations And Non-Goals

QR requires `m >= n`. Rank-deficient solves are refused. The implementation does
not produce truncated, pseudoinverse, regularized, or minimum-norm solutions.

## One-Sided Jacobi SVD Solver

### Purpose

`factorize_svd(...)` and `solve_svd(...)` provide a generic square/tall dense
solver with singular-value and numerical-rank diagnostics. In the conditional
Vandermonde adapter, SVD is an explicit direct backend for the complete weighted
system.

### Mathematical Problem

The factorization is

\[
A = U\Sigma V^T,
\]

for `A` with `m >= n`. The implementation stores thin left singular vectors
`U` with dimensions `m x n`, a length-`n` singular-value vector in descending
order, and full right singular vectors `V` with dimensions `n x n`.

### Algorithm

The implementation uses dependency-free one-sided Jacobi rotations to
orthogonalize columns. It applies one global power-of-two scale before the
sweeps to keep the working matrix in a safe magnitude range, then rescales the
singular values afterward. After convergence, singular values are computed from
the orthogonalized column norms and singular triplets are sorted consistently.

### Supported Shapes And Right-Hand Sides

SVD supports square and tall matrices with `m >= n`. Square full-column-rank
systems are ordinary linear solves. Tall full-column-rank systems produce
least-squares solutions. Factorizations can be reused for multiple RHS solves,
and each solve may contain multiple RHS columns.

### Numerical Characteristics

SVD exposes the singular-value spectrum, reports numerical rank, and has an
explicit nonconvergence status for the Jacobi iteration. It performs more
computational work than the simpler direct kernels, but this document does not
make performance or accuracy rankings.

### Failure Semantics

Numerical rank is determined from sorted singular values. A singular value is
counted as active when

\[
\sigma_j > 10^{-14}\sigma_0.
\]

The current relative rank tolerance is the internal default
`relative_rank_tolerance = 1.0e-14`. The default Jacobi convergence tolerance is
resolved as

\[
\tau_{\mathrm{Jacobi}}
=
\max(10^{-14},10\epsilon\max(m,n)),
\]

and the default sweep limit is

\[
N_{\mathrm{sweeps}}
=
\max(32,8n).
\]

Rank-deficient and nonconverged factorizations remain inspectable, but
`solve_svd(...)` refuses them. SVD rank and Jacobi controls are internal
defaults and are not configurable through `.qmom` files.

### Current QMOMKit Role

SVD is available to recursive CQMOM only through the explicit `svd`
Vandermonde policy. It solves the complete weighted conditional matrix
directly. `automatic` does not select SVD, and SVD failures do not trigger
dense fallback.

### Limitations And Non-Goals

SVD requires `m >= n`. Rank-deficient and nonconverged solves are refused. The
implementation does not produce pseudoinverse, truncated, regularized, or
minimum-norm recovery solutions.

## Conditional Vandermonde Adapter

### Weighted Conditional System

Conditional moment reconstruction solves

\[
A C = M,
\]

with

\[
A_{p,i}=w_i x_i^p,
\qquad
A=V^T\operatorname{diag}(w).
\]

Rows of `A` and `M` are parent moment orders, columns of `A` are parent nodes,
rows of `C` are parent nodes, and columns of `C` and `M` are target moment
orders or flattened branch-local moment indices.

### Backend Responsibilities

The dense, QR, and SVD backends build the complete weighted matrix `A` and solve
`A C = M` directly. The BP backend uses the unweighted BP kernel and then lets
the adapter recover conditional moments by dividing by weights.

The adapter owns dimension and finite-input validation, policy resolution,
weighted matrix construction, weight-singularity protection for BP, residual
acceptance, fallback, diagnostics, and performance counters.

### Björck-Pereyra Weight Recovery

For BP, the adapter solves

\[
V^T Y = M
\]

with the unweighted BP kernel. It then recovers conditional moments with

\[
C_{i,q}=Y_{i,q}/w_i.
\]

Before this division, the adapter rejects weights satisfying

\[
|w_i| \le \tau_{\mathrm{sing}}\max(1,\max_j |w_j|).
\]

### Diagnostics And Performance Counters

Adapter results record the requested solver, selected solver after policy
resolution, algorithm used internally, and whether dense fallback occurred.
Recursive CQMOM propagates the public diagnostic fields it exposes and records
branch-local residual and condition diagnostics when detailed diagnostics are
requested. Performance counters distinguish Vandermonde solves, dense solves,
BP solves, QR solves, SVD solves, residual checks, condition estimates, and
dense Vandermonde fallbacks when instrumentation is enabled.

## Solver-Selection Policy

`VandermondeSolver` is the user-facing policy for CQMOM conditional
Vandermonde solves. User-facing syntax and examples belong in
[configuration.md](../user/configuration.md).

| Requested policy | Backend sequence |
|---|---|
| `dense` | Dense solve |
| `pivoted_qr` | Pivoted QR solve |
| `svd` | SVD solve |
| `bjorck_pereyra` | Björck-Pereyra solve |
| `bjorck_pereyra_with_dense_fallback` | Björck-Pereyra solve, then dense retry if needed |
| `automatic` | Björck-Pereyra solve |

`automatic` resolves to `bjorck_pereyra`. QR and SVD are explicit policies.
The conditional Vandermonde adapter owns the dense retry for
`bjorck_pereyra_with_dense_fallback`.

## Failure, Rank, Residual, And Tolerance Semantics

### Solver Failures

Dense failures are singular-pivot failures. BP failures are node-separation,
invalid-input, structured-solve, or non-finite-output failures. QR failures
include invalid inputs and numerical rank deficiency. SVD failures include
invalid inputs, numerical rank deficiency, and Jacobi nonconvergence.

Residual rejection is an adapter-level failure. QMOMKit does not silently
replace these failures with pseudoinverse, truncated, regularized, or
minimum-norm solutions.

### Residual Validation

Explicit BP, QR, and SVD policies use residual validation when requested by
`VandermondeSolveOptions::validate_residual`. BP-with-dense-fallback requires
BP residual acceptance before retaining the BP result. Residual rejection in
that fallback-capable policy may lead to a dense retry.
Explicit QR and SVD residual rejection returns failure without fallback.
Recursive diagnostics do not themselves request residual validation; they
compute diagnostic residuals after successful branch-local solves when detailed
diagnostics are enabled.

The adapter computes the residual magnitude with `residual_norm(...)`, which is
a square-root sum of squared residual entries, and scales it with infinity
norms. The acceptance criterion is

\[
\|AC-M\|_F
\le
\tau_{\mathrm{res}}
\max\left(
1,
\|A\|_\infty\|C\|_\infty,
\|M\|_\infty
\right).
\]

### Tolerance Inventory

| Tolerance | Default | Meaning | Configurable through `.qmom` |
|---|---:|---|---|
| Vandermonde singular tolerance | \(10^{-14}\) | Dense pivots, BP node separation, and BP weight protection through their respective scaled criteria | No documented `.qmom` key; available through direct C++ conditional options |
| Residual tolerance | \(10^{-10}\) | Conditional-system residual acceptance in the adapter | No |
| QR relative-rank tolerance | \(10^{-14}\) | Threshold relative to \(|R_{00}|\) | No |
| SVD relative-rank tolerance | \(10^{-14}\) | Threshold relative to \(\sigma_0\) | No |
| SVD Jacobi tolerance | \(\max(10^{-14},10\epsilon\max(m,n))\) | Column-orthogonality convergence | No |
| SVD maximum sweeps | \(\max(32,8n)\) | Jacobi iteration ceiling | No |

Similarly named tolerances should not be conflated. The Vandermonde singular
tolerance is not the QR or SVD rank tolerance, and SVD Jacobi convergence
controls are independent of conditional-system residual acceptance.

## Verification And Testing

Linear-algebra coverage should remain layered with the implementation:

- direct unit tests cover each numerical kernel and its failure paths;
- conditional-adapter tests use weighted manufactured systems and policy checks;
- recursive CQMOM tests verify policy propagation and status mapping through
  branch-local solves;
- configuration-driven verification checks the public solver policy values
  against analytical CQMOM references and reconstructed moments;
- independent analytical references and cross-solver comparisons provide
  correctness evidence without relying only on agreement with one backend;
- instrumented and non-instrumented builds preserve performance-counter
  behavior without making counters part of solver correctness;
- benchmarks may compile and measure the solvers, but benchmark timings are not
  correctness evidence and are not summarized as documentation claims here.

### Borderline-Conditioning Example

[`examples/cqmom_2d_borderline_conditioning/`](../../examples/cqmom_2d_borderline_conditioning/)
is a worked example for studying conditional Vandermonde
conditioning, explicit solver policies, fallback, residuals, and fixed-count
failure alongside lower-node Wheeler reconstructions. It uses a deterministic
four-atom population and diagnostics reported by CQMOM; its local 2x2 matrix
analysis is example-side diagnostic data, not production QR or SVD metadata.

The example produces two complementary outputs. `output/report.md` is the
human-readable comparison of solver outcomes, failed stages, the dense retry,
and lower-node reconstructions. `output/diagnostics.txt` is the complete flat dataset for
scripts, plotting, residual analysis, condition-estimate inspection, and
detailed numerical comparison. Start with the report for interpretation and use
the flat file when individual rows or diagnostic values are needed.

The four direct solvers are compared to determine whether a failed
reconstruction depends on the linear solver. Dense, Björck–Pereyra, pivoted QR,
and SVD take different numerical routes through the same formulation.
Disagreement can expose solver-specific sensitivity. If all four fail the same
two-node reconstruction, the formulation has lost numerical independence and
changing the solver cannot recover the missing information. The BP-with-dense-
fallback policy is useful for identifying a BP-specific failure that a dense
retry can accept; it does not repair rank loss in the represented system.

The example keeps three metrics conceptually separate:

- the relative singular value measures column independence;
- the 2-norm condition number measures sensitivity; and
- the reconstructed-moment error measures accuracy for the current data.

A residual is a backward-error-style consistency check, not a robustness
metric. The Markdown tables compare the four direct solvers and omit BP+dense
as a separate column; the flat diagnostics file retains every requested policy
and every reported field.

Use the example as a workflow for evaluating another CQMOM formulation:

1. Run the example unchanged and inspect the grouped metrics in
   `output/report.md`.
2. Compare `x_then_y` and `y_then_x` at the same separation.
3. Compare explicit solver and fallback outcomes at the same separation.
4. Inspect failures and fallback rows, then use `output/diagnostics.txt` for
   residual, condition-estimate, and moment-error detail.
5. Change one formulation parameter at a time and repeat the comparison.
6. Enable adaptation separately to determine whether a lower-node quadrature
   remains valid.

The case can show sensitivity to conditioning order, weak parent systems,
solver-policy and fallback outcomes, reported failure location, and successful
lower-node reconstructions. It cannot directly show production QR diagonal
values, production SVD singular values or sweep counts, or internal algorithm
metadata beyond the selected policy reported by the example. It is a focused
workflow, not a generic production diagnostic API, and its threshold values
must not be copied blindly into other formulations.

See the [example README](../../examples/cqmom_2d_borderline_conditioning/README.md)
for the prescribed population, interpretation, and standalone build command.

## Current Limitations And Deferred Features

The following behaviors are not currently implemented:

- automatic selection of QR or SVD;
- QR or SVD fallback chains;
- pseudoinverse, truncated-SVD, regularized, or minimum-norm recovery solves;
- underdetermined systems;
- `.qmom` configuration of QR or SVD rank tolerances;
- public singular-value or numerical-rank diagnostics for CQMOM users;
- adaptive CQMOM node-count changes based on QR or SVD numerical rank.

These are current non-features, not roadmap commitments.
