# Internal Linear Algebra

This page documents the internal dense linear-algebra components used by
QMOMKit. These headers live under `src/internal/` and are implementation
details, not public APIs. Public workflows should continue to use the
documented QMOM, CQMOM, transport, and configuration interfaces.

## Scope

The internal linear-algebra layer provides small dense numerical building
blocks for algorithms that already own their mathematical inputs. It does not
own phase-space hierarchy, configuration parsing, recursive CQMOM policy,
application physics, or public solver selection beyond the internal adapters
described here.

These components are intentionally narrow:

- `linear_algebra.*`: generic dense matrix infrastructure and baseline dense
  solves;
- `bjorck_pereyra.*`: an unweighted transposed Vandermonde kernel;
- `vandermonde.*`: conditional Vandermonde system construction and policy;
- `householder_qr.*`: generic dense pivoted Householder QR;
- `jacobi_svd.*`: generic dense one-sided Jacobi SVD.

## Dense Matrix Infrastructure

`DenseMatrix` is the internal row-major dense matrix container used by the
core numerical kernels. It stores dimensions and finite-sized value arrays, and
is designed for straightforward local algorithms rather than for exposing a
general matrix package.

`solve_dense(...)` solves square systems with one or more right-hand-side
columns using Gaussian elimination with pivoting. It is the simple generic
baseline for dense square systems, residual checks, condition estimates, and
comparisons in tests. It is deliberately independent of Vandermonde structure,
CQMOM recursion, and solver fallback policy.

## Björck-Pereyra

`solve_bjorck_pereyra(...)` solves the unweighted transposed Vandermonde
system

```text
V^T Y = B,
V_{i,p} = x_i^p.
```

Rows of `B` correspond to powers, rows of `Y` correspond to nodes, and each
right-hand-side column is solved independently. This kernel is
structure-specific and owns only the unweighted Björck-Pereyra transformation
and node-separation checks.

The kernel does not know about quadrature weights, conditional moments,
recursive CQMOM traversal, residual validation, dense fallback, or public
Vandermonde solver policy.

## Conditional Vandermonde Adapter

Conditional moment reconstruction solves

```text
A C = M,
A = V^T diag(w),
A[p, i] = w_i x_i^p.
```

The helper `make_conditional_vandermonde_matrix(...)` constructs this complete
weighted conditional matrix. The adapter
`solve_conditional_vandermonde_system(...)` owns the internal policy for this
system:

- the dense backend builds `A` and solves `A C = M` directly;
- the pivoted-QR backend builds `A` and solves `A C = M` directly with the
  internal `solve_qr(...)` entry point;
- the Björck-Pereyra backend solves `V^T Y = M`, then recovers
  `C[i, q] = Y[i, q] / w_i`;
- the adapter owns solver selection, fallback, residual checks, counters, and
  diagnostics.

Weights and conditional reconstruction policy therefore belong to
`vandermonde.*`, not to the Björck-Pereyra kernel. Recursive CQMOM uses this
adapter for branch-local conditional solves; see `cqmom_recursive.md` for the
full recursive construction and index conventions.

## Pivoted Householder QR

`householder_qr.*` provides `factorize_qr(...)` and `solve_qr(...)` for a
generic dense pivoted Householder QR factorization

```text
A P = Q R.
```

It supports square and tall matrices with `m >= n`, column pivoting,
full-column-rank square solves, overdetermined least-squares solves, reusable
factorizations, and multiple right-hand-side columns. The column permutation is
reported as

```text
column_permutation[pivoted_column] = original_column.
```

The factorization stores `R` on and above the diagonal. Entries below each
diagonal store the tail of the Householder vector, whose leading one is
implicit, and one Householder coefficient is stored per reflector.

Numerical rank is determined with the relative threshold

```text
relative_rank_tolerance * abs(R(0, 0)).
```

The current QR default is `relative_rank_tolerance = 1.0e-14`, giving the rank
threshold

\[
\tau_{\mathrm{rank}} = 10^{-14}|R_{00}|.
\]

A `rank_deficient` factorization is inspectable for diagnostics, but solves
are refused unless the factorization is full column rank. QR does not produce
truncated, pseudoinverse, regularized, or minimum-norm solutions for
rank-deficient systems.

Production QR solves do not form explicit `Q`, do not form normal equations,
and do not recompute the factorization during solve-time metadata validation.
Tests may form explicit `Q` only to verify the stored factorization convention.

## One-Sided Jacobi SVD

`jacobi_svd.*` provides `factorize_svd(...)` and `solve_svd(...)` for a
generic dense singular value decomposition

```text
A = U Sigma V^T.
```

It supports square and tall matrices with `m >= n`, thin left singular vectors
`U` with dimensions `m x n`, a length-`n` singular-value vector, and full right
singular vectors `V` with dimensions `n x n`. Singular values are sorted in
descending order, with the matching columns of `U` and `V` permuted
consistently. Columns of `U` associated with exact zero singular values are
stored as zero columns rather than as an orthonormal completion; rank-deficient
solves are refused, so inactive left singular vectors are diagnostic only.

The implementation uses dependency-free one-sided Jacobi sweeps. It maintains
`B = A V`, starts from `B = A` and `V = I`, and applies right rotations to
column pairs of `B` and to the same columns of `V`. A pair is converged when

\[
|b_p^T b_q| \leq \tau_{\mathrm{Jacobi}}\|b_p\|_2\|b_q\|_2.
\]

The default Jacobi tolerance is resolved internally from machine precision and
matrix size, and the default sweep limit scales with the number of columns.
Before sweeps, one global power-of-two scale is applied to keep the working
matrix in a safe magnitude range; singular values are rescaled afterward.
Column-wise scaling is not used because it would change the singular-vector
problem.

Numerical rank is determined from the sorted singular values with

```text
relative_rank_tolerance * sigma[0].
```

The current default is `relative_rank_tolerance = 1.0e-14`. A `rank_deficient`
factorization remains inspectable but is refused by `solve_svd(...)`. If the
Jacobi sweeps do not converge within the resolved limit, the factorization is
returned with `no_convergence` and is also refused by `solve_svd(...)`.

For successful full-column-rank factorizations, `solve_svd(...)` computes

```text
X = V Sigma^{-1} U^T B
```

without forming a pseudoinverse matrix, normal equations, or explicit diagonal
matrices. Square systems are ordinary linear solves, and tall systems produce
full-column-rank least-squares solutions. SVD does not produce pseudoinverse,
truncated, regularized, or minimum-norm solutions for rank-deficient systems.

The SVD component is internal only. Recursive CQMOM can use it only through the
explicit `VandermondeSolver::svd` policy, which solves the complete weighted
conditional Vandermonde system directly through
`solve_conditional_vandermonde_system(...)`. The SVD path uses the internal SVD
defaults: `relative_rank_tolerance = 1.0e-14`, the resolved Jacobi convergence
tolerance, and the resolved maximum sweep count. These SVD tolerances are not
configurable through `.qmom` files and are not taken from
`VandermondeSolveOptions::singular_tolerance`.

`automatic` does not select SVD. SVD is not a BP-to-SVD fallback path and does
not fall back to dense, QR, or BP. Rank deficiency, Jacobi nonconvergence, or
residual-validation failure is reported as a failed conditional Vandermonde
solve rather than as a partial solution, pseudoinverse, truncated solve,
regularized solve, or dense fallback.

## Current Production Status

Pivoted Householder QR is available to recursive CQMOM only through the
explicit `VandermondeSolver::pivoted_qr` policy. It solves the complete
weighted conditional Vandermonde system directly through
`solve_conditional_vandermonde_system(...)` and uses the QR solver's internal
default relative rank tolerance. The QR relative rank tolerance is currently an
internal default, is not configurable through `.qmom`, and is not taken from
`VandermondeSolveOptions::singular_tolerance`. Exposing a user-configurable QR
rank tolerance is deferred until a demonstrated need exists.

`automatic` does not select QR. QR is not a BP-to-QR fallback path, and a QR
rank-deficiency or residual-validation failure is reported as a failed
conditional Vandermonde solve rather than as a partial solution or dense
fallback.
