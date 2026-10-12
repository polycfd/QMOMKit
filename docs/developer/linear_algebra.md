# Linear algebra

QMOMKit has no mandatory numerical backend. The few dense operations it needs
are small internal kernels under [`src/internal/`](../../src/internal/), none
of them public. The matrices are $n\times n$ with $n$ the node count of one
coordinate, so plain $O(n^3)$ algorithms suffice. `DenseMatrix`
([`src/internal/linear_algebra.hpp`](../../src/internal/linear_algebra.hpp)) is
row-major, with bounds asserted in Debug builds.

Every kernel reports failure through a flag or status; a solution returned with
a failure is diagnostic. Every rank, singularity and convergence test is
relative to the matrix it is made on, so scaling the whole matrix does not
change it. The conditional Vandermonde systems are also scaled by powers of
two, which is exact, so that their decisions depend on neither the number
density nor the unit of the coordinate.

## Symmetric eigensolver

`symmetric_jacobi_eigen` in
[`src/internal/linear_algebra.cpp`](../../src/internal/linear_algebra.cpp) is
the classical Jacobi method. It diagonalizes the Jacobi matrix of the Wheeler
inversion ([QMOM](qmom.md)) and the moment matrices of the realizability
checks ([Realizability](realizability.md)).

The matrix is first scaled by the power of two that places its largest entry in
$[1,2)$. Rotation parameters are ratios of entries, so the rotations are those
of the unscaled matrix, while the stopping test becomes relative to its largest
entry. Each step annihilates the largest off-diagonal entry. The iteration
stops when that entry is at most `JacobiOptions::tolerance`, default
$\varepsilon$, or after `max_rotations` rotations, default $10n^2$, a generous
multiple of what tridiagonal and Hankel matrices need. The quadratic
convergence of the last rotations makes the machine-precision tolerance nearly
free.

An exhausted budget sets `converged = false`; Wheeler then fails with
`not_realizable` and the realizability check reports
`numerically_indeterminate`. Eigenvalues are returned in ascending order with
their eigenvectors as columns.

## Dense solvers

| Method | Function, file | Used by | Failure test | Tolerance |
| --- | --- | --- | --- | --- |
| Gaussian elimination, partial pivoting | `solve_dense`, [`linear_algebra.cpp`](../../src/internal/linear_algebra.cpp) | `dense` solver, condition estimate | largest pivot candidate $\le\tau\,\lVert A\rVert_\infty$ of the input matrix | $\tau$ = `singular_tolerance` |
| Björck–Pereyra | `solve_bjorck_pereyra`, [`bjorck_pereyra.cpp`](../../src/internal/bjorck_pereyra.cpp) | `bjorck_pereyra`, `automatic` | $n\,\beta\,\tau\ge1$, $\beta$ the Gautschi bound below | $\tau$ = `singular_tolerance` |
| Householder QR, column pivoting | `factorize_qr`, `solve_qr`, [`householder_qr.cpp`](../../src/internal/householder_qr.cpp) | `pivoted_qr` | numerical rank $<n$: $\lvert R_{kk}\rvert\le10^{-14}\lvert R_{00}\rvert$ | `relative_rank_tolerance` = `1e-14`, [`householder_qr.hpp`](../../src/internal/householder_qr.hpp) |
| One-sided Jacobi SVD | `factorize_svd`, `solve_svd`, [`jacobi_svd.cpp`](../../src/internal/jacobi_svd.cpp) | `svd` | numerical rank $<n$: $\sigma_j\le10^{-14}\sigma_0$; no convergence | `SvdOptions::relative_rank_tolerance` = `1e-14`, [`jacobi_svd.hpp`](../../src/internal/jacobi_svd.hpp) |
| Truncated Moore–Penrose | `solve_svd_pseudoinverse`, [`jacobi_svd.cpp`](../../src/internal/jacobi_svd.cpp) | `moore_penrose` | no convergence; rank deficiency is truncated, not a failure | keeps $\sigma_j>\texttt{rcond}\,\sigma_0$ |

All five also fail on non-finite input or output. QR and SVD accept tall
matrices and then return the least-squares solution, without normal
equations; production systems are square. The QR and SVD rank tolerances are
their own constants and do not follow `singular_tolerance`.

`condition_estimate` returns
$\lVert A\rVert_\infty\lVert A^{-1}\rVert_\infty$, with $A^{-1}$ from
`solve_dense` on the identity, and infinity when that solve fails.
`residual_norm` returns $\lVert AX-B\rVert_F$.

### Björck–Pereyra

`solve_bjorck_pereyra` solves the transposed Vandermonde system $V^TY=B$,
$V_{i,p}=\xi_i^p$, for every column of $B$ in $O(n^2)$ operations per column
without forming $V$ [[Björck & Pereyra, 1970](https://doi.org/10.2307/2004623)].
Its singularity test divides the nodes by the power of two $2^e$ above
$\max_i|\xi_i|$, so $\tilde\xi_i=\xi_i/2^e\in(-1,1)$ and
$\lVert V^T\rVert_\infty\le n$, and bounds the inverse by
[[Gautschi, 1962](https://doi.org/10.1007/BF01386302)]

```math
\lVert V^{-T}\rVert_\infty\le\beta=\max_i\prod_{j\ne i}\frac{1+|\tilde\xi_j|}{|\tilde\xi_i-\tilde\xi_j|},
```

with equality when the nodes share a sign. The system is singular when the
estimated condition number $n\beta$ reaches $1/\tau$. The estimate is
independent of the unit, covers clusters of close nodes as well as a single
close pair, and is infinite for coincident nodes.

### Jacobi SVD

`factorize_svd` scales the matrix by a power of two, then sweeps cyclically
over column pairs. A pair $u$, $v$ counts as orthogonal when
$|u^Tv|\le\tau_J\lVert u\rVert_2\lVert v\rVert_2$, with the products
accumulated in `long double`.
Defaults are $\tau_J=\max(10^{-14},\,10\,\varepsilon\max(m,n))$ and
$\max(32,8n)$ sweeps. Singular values come back in descending order.
`solve_svd` refuses a rank-deficient or unconverged factorization;
`solve_svd_pseudoinverse` accepts a converged rank-deficient one and zeroes
the discarded directions.

## Conditional Vandermonde systems

At a CQMOM branch with parent nodes $\xi_i$ and weights $w_i$, $i<n$, the branch
moments of parent orders $p<n$ determine the conditional moments $C_{i,q}$ of
the remaining coordinates at each parent node ([CQMOM](cqmom.md)):

```math
M_{p,q}=\sum_{i=0}^{n-1}w_i\,\xi_i^p\,C_{i,q},
\qquad\text{i.e.}\qquad
AC=M,\quad A=V^T\operatorname{diag}(w),\quad A_{p,i}=w_i\,\xi_i^p .
```

Each column $q$ is one moment of the remaining coordinates; all columns are
solved together. `solve_conditional_vandermonde_system` in
[`src/internal/vandermonde.cpp`](../../src/internal/vandermonde.cpp) solves the
system; `reconstruct_child_slices` in
[`src/recursive_quadrature.cpp`](../../src/recursive_quadrature.cpp) builds
$M$ and calls it.

### Weights and node scale

No solver sees $A$ as given. Before any solve, a weight with
$|w_i|\le\tau\max_j|w_j|$ fails the system as singular; the test is relative
to the largest weight and applies to every solver. The direct solvers
(`dense`, `pivoted_qr`, `svd`, `moore_penrose`) use the factorization

```math
A=S\,\tilde V^T D,\qquad S=\operatorname{diag}(2^{ep}),\qquad D=\operatorname{diag}(w),
```

with $2^e$ the power of two above $\max_i|\xi_i|$ and $\tilde V$ the Vandermonde
matrix of $\xi_i/2^e$. They solve $\tilde V^TY=S^{-1}M$ and return $C=D^{-1}Y$.
$\tilde V^T$ has a first row of ones and entries in $[-1,1]$, so its pivot, rank
and truncation decisions depend on neither the number density nor the unit.
Björck–Pereyra solves $V^TY=M$ on the unscaled nodes, since its singularity test
scales them itself, and also returns $C=D^{-1}Y$. A non-finite $C$ is a failure.
The Moore–Penrose solution is minimum-norm in $Y$, not in $C$.

### Solver selection

`VandermondeSolver`
([`include/qmomkit/vandermonde_solver.hpp`](../../include/qmomkit/vandermonde_solver.hpp))
is set by `CqmomCase::vandermonde_solver` or `[vandermonde] solver`
([Reference](../user/7_reference.md#cqmom-and-vandermonde)). `automatic`,
the default, always selects `bjorck_pereyra`, the only solver that uses the
Vandermonde structure. Every
explicit choice runs that solver alone. The requested and selected solvers are
recorded in the branch diagnostics.

A solve that fails is reported: the branch gets
`ConditionalMomentStatus::singular_system` and the reconstruction
`CqmomStatus::conditional_moment_failed`. QMOMKit never retries with another
solver, regularizes, or substitutes a pseudoinverse. Truncation happens only
when `moore_penrose` is selected, where a rank-deficient system, for example
two coincident parent nodes, succeeds with the truncated solution.

### Tolerances

| Option | Default | Defined in | Effect |
| --- | --- | --- | --- |
| `singular_tolerance` | `1e-14` | `ConditionalMomentOptions` | Weight test, dense pivot test, Björck–Pereyra bound, condition estimate. |
| `conditioning_warning_threshold` | `1e12` | `ConditionalMomentOptions` | Sets `conditioning_warning` when the condition estimate exceeds it. |
| `moore_penrose_rcond` | `1e-14` | `ConditionalMomentOptions` | Truncation of `moore_penrose`; `[vandermonde] rcond`. |

`ConditionalMomentOptions` is defined in
[`include/qmomkit/conditional_moments.hpp`](../../include/qmomkit/conditional_moments.hpp)
and set through `CqmomCase::conditional_options`; only `rcond` is a `.qmom`
key. `singular_tolerance` and `conditioning_warning_threshold` must be finite
and positive, `rcond` in $[0,1)$; otherwise the case is invalid input.

### Diagnostics

At `CqmomDiagnosticLevel::detailed`,
`CqmomDiagnostics::conditional_moment_diagnostics` describes a successful solve
on the system the solvers decide on, $\tilde V^TY=S^{-1}M$ with $Y=DC$, whatever
solver ran:

- `condition_estimate` is the condition estimate of $\tilde V^T$;
- `residual_norm` is the relative residual
  $\lVert\tilde V^TY-S^{-1}M\rVert_F/\max\bigl(\lVert\tilde V^T\rVert_\infty\lVert Y\rVert_\infty,\lVert S^{-1}M\rVert_\infty\bigr)$;
- `conditioning_warning` is set when the estimate exceeds the threshold.

Neither value depends on the number density, nor on a change of the parent
coordinate's unit by a power of two; another change of unit moves the estimate
by less than a factor $2^{n-1}$. The estimate costs a full inverse, and none of
the three changes the result. No solver checks its residual: a solve fails only
by its own test in the table above. `make_conditional_vandermonde_matrix`,
which forms $A$ explicitly as an independent reference, lives with the tests in
[`tests/conditional_vandermonde_matrix.hpp`](../../tests/conditional_vandermonde_matrix.hpp).

Coverage is in `tests/unit/linear_algebra.cpp`,
`tests/unit/bjorck_pereyra.cpp`, `tests/unit/householder_qr.cpp`,
`tests/unit/jacobi_svd.cpp`, `tests/unit/vandermonde.cpp` and
`tests/verification/vandermonde/`.
