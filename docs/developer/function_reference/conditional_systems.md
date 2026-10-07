# Conditional systems and numerical kernels

## Mathematical and implementation boundary

Recursive conditioning solves `A*C=R`, where `A[p,i]=w[i]*x[i]^p`. Rows are
parent orders; each RHS column is a trailing moment index. `C[i,q]` contains
the child moment at parent node `i`. The [linear algebra
design](../linear_algebra.md) and [recursive design](../cqmom_recursive.md)
explain this correspondence. The private kernels below are `qmomkit::internal`;
file-local helpers have anonymous namespace scope. Public one-level conditional
reconstruction is a separate entry point described below.

## A — conditional Vandermonde adapter

`make_conditional_vandermonde_matrix` and
`solve_conditional_vandermonde_system` are defined in
[`src/internal/vandermonde.cpp`](../../../src/internal/vandermonde.cpp) and
declared in [`vandermonde.hpp`](../../../src/internal/vandermonde.hpp).
Recursive `reconstruct_child_slices` calls the adapter. Matrix construction
loops over parent powers and nodes to form `w*x^p`, with no coordinate
normalization here.

`resolve_solver` retains explicit dense/QR/SVD/Moore-Penrose/BP requests.
`select_automatic_solver` unconditionally returns BP; dimensions and RHS count
are unused. There is **no solver fallback chain**. An unknown internal enum
falls through to dense; the configuration parser rejects unknown textual
values. Returned requested solver, selected solver, and algorithm identify
different stages of that decision.

- **Dense:** check matching dimensions, build `A`, call `solve_dense` with
  `singular_tolerance`. The adapter does not apply residual acceptance.
- **Pivoted QR / SVD:** check dimensions, build `A`, solve using each kernel's
  default rank controls. `singular_tolerance` is not passed as their rank
  tolerance. Rank deficiency maps to failure, not a lower-rank conditional rule.
  Optional `validate_residual` tests the computed result; recursive callers set
  this false regardless of diagnostic level.
- **Moore-Penrose:** check dimensions, build `A`, factor it with the Jacobi
  SVD, then retain precisely `sigma > rcond * sigma_max`. A converged
  rank-deficient factorization is accepted; discarded singular directions are
  zeroed. Invalid `rcond` and SVD nonconvergence fail. Optional residual
  validation applies as for QR and SVD.
- **BP:** reject empty/mismatched/non-finite inputs or nonpositive tolerance.
  Require every `abs(w)>singular_tolerance*max(abs(w))`. Solve the
  unweighted transposed system `V^T*Y=R`, then divide every solution row by its
  corresponding weight. Reject non-finite recovered `C`. Negative finite
  weights are not excluded by this algebraic adapter; quadrature callers supply
  their own accepted weights.

A successful dense, QR, SVD, Moore-Penrose, or BP solve always returns finite
solution values. Residual validation is a separate, optional policy where
supported; rank, pivot, and convergence decisions remain specific to the
selected solver. Success does not by itself establish favorable conditioning or
a requested accuracy.

`residual_is_accepted` requires finite Frobenius residual and compares it to
`residual_tolerance*max(1,||A||inf*||C||inf,||R||inf)`. This is an implemented
mixed-norm acceptance measure, not a forward-error bound. The BP path does not
apply this optional flag. Detailed recursive diagnostics separately compute
residual and condition information without turning warnings into rejection.

The recursive caller supplies these controls as `ConditionalMomentOptions` and
reports the outcome as `ConditionalMomentStatus` with
`ConditionalMomentDiagnostics`, surfaced per branch on `CqmomDiagnostics` and
declared in
[`conditional_moments.hpp`](../../../include/qmomkit/conditional_moments.hpp).
There is no public one-level conditional-moment entry point, so the solver
policy, tolerances and diagnostics described here are the ones that apply to
every conditional solve QMOMKit performs.

Evidence: [`tests/unit/vandermonde.cpp`](../../../tests/unit/vandermonde.cpp)
explicitly checks automatic BP at small/wide and larger sizes, forced modes,
small/zero weights, rank deficiency, zero RHS columns, Moore-Penrose
truncation, and optional QR/SVD/Moore-Penrose residual rejection. Tests
comparing to dense are parity evidence, while
[`solver_selection.cpp`](../../../tests/verification/vandermonde/solver_selection.cpp)
and
[`clustered_solver_policy.cpp`](../../../tests/verification/vandermonde/clustered_solver_policy.cpp)
provide workflow/manufactured comparisons. Cluster separation, weight ratios,
translation, moment powers, and policy thresholds affect the solve; a low
residual alone does not establish accurate child moments.

## A — `solve_bjorck_pereyra`

Definition [`bjorck_pereyra.cpp`](../../../src/internal/bjorck_pereyra.cpp),
declaration [`bjorck_pereyra.hpp`](../../../src/internal/bjorck_pereyra.hpp).
The adapter calls it with nodes and a copied RHS; it knows nothing about
weights. Inputs must be finite, nonempty, dimensionally compatible, and have
positive finite singular tolerance. Zero RHS columns are allowed.

The first outer loop advances the transform index `k`. Its descending row loop
updates `rhs[i+1]-=x[k]*rhs[i]` for each RHS column, preserving the previous
row value needed by the transform. The second outer loop descends `k`; its
ascending row loop divides by `x[i]-x[i-k-1]`, subtracts the quotient from the
previous row, and installs it in the current row. Changing these iteration
directions would change which already transformed values are consumed.

Each denominator must exceed
`singular_tolerance*max(1,abs(lhs_node),abs(rhs_node))` in magnitude. Repeated
or insufficiently separated nodes fail immediately; there is no sorting or
merging. A final finite-matrix check rejects arithmetic overflow. Success
returns `Y` in the original node order. Tests in
[`tests/unit/bjorck_pereyra.cpp`](../../../tests/unit/bjorck_pereyra.cpp)
manufacture `R` from known solutions, cover multiple RHS, mixed-sign/clustered
nodes, repetition, and the transform update count. Small denominator
amplification and subtractive cancellation remain numerical sensitivities.

## A — `solve_dense`, residual and condition calculations

Defined in [`linear_algebra.cpp`](../../../src/internal/linear_algebra.cpp),
declared in [`linear_algebra.hpp`](../../../src/internal/linear_algebra.hpp).
`solve_dense` performs partial row pivoting: each pivot column scans remaining
rows for largest absolute entry, rejects a pivot `<=
tolerance*max(1,||A||inf)`, swaps full matrix and RHS rows, and eliminates the
lower rows across all RHS columns. Reverse row traversal performs triangular
back substitution.

The scale is fixed from the original matrix. The kernel assumes a square
matrix, compatible RHS, finite entries, and a usable tolerance; callers supply
these preconditions. After back-substitution, the kernel rejects a non-finite
final solution and returns failure. Finite inputs do not guarantee finite
arithmetic, so this final gate is required. No rank-reduced or least-squares
result is returned on singularity.

`matrix_infinity_norm` takes maximum absolute row sum; `residual_norm` nests
output row, RHS column, and shared matrix index to compute
`sqrt(sum((A*C-R)^2))`. Both use ordinary double sums. `make_identity` supplies
the RHS for `condition_estimate`, which solves for an explicit inverse and
returns `||A||inf*||A^-1||inf`, or infinity on failed inversion. This is not
the SVD 2-norm condition number. Overflow of norm accumulation and cancellation
in residuals remain numerical limitations. Conditional/Vandermonde tests
exercise these kernels indirectly; no dedicated linear-algebra test file is
present.

## A — pivoted Householder QR

`factorize_qr`, `solve_qr` overloads, and helpers are defined in
[`householder_qr.cpp`](../../../src/internal/householder_qr.cpp), declared in
[`householder_qr.hpp`](../../../src/internal/householder_qr.hpp). The
factorization requires nonempty finite `rows>=cols` and a positive relative
rank tolerance. Each column step recomputes remaining tail norms with repeated
`hypot`, pivots the largest column, and records the column permutation. It
forms a reflector with the sign of the leading entry chosen to avoid the
obvious subtraction cancellation, stores the reflector tail below the diagonal,
and applies it to all trailing columns. A zero norm installs a zero reflector.

`numerical_rank_from_factor` counts the initial diagonal prefix satisfying
`abs(Rii)>relative_rank_tolerance*abs(R00)`; it stops at the first failure.
`solve_qr(factorization,rhs)` validates stored shapes, coefficients,
permutation, status and rank consistency, requires full column rank, applies
reflectors to all RHS columns (`apply_householder_to_rhs`), back-substitutes,
then undoes the column permutation. Tall input therefore supports
full-column-rank least squares; deficient input is refused, not regularized.
The matrix overload factorizes then delegates. The final de-pivoted solution
must be finite before the solve reports success, matching the BP adapter's
finite-output contract.

[`tests/unit/householder_qr.cpp`](../../../tests/unit/householder_qr.cpp)
checks factorization identity, pivot restoration, reusable factors,
exact/inconsistent tall systems, near-dependent rank decisions, scaling, and
malformed factor metadata. These tests justify the listed roles, not arbitrary
dynamic-range safety of reflector construction or back substitution.

## A — one-sided Jacobi SVD

Definitions [`jacobi_svd.cpp`](../../../src/internal/jacobi_svd.cpp),
declarations [`jacobi_svd.hpp`](../../../src/internal/jacobi_svd.hpp).
`factorize_svd` validates finite nonempty `rows>=cols` and options. It scales
the working matrix by a bounded power of two based on the largest entry and
starts right singular vectors at identity. Each sweep visits every column pair.
`rotate_pair` uses column squared norms and cross product (accumulated in `long
double`) to test `abs(gamma)<=jacobi_tolerance*sqrt(alpha*beta)`. Zero norm/dot
also terminates a pair. `jacobi_tangent` computes a stable rotation parameter,
with explicit large/non-finite-zeta branches; both working and right-vector
columns rotate.

A sweep with no rotations establishes convergence. The limit defaults to
`max(32,8*cols)` sweeps; exhausting it reports `no_convergence`, as the
symmetric eigenkernel reports an exhausted rotation budget.
`completed_factorization` uses stable column norms for singular values,
normalizes left vectors, restores scale with `ldexp`, and calls
`sort_singular_triplets` for stable descending ordering. Rank uses a prefix
above `relative_rank_tolerance*sigma_max`. Non-finite completed output is
rejected. The zero matrix returns rank deficient.

`solve_svd` checks factor metadata/status, requires full column rank and
convergence, then loops to compute `U^T*R`, divide by singular values, and
multiply by `V`. `solve_svd_pseudoinverse` uses the same factorization after
convergence, accepts rank-deficient metadata, and applies the strict `rcond`
cutoff before the same `V` multiplication. The final assembled solution must be
finite before either solve reports success. Sensitivities include nearly
parallel columns, cross-product underflow, rank cutoffs, and amplification by
small accepted singular values.

[`tests/unit/jacobi_svd.cpp`](../../../tests/unit/jacobi_svd.cpp) checks known
singular values, orthogonality/reconstruction, repeated values, scale changes,
tall systems, rank deficiency, and forced nonconvergence. These are stronger
than agreement with a previous solver, but do not constitute an envelope.

## Shared symmetric eigenkernel

**A:** `symmetric_jacobi_eigen` and local `rotate` in
`src/internal/linear_algebra.cpp` serve Wheeler and realizability, both with
the default `JacobiOptions`. `power_of_two_shift` first finds the power of two
that places the largest entry in [1, 2), and the matrix is scaled by it. That
scaling is exact, apart from entries so far below the largest that they leave
the normal range, and every rotation parameter is a ratio of entries, so the
rotations are those of the unscaled matrix while the tolerance, machine epsilon
by default, becomes relative to the largest entry. Starting from identity
eigenvectors, each iteration scans the upper triangle for the largest
off-diagonal entry and rotates that pair, updating symmetric matrix entries and
eigenvector columns. Iteration stops when the largest entry is within the
tolerance, which sets `converged`, or when `max_rotations` single rotations
are exhausted; zero selects `10 n^2`, five times the under `2 n^2` measured for
tridiagonal Jacobi and dense Hankel matrices.

Finally sort diagonal values ascending, undo the scaling of the eigenvalues,
and permute eigenvector columns with them. Square, symmetric, finite matrices
are upstream assumptions; a non-finite entry leaves the matrix unscaled.
`tests/unit/linear_algebra.cpp` checks the kernel directly against closed-form
2 x 2 and Gauss-Legendre eigenpairs at scales from `1e-20` to `1e20`, exact
results under power-of-two scaling, polynomial exactness up to forty nodes, and
a reported exhausted budget; Wheeler and realizability tests exercise it
through their callers.

## B/C — supporting helpers

| Functions | Responsibility / consequential semantics | Evidence |
|---|---|---|
| `DenseMatrix` constructor, `rows`, `cols`, `operator()` | Zero-initialized row-major work storage; element indexing uses debug assertions, disabled by `NDEBUG` | `src/internal/linear_algebra.hpp`, numerical kernel callers |
| `matrix_all_finite`, `finite_matrix`, `finite_values` | Reject non-finite entries at the particular caller's validation boundary | Kernel invalid-input tests |
| QR `stable_column_norm`, `swap_columns`, `valid_permutation`, `valid_qr_status`, `rank_metadata_is_consistent`, `factorization_metadata_is_consistent`, `diagnostic_rank_from_metadata` | Stable norms, pivot tracking, and rejecting inconsistent reusable factors; these do not prove the supplied factors represent an original matrix | QR malformed-factor tests |
| SVD `stable_column_norm`, `stable_column_dot`, `max_abs_entry`, `scale_matrix`, `rotate_columns`, `copy_column` | Norm/dot arithmetic, normalization and aligned vector movement for the pair algorithm | SVD identity/scaling tests |
| SVD `resolved_jacobi_tolerance`, `resolved_max_sweeps`, `resolved_scale_exponent` | Resolve zero/default controls and bounded binary scaling | SVD convergence/scaling tests |
| SVD rank, ordering and metadata guards | `numerical_rank_from_singular_values`, `singular_values_are_nonnegative_descending`, `full_rank_singular_values_are_safe`, `valid_svd_status`, `rank_metadata_is_consistent`, `factorization_metadata_is_consistent`, `diagnostic_rank_from_metadata` validate reusable solve inputs | SVD malformed metadata tests |
| `qr_linear_solve_result`, `svd_linear_solve_result`, `vandermonde_result` | Convert kernel status/payload to adapter result; rank failure stays failure | Forced-policy tests |
| `failure_result`, `qr_failure_result`, `factorization_failure_result`, `solve_failure_result`, `invalid_result`, `valid_options` | Construct status-specific empty/zero diagnostic payloads and validate relevant options | Relevant failure tests; no assertion of valid numerical payload on failure |
| BP `node_scale`, `denominator_is_singular` | Relative node-difference threshold with absolute floor one | Repeated/clustered-node tests |
