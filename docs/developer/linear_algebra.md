# Internal linear algebra

This document defines the architecture, mathematical conventions, numerical
behavior, and current QMOMKit roles of the internal dense linear-algebra
kernels and the conditional Vandermonde adapter.

The functionalities described here live under `src/internal/` and are
implementation details, not public APIs. This page is for QMOMKit developers
who need to maintain the numerical kernels or their CQMOM integration.
User-facing configuration syntax belongs in
[configuration.md](../user/configuration.md); CQMOM workflow and branch
traversal details belong in [cqmom_recursive.md](cqmom_recursive.md).

## Contents

- [Scope And Design Principles](#scope-and-design-principles)
- [Internal Component Map](#internal-component-map)
- [Common Matrix And Solve Conventions](#common-matrix-and-solve-conventions)
- [Solver Comparison](#solver-comparison)
- [Björck-Pereyra Solver](#björck-pereyra-solver)
- [Dense Gaussian-Elimination Solver](#dense-gaussian-elimination-solver)
- [Pivoted Householder QR Solver](#pivoted-householder-qr-solver)
- [One-Sided Jacobi SVD Solver](#one-sided-jacobi-svd-solver)
- [Vandermonde Adapter](#vandermonde-adapter)
- [Solver-Selection Policy](#solver-selection-policy)
- [Failure, Rank, Residual, And Tolerance Semantics](#failure-rank-residual-and-tolerance-semantics)
- [Verification And Testing](#verification-and-testing)

## Scope and design principles

The internal linear-algebra layer provides dependency-minimal dense numerical
building blocks for small systems whose mathematical inputs are owned by higher
level QMOMKit algorithms. It prioritizes mathematical correctness, explicit
matrix conventions, clear separation between numerical kernels and adapter
policy, and algorithms proportional to the small dense systems used by CQMOM.

Solver failures remain visible. QMOMKit does not silently change the
mathematical problem by substituting pseudoinverse, truncated, regularized, or
minimum-norm recovery paths when a configured solve is not accepted. The
explicit `moore_penrose` policy is the sole requested truncated pseudoinverse
path.

This layer does not own phase-space hierarchy, configuration parsing, recursive
CQMOM traversal, application physics, or public workflow design. It also does
not provide a general-purpose matrix package.

## Internal component map

### Numerical kernels

The numerical kernels own their algorithms, supported matrix shapes,
singularity or numerical-rank detection, and direct solve semantics:

```text
solve_dense(...)
solve_bjorck_pereyra(...)
factorize_qr(...) / solve_qr(...)
factorize_svd(...) / solve_svd(...) / solve_svd_pseudoinverse(...)
```

`solve_dense(...)`, `factorize_qr(...)`, `solve_qr(...)`, `factorize_svd(...)`,
and `solve_svd(...)` are generic dense kernels. `solve_svd_pseudoinverse(...)`
forms an explicitly requested truncated SVD solution.
`solve_bjorck_pereyra(...)` is a structured, unweighted transposed Vandermonde
kernel.

### Conditional Vandermonde adapter

The conditional Vandermonde adapter owns the weighted production system used by
CQMOM conditional reconstruction:

```text
make_conditional_vandermonde_matrix(...)
solve_conditional_vandermonde_system(...)
```

The adapter checks the inputs, resolves the selected `VandermondeSolver`
policy, separates the weights and the node scale from the system, dispatches to
the numerical kernels, recovers the conditional moments by dividing by the
weights, validates QR, SVD, and Moore-Penrose residuals when requested, records
performance counters, and returns solver diagnostics.

### CQMOM integration

CQMOM constructs branch-local conditional systems, propagates the configured
`VandermondeSolver` policy to the adapter, maps failed adapter solves to CQMOM
status values, and records branch diagnostics. It does not implement the dense,
BP, QR, or SVD algorithms; explicitly requested Moore-Penrose recovery remains
adapter-owned. See [cqmom_recursive.md](cqmom_recursive.md) for the CQMOM
algorithm.

## Common matrix and solve conventions

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

## Solver comparison

Five solver policies are available: a specialized Björck-Pereyra solver
[[Björck & Pereyra, 1970](https://doi.org/10.2307/2004623)], a Moore-Penrose
solver that relies on a pseudoinverse [[Marco & Martinez,
2019](https://www.doi.org/10.1016/j.cam.2018.10.009)], as well as generic
solvers based on Gauss eliminiation, QR-factorization, and singular-value
decomposition [[Press,
1992](https://books.google.ca/books/about/Numerical_Recipes_in_FORTRAN_77_Volume_1.html?id=gn_4mpdN9WkC&redir_esc=y)].

| Solver | Mathematical system | Type | Supported shapes | Multiple RHS | Rank diagnostics | Current Vandermonde policy |
|---|---|---|---|---|---|---|
| Björck-Pereyra | $V^T Y=B$ | Structured, unweighted Vandermonde | Square Vandermonde | Yes | No; reports singularity from a condition estimate | `bjorck_pereyra`; backend for `automatic` |
| Dense Gaussian elimination | $AX=B$ | Generic dense | Square | Yes | No; reports singular-pivot failure | `dense` |
| Pivoted Householder QR | $AP=QR$, solve $AX\approx B$ | Generic dense | Square and tall, $m\ge n$ | Yes | Numerical rank | `pivoted_qr` |
| One-sided Jacobi SVD | $A=U\Sigma V^T$, solve $AX\approx B$ | Generic dense | Square and tall, $m\ge n$ | Yes | Singular values and numerical rank | `svd` |
| Truncated Moore-Penrose SVD | $A^+B$ | Generic dense | Square and tall, $m\ge n$ | Yes | Retained singular directions | `moore_penrose` |


## Björck-Pereyra solver

### Purpose

`solve_bjorck_pereyra(...)` exploits Vandermonde structure for the conditional
Vandermonde pathway.

### Mathematical problem

For nodes `x_i`, define

$$
V_{i,p}=x_i^p.
$$

The BP kernel solves the unweighted transposed Vandermonde system

$$
V^T Y = B.
$$

Rows of `B` correspond to powers, rows of `Y` correspond to nodes, and columns
of `B` are independent right-hand sides. The kernel is unweighted; quadrature
weights belong to the conditional Vandermonde adapter.

### Algorithm

The implementation applies structured $\Theta(n^2 r)$ Vandermonde transformations
equivalent to divided-difference and Newton-basis operations. It avoids forming
or factoring a generic dense matrix and applies the same structured operations
to every right-hand-side column, where $n$ is the Vandermonde size and $r$
is the right-hand-side width. For fixed $r$, the kernel work is
$\Theta(n^2)$. Optional internal instrumentation counts the executed
transform updates exactly as $n(n-1)r$; this is a deterministic work-structure
check, not a wall-clock performance guarantee.

### Supported shapes and right-hand sides

The BP kernel supports square Vandermonde systems with one node per unknown and
one RHS row per power. It supports multiple right-hand-side columns. It is not
a generic rectangular solver.

### Numerical characteristics

BP exploits exact Vandermonde structure and avoids the overhead of a generic
dense factorization. Its behavior depends on the structured problem, the node
values, and node ordering. Its accuracy degrades with the conditioning of the
node set; it estimates the condition number for its singularity decision but
provides no numerical-rank estimate.

### Failure semantics

The kernel rejects empty systems, inconsistent dimensions, invalid tolerances,
non-finite inputs, non-finite outputs, and node sets whose system is too
ill-conditioned to solve. With the nodes divided by the power of two above the
largest $|x_i|$, written $\xi_i$, every entry of $V^T$ lies in $[-1,1]$, so
$\|V^T\|_\infty\le n$, and Gautschi's bound
[[Gautschi, 1962](https://doi.org/10.1007/BF01386302), Theorem 1]

$$
\|V^{-T}\|_\infty
\le
\max_i\prod_{j\ne i}\frac{1+|\xi_j|}{|\xi_i-\xi_j|}
$$

holds with equality when the nodes share a sign. The system is treated as
singular when $n$ times this bound reaches $1/\tau_{\mathrm{sing}}$, the
threshold of the dense pivot test and of the QR and SVD rank tests. The
estimate does not depend on the unit of the coordinate, covers several nodes
that are close together as well as a single close pair, and is infinite for
coincident nodes. The kernel reports this as a transposed Vandermonde singular
failure. It does not perform weighted residual validation.

### Current QMOMKit role

BP is the **default** solver, used by the explicit `bjorck_pereyra` policy and
by `automatic`.

### Limitations and non-goals

The BP kernel is Vandermonde-specific, unweighted, and square-only. It does not
construct conditional systems, divide by quadrature weights, or provide rank
diagnostics.


## Dense Gaussian-elimination solver

### Purpose

`solve_dense(...)` is the simple generic square dense solver. It provides a
baseline for dense square systems and the explicit `dense` Vandermonde policy.

### Mathematical problem

The dense solver computes `X` in

$$
A X = B,
\qquad
A \in \mathbb{R}^{n \times n}.
$$

`B` may contain one or more right-hand-side columns.

### Algorithm

The implementation uses Gaussian elimination with row pivoting, followed by
triangular back-substitution. All right-hand-side columns are advanced through
the same elimination, so multiple RHS columns are solved simultaneously.

### Supported shapes and right-hand sides

The dense solver supports square matrices only. It accepts multiple
right-hand-side columns. It does not support rectangular least-squares
problems.

### Numerical characteristics

The solver is general-purpose and does not exploit Vandermonde structure. It
does not produce a numerical-rank decomposition or singular values. Its main
role is a straightforward independent baseline and a robustness reference for
the conditional Vandermonde adapter.

### Failure semantics

At elimination step `k`, the pivot is the largest absolute entry in column `k`
over rows `k..n-1`. The solve fails when the selected pivot satisfies

$$
\max_{i \ge k} |a_{ik}|
\le
\tau_{\mathrm{sing}}\|A\|_\infty,
$$

where the infinity norm is computed from the original matrix passed into the
solve, so scaling the whole matrix does not change the decision. Failure
returns the caller-provided singular message and no accepted solution.

### Current QMOMKit role

The dense solver is used by the explicit `dense` Vandermonde policy and by the
condition estimate of detailed CQMOM diagnostics. The `dense` policy solves the
scaled system described in [Weight and scale separation](#weight-and-scale-separation),
not the weighted matrix as given. The condition estimate inverts the weighted
matrix itself.

### Limitations and non-goals

The dense solver is square-only. It does not report numerical rank, compute
least-squares solutions, regularize singular systems, or expose a public matrix
API.

## Pivoted Householder QR solver

### Purpose

`factorize_qr(...)` and `solve_qr(...)` provide a generic square/tall dense
solver with numerical-rank detection. In the conditional Vandermonde adapter,
QR is an explicit direct backend for the scaled conditional system.

### Mathematical problem

The factorization is

$$
A P = Q R,
$$

for `A` with `m >= n`. The stored permutation satisfies

```text
column_permutation[pivoted_column] = original_column.
```

### Algorithm

The implementation uses Householder reflectors with column pivoting. Solves
apply $Q^T$ implicitly to the right-hand side and then use triangular
back-substitution in the pivoted coordinates before unpermuting the solution.
The factorization stores reflector data compactly in the factor matrix and a
coefficient vector.

### Supported shapes and right-hand sides

QR supports square and tall matrices with `m >= n`. Square full-column-rank
systems are ordinary linear solves. Tall full-column-rank systems produce
least-squares solutions. Factorizations can be reused for multiple RHS solves,
and each solve may contain multiple RHS columns.

### Numerical characteristics

QR uses orthogonal transformations and column pivoting. It reports a numerical
rank estimate and does not exploit Vandermonde structure.

### Failure semantics

Numerical rank is determined from the diagonal of `R` with the threshold

$$
\tau_{\mathrm{rank}} = 10^{-14}|R_{00}|.
$$

The current relative rank tolerance is the internal default
`relative_rank_tolerance = 1.0e-14`. Rank-deficient factorizations remain
inspectable, but `solve_qr(...)` refuses them. The QR rank tolerance is not
configurable through `.qmom` files.

### Current QMOMKit role

QR is available to CQMOM through the explicit `pivoted_qr` policy. It solves
the scaled system described in [Weight and scale separation](#weight-and-scale-separation),
so its rank decision depends on neither the number density nor the unit of the
coordinate.

### Limitations and non-goals

QR requires `m >= n`. Rank-deficient solves are refused. The implementation
does not produce truncated, pseudoinverse, regularized, or minimum-norm
solutions.

## One-sided Jacobi SVD solver

### Purpose

`factorize_svd(...)` and `solve_svd(...)` provide a generic square/tall dense
solver with singular-value and numerical-rank diagnostics. In the conditional
Vandermonde adapter, SVD is an explicit direct backend for the scaled
conditional system.

### Mathematical problem

The factorization is

$$
A = U\Sigma V^T,
$$

for `A` with `m >= n`. The implementation stores thin left singular vectors
`U` with dimensions `m x n`, a length-`n` singular-value vector in descending
order, and full right singular vectors `V` with dimensions `n x n`.

### Algorithm

The implementation uses dependency-free one-sided Jacobi rotations to
orthogonalize columns. It applies one global power-of-two scale before the
sweeps to keep the working matrix in a safe magnitude range, then rescales the
singular values afterward. After convergence, singular values are computed from
the orthogonalized column norms and singular triplets are sorted consistently.

### Supported shapes and right-hand sides

SVD supports square and tall matrices with $m \geq n$. Square full-column-rank
systems are ordinary linear solves. Tall full-column-rank systems produce
least-squares solutions. Factorizations can be reused for multiple RHS solves,
and each solve may contain multiple RHS columns.

### Numerical characteristics

SVD exposes the singular-value spectrum, reports numerical rank, and has an
explicit nonconvergence status for the Jacobi iteration. It performs more
computational work than the simpler direct kernels.

### Failure semantics

Numerical rank is determined from sorted singular values. A singular value is
counted as active when

$$
\sigma_j > 10^{-14}\sigma_0.
$$

The current relative rank tolerance is the internal default
`relative_rank_tolerance = 1.0e-14`. The default Jacobi convergence tolerance
is resolved as

$$
\tau_{\mathrm{Jacobi}}
=
\max(10^{-14},10\epsilon\max(m,n)),
$$

and the default sweep limit is

$$
N_{\mathrm{sweeps}}
=
\max(32,8n).
$$

Rank-deficient and nonconverged factorizations remain inspectable, but
`solve_svd(...)` refuses them. `solve_svd_pseudoinverse(...)` accepts a
converged full-rank or rank-deficient factorization and retains exactly
$\sigma_j > rcond\,\sigma_0$. It zeros discarded directions, returns a finite
zero solution when none are retained, and reports the retained rank. It rejects
nonconverged factorizations and invalid `rcond` values. SVD rank and Jacobi
controls are internal defaults and are not configurable through `.qmom` files.

### Current QMOMKit role

SVD is available to CQMOM through the explicit `svd` policy. It solves the
scaled system described in [Weight and scale separation](#weight-and-scale-separation),
so its rank decision depends on neither the number density nor the unit of the
coordinate. The separately explicit `moore_penrose` policy reuses the Jacobi
factorization of the same scaled system and applies the configured
`vandermonde.rcond` cutoff to its singular values. Its truncated solution is
therefore the minimum-norm solution in the scaled unknowns `Y` described below,
not in the conditional moments `C`.

### Limitations and non-goals

SVD requires $m \geq n$. Regular SVD solves and nonconverged factorizations are
refused. The implementation does not provide regularization or an implicit
recovery path; truncated minimum-norm recovery is available only through the
explicit Moore-Penrose policy.

## Vandermonde adapter

### Weighted conditional system

Conditional moment reconstruction solves

$$
A C = M,
$$

with

$$
A_{p,i}=w_i x_i^p,
\qquad
A=V^T\operatorname{diag}(w).
$$

Rows of `A` and `M` are parent moment orders, columns of `A` are parent nodes,
rows of `C` are parent nodes, and columns of `C` and `M` are target moment
orders or flattened branch-local moment indices.

### Weight and scale separation

No backend solves `A C = M` with the weighted matrix as given. For BP, the
adapter solves

$$
V^T Y = M
$$

with the unweighted BP kernel, whose singularity test does not depend on the
node magnitudes. For the dense, QR, SVD, and Moore-Penrose backends, it factors

$$
A = S\,\tilde V^T D,
\qquad
S=\operatorname{diag}(s^p),
\qquad
D=\operatorname{diag}(w),
$$

where $\tilde V$ is the unweighted Vandermonde matrix of the scaled nodes
$x_i/s$ and $s$ is the power of two above $\max_i |x_i|$. The backend solves

$$
\tilde V^T Y = S^{-1} M.
$$

Every entry of $\tilde V^T$ lies in $[-1,1]$ and its first row is all ones, and
both scalings are exact powers of two. The singularity, rank, truncation, and
residual decisions of these backends therefore depend on neither the number
density nor the unit of the coordinate.

Every policy then recovers the conditional moments with

$$
C_{i,q}=Y_{i,q}/w_i,
$$

and fails if the division overflows. Before solving, the adapter rejects
weights satisfying

$$
|w_i| \le \tau_{\mathrm{sing}}\max_j |w_j|.
$$

The test is relative to the largest weight, so it does not depend on the
number density that the weights carry. It applies to every policy, including
Moore-Penrose.

The adapter owns dimension and finite-input validation, policy resolution,
weight and scale separation, weight-singularity protection, residual
acceptance, diagnostics, and performance counters.

### Diagnostics and performance counters

Adapter results record the requested solver, selected solver after policy
resolution, and algorithm used internally.
Recursive CQMOM propagates the public diagnostic fields it exposes and records
branch-local residual and condition diagnostics when detailed diagnostics are
requested. Performance counters distinguish Vandermonde solves, dense solves,
BP solves, QR solves, SVD solves, residual checks, and condition estimates when
instrumentation is enabled.

## Solver-selection policy

`VandermondeSolver` is the user-facing policy for CQMOM conditional
Vandermonde solves. User-facing syntax and examples belong in
[configuration.md](../user/configuration.md).

| Requested policy | Backend sequence |
|---|---|
| `dense` | Dense solve |
| `pivoted_qr` | Pivoted QR solve |
| `svd` | SVD solve |
| `moore_penrose` | Truncated SVD pseudoinverse solve |
| `bjorck_pereyra` | Björck-Pereyra solve |
| `automatic` | Björck-Pereyra solve |

`automatic` resolves to `bjorck_pereyra`.

## Failure, rank, residual, and tolerance semantics

### Solver failures

Dense failures are singular-pivot failures. BP failures are condition-estimate,
invalid-input, structured-solve, or non-finite-output failures. QR failures
include invalid inputs and numerical rank deficiency. SVD failures include
invalid inputs, numerical rank deficiency, and Jacobi nonconvergence. The
explicit Moore-Penrose policy accepts converged rank-deficient SVDs but still
fails for invalid input, invalid `rcond`, and Jacobi nonconvergence.

Residual rejection is an adapter-level failure. QMOMKit does not silently
replace these failures with pseudoinverse, truncated, regularized, or
minimum-norm solutions; the Moore-Penrose path is used only when explicitly
requested.

### Residual validation

Explicit QR, SVD, and Moore-Penrose policies use residual validation when
requested by `VandermondeSolveOptions::validate_residual`. Residual rejection
returns failure. Recursive diagnostics do not themselves request residual
validation; they compute diagnostic residuals after successful branch-local
solves when detailed diagnostics are enabled.

The adapter computes the residual magnitude of the scaled system with
`residual_norm(...)`, which is a square-root sum of squared residual entries,
and scales it with infinity norms. The acceptance criterion is

$$
\|\tilde V^T Y-S^{-1}M\|_F
\le
\tau_{\mathrm{res}}
\max\left(
\|\tilde V^T\|_\infty\|Y\|_\infty,
\|S^{-1}M\|_\infty
\right).
$$

The bound has no absolute floor, so the decision does not depend on the
magnitude of the moments.

### Tolerance inventory

| Tolerance | Default | Meaning | Configurable through `.qmom` |
|---|---:|---|---|
| Vandermonde singular tolerance | $10^{-14}$ | Dense pivots, the BP condition estimate, and weight protection for every policy through their respective scaled criteria | No documented `.qmom` key; available through direct C++ conditional options |
| Residual tolerance | $10^{-10}$ | Conditional-system residual acceptance in the adapter | No |
| QR relative-rank tolerance | $10^{-14}$ | Threshold relative to $|R_{00}|$ | No |
| SVD relative-rank tolerance | $10^{-14}$ | Threshold relative to $\sigma_0$ | No |
| SVD Jacobi tolerance | $\max(10^{-14},10\epsilon\max(m,n))$ | Column-orthogonality convergence | No |
| SVD maximum sweeps | $\max(32,8n)$ | Jacobi iteration ceiling | No |
| Moore-Penrose `rcond` | $10^{-14}$ | Strict relative singular-value cutoff for the explicit pseudoinverse | Yes, `vandermonde.rcond` |

Similarly named tolerances should not be conflated. The Vandermonde singular
tolerance is not the QR or SVD rank tolerance, and SVD Jacobi convergence
controls are independent of conditional-system residual acceptance.

## Verification and testing

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
