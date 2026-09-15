# Realizability and runtime acceptance

All non-inline definitions below are in
[`src/realizability.cpp`](../../../src/realizability.cpp), public declarations
in [`realizability.hpp`](../../../include/qmomkit/realizability.hpp). Numerical
helper functions are in a local anonymous namespace. This module classifies a
supplied raw sequence; it neither reconstructs quadrature nor chooses active
node counts. The [existing theory/design page](../realizability.md) describes
the intended Hamburger, Stieltjes, and Hausdorff correspondence. The
distinctions and concern below qualify what the implementation actually checks.

## A — `check_realizability` and `validate_common_inputs`

The public overloads accept an inclusive sequence `m0..mN` and either a support
enum or explicit interval. They are called directly by applications and by the
runtime policy adapters. Common validation rejects empty/non-finite sequences
and nonpositive/non-finite tolerance. Negative mass is `not_realizable`;
exactly zero mass accepts only an exactly all-zero sequence, classified
`boundary`. There is no near-zero-mass cutoff in this classifier.

The support overload returns any non-strict common result **before**
dispatching support. Thus zero-measure or negative-mass early results can
precede rejection of a bounded-support enum without interval. The interval
overload checks finite `lower < upper` before returning the common
classification. This validation precedence is visible in source; the
zero-measure/missing-interval combination has no targeted test identified in
the relevant realizability tests.

Positive-mass valid inputs dispatch to `check_hamburger`, `check_stieltjes`, or
`check_hausdorff`. Their status is a classification, not an exception-based
numerical failure interface. Allocation/index exceptions are not generally
caught by the public wrapper. The all-finite requirement is for the whole
supplied sequence, unlike Wheeler's required-prefix check.

Evidence:
[`tests/unit/realizability.cpp`](../../../tests/unit/realizability.cpp),
`test_invalid_inputs`, `test_m0_policy`, `test_status_helpers`; support-family
verification in
[`realizability_analytical.cpp`](../../../tests/verification/realizability/realizability_analytical.cpp).

## A — matrix builders and support dispatch

The matrix entries explain the mathematical object being tested:

| Function | Loop/index meaning and responsibility |
|---|---|
| `hankel_matrix` | Row/column loops assign `m[row+col+shift]` to an `(order+1)` square matrix; shift zero tests squared polynomials, shift one adds an `x` factor |
| `interval_quadratic_localizing_matrix` | Assigns `(a+b)m[k+1]-ab*m[k]-m[k+2]`, `k=row+col`, for the polynomial `(x-a)(b-x)` |
| `interval_linear_localizing_matrix` | Selects `m[k+1]-a*m[k]` or `b*m[k]-m[k+1]` for the lower/upper endpoint polynomial |
| `extension_column` | Row loop gathers `m[row+matrix_order+1+shift]`, the next column requiring range consistency |

The builders assume the dispatcher's parity/order calculations make all
required indices available; `MomentSet::operator[]` would throw on an invalid
index. They do not scale coordinates or guard every arithmetic operation;
non-finite constructed matrices are detected by `analyze_matrix`.

`check_hamburger` analyzes `H_floor(N/2)`. Odd `N` checks its next extension
column against the matrix range. For a numerically singular even `H_n`, it also
analyzes `H_(n-1)` and requires equal numerical ranks: the singular truncated
matrix must be a flat extension. This rejects PSD-but-incompatible data such as
`[1,0,0,0,1]` without changing the positive-definite or odd-order paths.

`check_stieltjes` treats order zero separately; otherwise it analyzes `H_k`.
Odd `N=2k+1` also analyzes shifted `S_k` and checks the ordinary Hankel next
column. Even `N=2k` analyzes `S_(k-1)` and checks its shifted extension. Matrix
failure stops the sequence of checks immediately.

`check_hausdorff` treats order zero with `H0`. Odd truncation analyzes lower
and upper linear localizers. Even truncation analyzes `H_k` and, for `k>0`, the
quadratic interval localizer. These functions establish the tested matrix
conditions; their names should not replace independent mathematical validation.

`verify_hamburger`, `verify_stieltjes`, and `verify_hausdorff` in the
analytical verification file use Gaussian/discrete, gamma/log-normal/discrete,
and beta/uniform/endpoint cases, respectively. Singular next-column failures
and atomic boundary measures and the incompatible even-order sequence remain
covered. Conditioned matrices, cancellation in localizers, interval
translation/width, and truncation parity remain numerical sensitivities.

## A — `analyze_matrix` and local `symmetric_jacobi_eigen`

`analyze_matrix` first scans all entries for finiteness and the maximum
absolute entry. A second nested loop computes the infinity norm relative to
that maximum, avoiding an overflowing raw row sum in the relative-tolerance
calculation. The effective user threshold is either absolute `tolerance` or
`tolerance*entry_scale*normalized_infinity_norm`.

The local eigen wrapper divides a temporary matrix by its entry scale before
calling the shared symmetric Jacobi kernel (`1e-14`, at most `80*max(n,1)`
rotations), then restores eigenvalue magnitudes. Positive scalar scaling leaves
the intended eigenvectors unchanged. There is no convergence status returned by
that kernel. The classification branches, in order, are:

1. non-finite eigenvalues → indeterminate;
2. minimum eigenvalue below minus the user threshold → not PSD;
3. absolute minimum eigenvalue within
   `1000*epsilon*entry_scale*normalized_infinity_norm` → PSD boundary;
4. minimum eigenvalue below the user threshold → indeterminate;
5. otherwise → positive definite.

Step four includes **small positive** eigenvalues. Its diagnostic therefore
says that an eigenvalue is within the numerical tolerance, without assigning a
sign. The roundoff allowance and user threshold are separate quantities, not
interchangeable acceptance tolerances. Matrix scaling, eigen gaps, and
relative/absolute mode affect this classification.
`verify_near_boundary_and_scaling` checks both signs near zero, low-density
non-diagonal matrices, and an overflowing raw row sum; it does not characterize
arbitrary rank or dimension.

## A — `check_range_consistency`, `incorporate_range_consistency`

For a matrix not classified positive definite, the extension must have small
projection onto every eigenvector whose eigenvalue is `<= analysis.tolerance`.
`projection_on_eigenvector` sums a column dot product; the outer loop retains
the maximum absolute projection. The threshold is absolute tolerance or the
maximum of matrix tolerance and `options.tolerance*vector_scale(extension)`.

Projection above that threshold rejects the sequence; above `0.1*threshold`
reports indeterminate; otherwise the extension is accepted as boundary. The
factor `0.1` creates an uncertainty band. Its purpose is apparent from the
branch, but a quantitative justification for its exact width was not
identified. This is not a least-squares correction: no input moment is changed.
The dispatcher checks matrix status first, so unsuccessful matrix analyses do
not reach the normal range path.

Evidence: singular ordinary/shifted extension examples in the realizability
analytical verification. The `0.1` projection band is a fixed implementation
policy, not a quantitative accuracy bound.

## A/B — policy adapters

Public inline `accepts_realizability_policy`,
`make_realizability_acceptance_predicate`, and
`make_realizability_stage_acceptance_predicate` live in
[`realizability_policy.hpp`](../../../include/qmomkit/realizability_policy.hpp).
A disabled policy accepts without calling the converter. Enabled evaluation
requires a converter, converts once, and accepts strict classification plus
boundary only when configured. Invalid/indeterminate/not-realizable results
reject. Converter exceptions propagate.

The final predicate factory returns a callback only for `FinalOnly`. The stage
factory returns none for disabled policies; for `FinalOnly` it ignores
intermediate stages; for `EveryStage` it checks every supplied stage. The
source and evolution wrappers compose these with caller predicates by logical
AND. CQMOM coordinate policies operate on pure marginals, not all mixed or
conditional moments; see [sources](moment_sources.md).
[`tests/unit/realizability_policy.cpp`](../../../tests/unit/realizability_policy.cpp)
and the coordinate-policy tests in
[`tests/unit/moment_source.cpp`](../../../tests/unit/moment_source.cpp) check
mode selection, status acceptance, and composition.

## B/C — classification bookkeeping

| Function(s) | Purpose and important semantics | Evidence |
|---|---|---|
| `update_limiting_matrix` | Retains the smallest normalized absolute limiting eigenvalue; `minimum_pivot` stores an eigenvalue here, not a Cholesky pivot | Matrix diagnostics tests |
| `incorporate_matrix` | Converts matrix classification to early failure or sets a persistent boundary flag | Support-family verification |
| `finish_classification` | Chooses boundary if any accepted matrix/range check was singular; otherwise strict | Same verification |
| `make_result`, `with_support_context`, `support_name` | Populate diagnostics and append support text when absent | Diagnostic verification helpers |
| `vector_scale`, `projection_on_eigenvector` | Maximum absolute extension entry and ordinary dot product used by range checks | Indirect singular-range coverage |
| `RealizabilityResult::realizable`, `on_boundary` | Exact status predicates; no fresh numerical check | `test_status_helpers` |
