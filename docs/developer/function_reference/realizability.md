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

Positive-mass valid inputs are first passed through `tested_moments` (see
below) and then dispatch to `check_hamburger`, `check_stieltjes`, or
`check_hausdorff`. Their status is a classification, not an exception-based
numerical failure interface. Allocation/index exceptions are not generally
caught by the public wrapper. The all-finite requirement is for the whole
supplied sequence, unlike Wheeler's required-prefix check.

Evidence:
[`tests/unit/realizability.cpp`](../../../tests/unit/realizability.cpp),
`test_invalid_inputs`, `test_m0_policy`, `test_status_helpers`; support-family
verification in
[`realizability_analytical.cpp`](../../../tests/verification/realizability/realizability_analytical.cpp).

## A — standardization: `tested_moments`, `standardized_moments`, `raw_moments`

`tested_moments` returns `standardized_moments` in relative-tolerance mode and
`raw_moments` otherwise. Both produce a `TestedMoments` value: a
`MomentSequence` of values with an absolute uncertainty per entry and an
information scale, and the centre and power-of-two scale of the coordinate
change. `raw_moments` keeps the moments and assigns zero uncertainty and an
infinite information scale; `standardized_moments` sets the information scale
to the standardized mass.

`standardized_moments` divides every moment by the power of two of `m0`,
centres at `c = m1/m0` with the binomial sum, and scales the order-`k` entry by
`2^(-k e)`, where `2^e` is the power of two above the spread. The spread comes
from the first even central moment that exceeds twice its uncertainty,
otherwise from `|c|`, or is one. Each entry's uncertainty is `(k+2)*epsilon`
times the sum of the absolute binomial terms: one `epsilon` for the input,
taken as exact to its last digit, and less than `(k+1)*epsilon` for the at
most `2k+1` roundings of each term in evaluating the sum. Dividing by the mass
first keeps the binomial sums finite for moments near the overflow limit.
`TestedMoments::bound` maps a support bound into the tested coordinate.

The congruence argument for invariance and the resulting formulas are in the
[design page](../realizability.md#numerical-classification).

## A — sequence builders and support dispatch

| Function | Loop/index meaning and responsibility |
|---|---|
| `hankel_matrix` | Row/column loops assign `seq[row+col]` to an `(order+1)` square matrix |
| `hankel_uncertainty` | Frobenius norm of the Hankel matrix of entry uncertainties; bounds the eigenvalue shift they can cause |
| `lower_localizing_sequence` | Entry `k` is `m[k+1]-a*m[k]`, the localizer of `x-a`; its uncertainty combines those of the two moments, and its information scale is `(1+|a|)` times that of the moments |
| `upper_localizing_sequence` | Entry `k` is `b*m[k]-m[k+1]`, the localizer of `b-x`; information scale `(1+|b|)` times that of the moments |
| `interval_localizing_sequence` | Entry `k` is `(a+b)m[k+1]-ab*m[k]-m[k+2]`, the localizer of `(x-a)(b-x)`; information scale `(1+|a+b|+|ab|)` times that of the moments |
| `extension_column` | Gathers `seq[row+matrix_order+1]` and its uncertainties, the next column requiring range consistency |

The builders assume the dispatcher's parity/order calculations make all
required indices available. Non-finite constructed matrices are detected by
`analyze_matrix`.

`check_hamburger` analyzes `H_floor(N/2)`. Odd `N` checks its next extension
column against the matrix range. For a numerically singular even `H_n`, it also
analyzes `H_(n-1)` and requires equal numerical ranks: the singular truncated
matrix must be a flat extension. This rejects PSD-but-incompatible data such as
`[1,0,0,0,1]` without changing the positive-definite or odd-order paths. A
rank has no margin against the data uncertainty, so a mismatch rejects only if
`H_n` is singular and the ranks still differ under `rank_as_given`, which
counts only the eigen-solve roundoff as zero; otherwise the result is
indeterminate.

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
and boundary measures of a few points and the incompatible even-order sequence
remain covered. `verify_shift_and_scale_preserve_classification` and
`verify_resolution_limit_never_reverses_classification` sweep location and
spread on every support.

## A — `analyze_matrix` and `classify_eigenvalue`

`analyze_matrix` builds the Hankel matrix of a sequence, scans all entries for
finiteness and the maximum absolute entry, and computes the infinity norm
relative to that maximum, avoiding an overflowing raw row sum. The user
threshold is either absolute `tolerance` or
`tolerance*entry_scale*normalized_infinity_norm`; the roundoff allowance is
`1000*epsilon*entry_scale*normalized_infinity_norm`. The matrix tolerance adds
the data uncertainty `delta` from `hankel_uncertainty` to the user threshold,
and the zero band adds it to the roundoff allowance.

The matrix is diagonalized by the [shared symmetric Jacobi
kernel](conditional_systems.md#shared-symmetric-eigenkernel) with its default
options, which scales the matrix itself. An unconverged decomposition,
non-finite eigenvalues or uncertainty, or an uncertainty that reaches the
information scale of the sequence make the matrix indeterminate. Otherwise
`classify_eigenvalue` assigns every eigenvalue to one class: negative below
minus the tolerance, zero within the zero band, positive at or above the
tolerance, ambiguous otherwise. Any negative eigenvalue makes the matrix not
PSD; otherwise any ambiguous one makes it indeterminate; otherwise any zero one
makes it a PSD boundary matrix; otherwise it is positive definite.

The ambiguous class includes **small positive** eigenvalues, so its diagnostic
says that an eigenvalue is within the numerical tolerance, without assigning a
sign; the same message names the input precision, which covers the matrices
whose uncertainty reaches the information scale. `numerical_rank` counts the positive class, and the range checks use the
zero class as the kernel, so the definiteness decision, the rank and the kernel
share one partition. `test_every_eigenvalue_is_classified` covers an
ambiguous eigenvalue that is not the smallest;
`verify_near_boundary_and_scaling` checks both signs near zero, low-density
non-diagonal matrices, and an overflowing raw row sum.

## A — `check_range_consistency`, `incorporate_range_consistency`

For a matrix not classified positive definite, the extension column must have
small projection onto every eigenvector of the zero class.
`projection_on_eigenvector` sums a column dot product; the outer loop retains
the maximum absolute projection. `kernel_projection_tolerance` sets the
threshold: absolute tolerance, or the maximum of the matrix tolerance and
`options.tolerance*vector_scale(extension)`, raised to ten times the data
uncertainty of the projection, which is the column's own uncertainty plus its
norm times the kernel rotation. The rotation is a Davis-Kahan bound, `delta`
over the gap between the zero and positive classes, with the gap estimated
from below as the smallest positive eigenvalue minus the zero band.

Projection above that threshold rejects the sequence; above `0.1*threshold`
reports indeterminate; otherwise the extension is accepted as boundary. The
factor `0.1` creates an uncertainty band, and the factor ten keeps the data
uncertainty below it. A quantitative justification for the exact width of the
band was not identified. This is not a least-squares correction: no input
moment is changed. The dispatcher checks matrix status first, so unsuccessful
matrix analyses do not reach the normal range path.

Evidence: singular ordinary/shifted extension examples in the realizability
analytical verification, and the shifted next-column cases of
`verify_shift_and_scale_preserve_classification`. The `0.1` projection band is
a fixed implementation policy, not a quantitative accuracy bound.

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
| `vector_scale`, `projection_on_eigenvector`, `euclidean_norm` | Maximum absolute extension entry, ordinary dot product and Euclidean norm used by range checks | Indirect singular-range coverage |
| `numerical_rank`, `rank_as_given` | Count of the positive class, and count of eigenvalues beyond the eigen-solve roundoff, compared by the flat-extension test | `test_flat_extension_within_data_uncertainty_is_indeterminate`, even-order kernel test |
| `RealizabilityResult::realizable`, `on_boundary` | Exact status predicates; no fresh numerical check | `test_status_helpers` |
