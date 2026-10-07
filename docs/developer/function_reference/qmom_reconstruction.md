# QMOM reconstruction and coordinate scaling

## Ownership and call context

Public `qmomkit::compute_qmom` is declared in
[`qmom.hpp`](../../../include/qmomkit/qmom.hpp) and defined in
[`public_api.cpp`](../../../src/public_api.cpp). It dispatches to local
`compute_qmom_with_automatic_scaling` → `compute_qmom_unscaled` → public
`compute_quadrature`, or the shared bounded-candidate selector.
`compute_quadrature` is declared in
[`wheeler.hpp`](../../../include/qmomkit/wheeler.hpp); its numerical helpers
are local to [`wheeler.cpp`](../../../src/wheeler.cpp). CQMOM calls the same
inversion for each branch through `internal::compute_conditional_quadrature`,
which passes the scale of the coordinate for conditional levels. Source evolution uses a separate
[current-count/growth selection](phase_space_evolution.md), not a new Wheeler
method. See [Wheeler design](../qmom_wheeler.md) for recurrence/eigenpair
correspondence.

## A — `compute_quadrature`

This function turns the first `2*max_nodes` raw moments into an accepted rule,
possibly with fewer nodes. `max_nodes` is a request ceiling, not a computed
numerical rank. The returned count is `weights.size()`; the low-level API has
no separate supported-count field.

It checks finite nonnegative `rmin` and `eabs`, finite positive `cutoff`,
positive node count, sufficient moments, and finiteness of the **required
prefix**. Extra moments are not checked here. Checked size arithmetic rejects a
`max_nodes` value for which `2*max_nodes` is unrepresentable before dependent
storage is constructed. The implementation then distinguishes:

1. Negative mass: `negative_mass_result` returns `not_realizable`, a diagnostic
   weight `abs(m0)`, and `m1/m0`. This helper does not check that ratio for overflow.
2. `m0 < cutoff`: failure with placeholder node `0` and weight `1`; this is not a
   mass-preserving rule.
3. One requested node: `one_node_result` accepts `w=m0`, `x=m1/m0` if the ratio
   is finite. `rmin` has no absolute-mass role.
4. Fixed reconstruction: one `attempt_quadrature_count` call; no order retry.
5. Adaptive reconstruction: descending candidate counts from the ceiling to two.
   A failed attempt is retried **only** when it reports a sufficiently negative
   sigma diagonal. Other failed attempts return immediately. Exhaustion accepts
   the one-node rule if its ratio is finite.

Thus adaptive success can discard information in higher moments, including a
negative-variance requested two-node sequence. It is not a full-sequence
realizability certificate. This behavior is intentional in the current tests,
not a newly inferred repair rule.

Evidence: [`tests/unit/wheeler.cpp`](../../../tests/unit/wheeler.cpp),
especially
`test_adaptive_requested_two_node_negative_variance_recovers_one_node`,
`test_adaptive_requested_three_node_high_rank_sigma_recovers_two_nodes`,
`test_fixed_count_negative_variance_remains_strict_failure`, and invalid-input
cases.
[`wheeler_analytical.cpp`](../../../tests/verification/qmom/wheeler_analytical.cpp)
and [`adaptive_qmom.cpp`](../../../tests/verification/qmom/adaptive_qmom.cpp)
compare analytical/discrete moments and recovered rules. These are finite
cases, not a characterization of all reduction boundaries.

## A — `compute_recurrence`, `make_jacobi_matrix`, `compute_weights`

`compute_recurrence` constructs Wheeler's sigma table from raw moments. The
first loop writes `sigma(1,i)=m[i-1]`; the outer `k` loop advances polynomial
degree, and the inner `l` loop computes the still-needed triangular strip from
the previous two rows. Each new `a` is a difference of sigma ratios; `b` is a
diagonal ratio. These divisions occur before the caller's finiteness/degeneracy
checks. A zero sigma denominator can therefore produce non-finite intermediate
coefficients. Adaptive handling deliberately examines sigma first and then
recomputes the reduced recurrence; fixed handling rejects non-finite recurrence
coefficients.

`make_jacobi_matrix` places `a` on the diagonal and `sqrt(b[i])` on both
adjacent diagonals, for the current candidate prefix. Its caller checks
negative `b` before entry. The [shared
eigenkernel](conditional_systems.md#shared-symmetric-eigenkernel) diagonalizes
it with its default options; an unconverged decomposition returns
`not_realizable` with the message "Jacobi matrix eigen solve did not
converge.". `compute_weights` loops over eigenvector columns and forms
`m0*v(0,i)^2`. Node sorting and corresponding eigenvector permutation belong to
the shared kernel. Exact arithmetic explains the intended mass sum; the source
has no final mass or full moment-residual acceptance check.

Candidate sensitivities are cancellation in sigma recurrence, tiny diagonal
ratios and raw moment dynamic range.
`test_adaptive_exact_degenerate_sequence_reduces_before_finite_check` and the
large-abscissa/roundoff variants in `tests/unit/wheeler.cpp` specifically cover
reduction before rejecting non-finite higher recurrence data. The analytical
Wheeler verification establishes selected exact rules, not eigen convergence at
every requested order.

## A — `attempt_quadrature_count`, `adaptive_sigma_tolerances`, `adaptive_accepts`

The attempt owns three distinct loops. First, adaptive sigma inspection
descends from the candidate ceiling to degree two. The threshold for
`sigma(k,k)` is `100*epsilon*max(M(k,k), m0*s^(2(k-1)))`: `M` runs the sigma
recurrence again with every term in absolute value, a bound on its rounding in
the same units, and `s` is the coordinate scale passed for conditional moments
(zero for data). A sigma below the
negative threshold marks an attempt for outer retry. An absolute sigma within
threshold sets a lower count, with an immediate successful one-node return at
count one. The recurrence is recomputed at the surviving count and must be
finite with nonnegative `b`.

Second, candidate Jacobi matrices descend from that surviving count to two,
using prefixes of the already computed recurrence. Non-finite eigenvalues or
weights fail immediately. Fixed mode accepts the first finite result. Adaptive
mode calls `adaptive_accepts`; exhaustion returns a successful one-node
payload. The one-node fallthrough does not increment the adaptive-reduction
counter on all paths; counters measure the implemented recording sites, not a
guaranteed one-event-per-count-change invariant.

Third, `adaptive_accepts` checks `min(w)/max(w) > rmin`. For candidate counts
of three or more it scans each unordered node pair once to recover minimum pair
distance and maximum pair distance and requires minimum distance/span `> eabs`;
a non-finite or non-positive span rejects. This is translation and
positive-scale invariant.

A two-node candidate has only one pairwise distance, so normalizing it by the
maximum pairwise distance (as above) is always exactly one and cannot
discriminate separation. Instead the two-node case normalizes the pairwise
distance by the larger abscissa magnitude and requires `|a0 - a1| / max(|a0|,
|a1|) > eabs`; a non-finite distance/magnitude or a non-positive magnitude
rejects. This form is positive-scale invariant but not translation invariant: a
pair with a fixed absolute separation is more likely to be rejected (collapsed
to one node) the larger its abscissas' magnitude, since the same absolute gap
becomes relatively smaller. The caller supplies at least two entries;
empty/singleton vectors are not valid inputs to this local helper.

Targeted tests in `tests/unit/wheeler.cpp` exercise offset versus symmetric
pairs, three-node span, small weight, duplicates, sigma retry, and threshold
one-node fallback, and the independence of the node count from the unit and
the number density. Number density affects only the absolute `m0 < cutoff`
gate.

## Numerical limits

Raw moments become poorly conditioned under large coordinate translations.
Automatic scaling rescales coordinates but does not center the distribution, so
it cannot recover shape information already lost to cancellation in translated
raw moments. Near-coincident support is likewise limited by finite precision in
the recurrence and Hankel information. `cutoff` is the absolute low-population
gate; `rmin` is a relative quadrature-weight criterion. These controls serve
different purposes and neither is a general realizability certificate.

## A — `internal::select_bounded_quadrature`

Defined in `src/wheeler.cpp`, declared in
[`bounded_quadrature.hpp`](../../../src/internal/bounded_quadrature.hpp), this
shared selector accepts nodes against `AbscissaBounds`. It is called by
ordinary QMOM, recursive CQMOM, and QMOM post-source selection.

The local `evaluate` lambda checks candidate count, available prefix, finite
ordered bounds, then inverts exactly that prefix. Success requires an exact
returned count and all nodes in `[lower-1e-12,upper+1e-12]`. It does not clip
nodes or invoke `check_realizability`. Adaptive Wheeler weight/separation
checks still apply because the supplied Wheeler options retain adaptivity.

Fixed selection evaluates the configured count once. Adaptive selection
descends to one, keeping the first rejected candidate for reduction
diagnostics. The first exact-count admissible rule wins. The descent starts at
the initial active/configured ceiling, or at the largest count the moment
sequence can feed when that ceiling is higher, since every count above half the
sequence length is rejected without an inversion; clamping selects the same
rule in a number of attempts bounded by the available moments rather than by
the requested count. At least one candidate is always evaluated, so a sequence
too short even for one node still reports the evaluator's failure. A clamped
start can therefore be reduced relative to the request and still be accepted on
its first attempt, with no rejected candidate to report. A permitted growth trial requires previous active metadata, adaptive
policies, no current count loss, and spare configured capacity; it evaluates
active+1 once. Rejected growth preserves the accepted current rule. Diagnostics
may describe a rejected candidate while returned quadrature describes the
accepted fallback; these must not be conflated. Both branches require a positive
initial and configured count; the selector rejects anything else as invalid
input before evaluating a candidate, so the adaptive descent never returns a
default-constructed selection whose quadrature would still read as successful.

Evidence: bounded tests in
[`tests/unit/recursive_quadrature.cpp`](../../../tests/unit/recursive_quadrature.cpp),
including
`test_recursive_bounded_support_returned_node_tolerance_and_no_clipping`,
weight/separation rejection, and rejected growth; public verification under
[`qmom/quadrature_bounds.cpp`](../../../tests/verification/qmom/quadrature_bounds.cpp)
and
[`cqmom/quadrature_bounds.cpp`](../../../tests/verification/cqmom/quadrature_bounds.cpp).
Bounds are acceptance conditions in the coordinates actually inverted, so
scaling changes the physical width corresponding to the absolute `1e-12`
margin.

## A — scaling preparation and restoration

[`reconstruction_scaling.cpp`](../../../src/reconstruction_scaling.cpp)
implements private declarations in
[`reconstruction_scaling.hpp`](../../../src/internal/reconstruction_scaling.hpp).
`make_qmom_scaling_context` and `make_cqmom_scaling_context` prepare temporary
coordinate changes; they do not divide all moments by mass.

Local `automatic_scale` prefers `sqrt(m2)/sqrt(m0)` for positive usable values,
then `abs(m1)/m0`, then `1`. CQMOM chooses one scale per **global coordinate
marginal**, not a separate scale for every conditional branch. Missing selected
scale moments supply zero to this selection. `scale_value` loops over integer
order dividing in `long double`; `scale_mixed_value` nests coordinate and power
loops to divide by the product of coordinate scales. Order zero is unchanged.
`scale_bounds` applies the same coordinate transformation. `convert_to_double`
rejects non-finite/out-of-double-range values and nonzero values becoming zero.
Failure produces `nullopt` for the complete preparation.

`restore_abscissas`, `restore_cqmom_abscissas`, and `restore_bounds` multiply
in extended precision and use the same representability checks. Restoration may
partially change its temporary argument before returning false. Public wrappers
then discard that temporary result and repeat unscaled reconstruction.
`can_restore_cqmom_result` checks copies of all terminal coordinates and
retained bounded diagnostics before public conversion. The fallback is for
unsafe transformation/restoration; a safely restored **failed numerical
inversion** is returned as failure, not retried unscaled merely to obtain
success.

Tests in [`tests/unit/public_api.cpp`](../../../tests/unit/public_api.cpp)
check physical restoration, available scale moments, covariance, and
unsafe-preparation fallback.
[`cqmom/automatic_scaling.cpp`](../../../tests/verification/cqmom/automatic_scaling.cpp)
uses prescribed ragged rules and rescaled tensors, including a case where
enabling scaling changes recovery. This does not establish shift invariance or
robustness for arbitrary marginal scales or long-double implementations.

## B/C — entry and payload helpers

| Functions | Purpose and important semantics | Evidence |
|---|---|---|
| `compute_qmom` overloads; `compute_qmom_unscaled`; `compute_qmom_with_automatic_scaling` | Case unpacking, optional transform, bounded/unbounded dispatch; output status maps without reclassifying realizability | `tests/unit/public_api.cpp` |
| `to_wheeler_options`, `to_qmom_status` | Preserve exposed thresholds/adaptivity and map enum values; automatic scaling belongs outside Wheeler | Definitions in `src/public_api.cpp` |
| `one_node_result`, `negative_mass_result` | Status-bearing diagnostic or accepted payloads as distinguished above | Negative mass/variance tests in `tests/unit/wheeler.cpp` |
| `recurrence_is_finite`, `has_negative_recurrence_b`, `required_moments_are_finite` | Numerical guards; the latter checks a prefix, not every stored moment | Same Wheeler tests |
| `finite_positive`, `convert_to_double` | Transform selection and representability guards, including nonzero underflow | Scaling tests in `tests/unit/public_api.cpp` |

The elementary `MomentSet` methods, defined in their own `src/moment_set.cpp`
translation unit, are covered in
[state and representations](state_and_representations.md).
