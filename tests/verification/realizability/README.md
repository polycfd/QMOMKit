# Realizability Analytical Verification

These checks verify one-dimensional truncated-moment realizability
classification for the public support choices:

- Hamburger moments on the real line;
- Stieltjes moments on the nonnegative half-line;
- Hausdorff moments on a bounded interval.

The expected statuses come from standard analytical examples: distributions
whose raw moments are strictly realizable on the requested support, simple
moment sequences that violate positivity on that support, and near-boundary
second-order moment sets that exercise the configured tolerance.

The verifier also checks the public diagnostic contract. For successful
Hamburger, Stieltjes, and Hausdorff cases it verifies the checked order,
limiting matrix diagnostic, minimum pivot, tolerance used, and support-aware
message. For non-realizable cases it verifies that diagnostics identify the
public mathematical condition that failed, such as negative `m0`, a Hankel
matrix failure, a Stieltjes shifted-matrix failure, or a Hausdorff localizing
matrix failure.

Further cases follow from two standard facts: `n` points inside the support
are strictly realizable up to order `2n - 1` and on the boundary from order
`2n` on, and a point at an end of the support puts the sequence on the
boundary once the localizing matrix for that end is larger than the number of
the other points. They are checked on every support and on the intervals
`[-1, 2]` and `[2, 5]` as well as `[0, 1]`. On an interval, a variance above
the Bhatia-Davis bound `(b - mean)(mean - a)` must be rejected and one below it
accepted.

Realizability is invariant under an affine change of coordinate, so strictly
realizable, boundary and non-realizable sequences built from a location `c`
and a spread `s` must keep their classification for every `s` from `1e-9` to
`1e9` and every `c` up to the number of spreads that double-precision input
resolves for that order. Beyond that, up to `1e8` spreads, the classification
may only give way to `boundary` or `numerically_indeterminate`, never to the
opposite verdict.

Boundary and near-boundary cases verify tolerance behavior. The zero measure
is reported as boundary-realizable with a `zero_measure` diagnostic. Under a
loose relative tolerance, perturbing the fourth moment of two points at `-1`
and `1` by `1e-10` is indeterminate and by `1e-6` decided in either direction,
while a variance of `1e-10` is decided at its own scale.

Run with:

```bash
ctest --test-dir build -R verification_realizability_analytical --output-on-failure
```
