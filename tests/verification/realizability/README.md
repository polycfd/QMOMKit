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
the other points. They are checked on every support, on the intervals `[-1, 2]` and
`[2, 5]` as well as `[0, 1]`, and for points shifted to centres 2 and 10, which
must keep their classification. On an interval, a variance above the
Bhatia-Davis bound `(b - mean)(mean - a)` must be rejected and one below it
accepted.

Boundary and near-boundary cases verify tolerance behavior. The zero measure
is reported as boundary-realizable with a `zero_measure` diagnostic. Small
positive and negative variance perturbations under a loose absolute tolerance
are checked to ensure classifications and diagnostic fields change
consistently.

Run with:

```bash
ctest --test-dir build -R verification_realizability_analytical --output-on-failure
```
