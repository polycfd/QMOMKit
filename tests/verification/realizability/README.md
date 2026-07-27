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

Boundary and near-boundary cases verify tolerance behavior. The zero measure
is reported as boundary-realizable with a `zero_measure` diagnostic. Small
positive and negative variance perturbations under a loose absolute tolerance
are checked to ensure classifications and diagnostic fields change
consistently.

Run with:

```bash
ctest --test-dir build -R verification_realizability_analytical --output-on-failure
```
