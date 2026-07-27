# Distribution Analytical Verification

These checks verify raw-moment generation for supported analytical
distributions and distribution-driven configuration workflows. They are
component-level verification: the expected values come from closed-form moment
formulas, algebraic identities between named distributions, mixture sums, and
independent-product tensor construction.

The cases exercise:

- one-dimensional analytical raw moments for supported distribution families;
- total number density scaling, including non-unit weights;
- invalid-parameter and overflow reporting for raw-moment generation;
- one-dimensional mixtures with absolute component weights;
- independent-product CQMOM tensors with one global number density;
- configuration parsing paths that construct QMOM and CQMOM moment states from
  distribution specifications.

Closed-form references are preferred. A small number of identity checks compare
equivalent distribution APIs, such as chi-squared as a gamma distribution and
Weibull with shape one as an exponential/gamma distribution. Those identity
checks are regression evidence for API consistency, not fully independent
analytical references.

Tolerances are near roundoff for low-order closed-form moments and relaxed only
where floating-point products, powers, or gamma-function evaluations amplify
roundoff. These tolerances should remain tied to the analytical formula and
moment scale, not to the behavior of a previous implementation.

Run with:

```bash
ctest --test-dir build -R verification_distributions_analytical --output-on-failure
```
