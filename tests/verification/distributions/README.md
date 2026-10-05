# Distribution Analytical Verification

These checks verify raw-moment generation for supported analytical
distributions and distribution-driven configuration workflows. They are
component-level verification: the expected values come from closed-form moment
formulas, algebraic identities between named distributions, mixture sums, and
independent-product tensor construction.

The cases exercise:

- one-dimensional analytical raw moments for supported distribution families,
  including chi and chi-squared against their closed forms for odd and even
  degrees of freedom, and gamma and beta with parameters up to fifty, to order
  ten;
- selected-order requests, out of sequence, against the same closed forms;
- total number density scaling, including non-unit weights;
- invalid-parameter and overflow reporting for raw-moment generation;
- one-dimensional mixtures with absolute component weights;
- independent-product CQMOM tensors with one global number density;
- configuration parsing paths that construct QMOM and CQMOM moment states from
  distribution specifications.
- a configuration-driven independent product with a normal marginal and a
  non-unit-total-mass lognormal-mixture marginal, analytical moment comparison,
  and public/default CQMOM reconstruction.

Closed-form references are preferred. A small number of identity checks compare
equivalent distribution APIs, such as chi-squared as a gamma distribution and
Weibull with shape one as an exponential/gamma distribution. Those identity
checks are regression evidence for API consistency, not fully independent
analytical references.

Tolerances are near roundoff for low-order closed-form moments and relaxed only
where floating-point products, powers, or gamma-function evaluations amplify
roundoff. The chi, chi-squared, large-parameter and selected-order cases compare
relatively, at `1e-13` for chi and chi-squared and `1e-12` for gamma and beta,
because log-gamma differences lose about `epsilon * |ln Gamma|` of relative
accuracy as the parameters grow. These tolerances should remain tied to the analytical formula and
moment scale, not to the behavior of a previous implementation.

The CTest registrations are:

- `verification_distributions_analytical_*`, one registration per distribution
  family, mixture, independent product, and configuration-driven initialization
  case.

Run all distribution analytical verification cases with:

```bash
ctest --test-dir build -R '^verification_distributions_analytical_' --output-on-failure
```
