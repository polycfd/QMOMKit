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
- gamma, chi-squared, beta and even chi moments for parameters up to `1e8`,
  against the same products multiplied out in error-free double-double
  arithmetic, and odd chi moments against values computed with mpmath;
- selected-order requests, out of sequence, against the same references;
- inversion of the moments of narrow gamma distributions, up to shape `1e7`;
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
roundoff. The gamma, beta, chi and chi-squared moments are products of one
simple factor per order, accurate to a few rounding errors whatever the
parameters, so their cases compare relatively: at `4e-15` against closed-form
products, and at `4 (k + 1) epsilon` for order `k` against the error-free
references. These tolerances should remain tied to the analytical formula and
moment scale, not to the behavior of a previous implementation.

The CTest registrations are:

- `verification_distributions_analytical_*`, one registration per distribution
  family, mixture, independent product, and configuration-driven initialization
  case.

Run all distribution analytical verification cases with:

```bash
ctest --test-dir build -R '^verification_distributions_analytical_' --output-on-failure
```
