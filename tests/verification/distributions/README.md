# Distribution moment verification

## Purpose

These cases verify the raw moments that `compute_raw_moments` and `.qmom`
initial conditions produce for named distributions, mixtures and independent
products, against the contracts of
[Distributions](../../../docs/developer/distributions.md): $m_0$ is the number
density and is not normalized; mixture components carry absolute number
densities; each marginal of an independent product is normalized by its own
$m_0$ and the product's number density is applied once; and gamma, beta, chi
and chi-squared moments are accurate to a few rounding errors for any
parameters.

## Cases

All are `verification_distributions_analytical_<case>`
(`distributions_analytical.cpp`, whose header lists the closed forms).

- **Closed forms of single families**: `standard_normal`,
  `shifted_scaled_normal`, `log_normal`, `beta`, `gamma`, `chi`, `weibull`;
  the products of one factor per order `chi_closed_form`,
  `chi_squared_closed_form`, `gamma_large_shape` and `beta_large_parameters`
  (parameters up to 50, order 10); and `selected_orders_closed_form`, with
  orders requested out of sequence.
- **Large parameters**: `large_parameter_moments_against_exact_products`
  compares gamma, chi-squared, beta and even chi moments for parameters up to
  $10^8$ with the same products formed in error-free double-double arithmetic.
  `chi_odd_orders_against_high_precision_values` compares odd chi moments with
  first moments computed with mpmath at 40 digits; the command is in the
  source.
- **Inversion**: `narrow_gamma_moments_invert` inverts gamma moments of shape
  up to $10^7$ with fixed-count Wheeler inversion and requires the rule to
  rebuild them; where the moments no longer determine the rule (shape $10^7$
  with three nodes, $10^5$ with four) it requires `not_realizable`.
- **Mixtures**, against the sum of component number densities times closed
  forms: `bimodal_normal_mixture`, `heterogeneous_gamma_log_normal_mixture`,
  `three_component_multimodal_mixture`, `strongly_unequal_component_densities`
  ($10^{12} : 1$), `normal_odd_moment_cancellation` (normals at $\pm 10$),
  `large_total_density_is_not_normalized`, and
  `zero_density_mixture_component_policy` (a zero-density component with
  invalid parameters is skipped).
- **Independent products**, against the number density times the product of
  normalized marginal closed forms, including the flattened index layout:
  `independent_product_{normal_gamma,beta_log_normal_weibull,signed_odd_moments}`
  and, configured,
  `configuration_independent_product_{two_dimensional,normal_gamma,beta_log_normal_weibull,dimension_sweep}`,
  the last in four to six dimensions.
  `configuration_independent_product_mixture_marginal` has a lognormal-mixture
  marginal with $m_0 = 4.75$, which must be divided out; its prepared case is
  then reconstructed with `compute_cqmom` and must rebuild its moments.
- **Failure reporting**: `maximum_order_zero_and_invalid_parameters`,
  `invalid_mixture_density`, `mixture_accumulation_overflow` and
  `independent_product_tensor_overflow` check the status and the failed
  component or order. `independent_product_zero_density` must return zeros
  without evaluating its invalid marginals.

The following are consistency checks between two APIs, not independent
references: `chi_squared` and `weibull_gamma_identity` (against the gamma
family), `one_component_mixture_equivalence`,
`component_order_permutation_invariance`,
`independent_product_one_dimensional_equivalence`, and
`configuration_distribution_initialization` and
`configuration_mixture_initialization` (configured moments against direct
`compute_raw_moments`).

## Tolerances

Errors are relative to $\max(|\text{reference}|, 1)$ unless stated otherwise.

- Single families: $10^{-14}$ to $10^{-13}$; observed $10^{-15}$ or below.
- Products of one factor per order and selected orders: purely relative
  $4 \times 10^{-15}$, observed below $10^{-15}$, since library and reference
  each add a few rounding errors.
- Error-free and high-precision references: purely relative
  $4(k + 1)\epsilon$ at order $k$, observed about $2\epsilon$.
- Mixtures, products and configured moments: $10^{-13}$ to $10^{-12}$, with
  $10^{-11}$ absolute for tensors with lognormal and Weibull factors; observed
  $10^{-15}$ or below. The odd moments of the normals at $\pm 10$ cancel to
  zero and are checked at $10^{-12}$ ($10^{-14}$ for $m_1$); they are observed
  exactly zero.
- Narrow gamma inversion: purely relative $10^{-13}$, observed $10^{-15}$.
- CQMOM reconstruction: $5 \times 10^{-10}$, a safety bound; observed errors
  are at rounding level.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_distributions_analytical_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing_and_validation.md#building-and-running).
