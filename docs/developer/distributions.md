# Distributions

The `compute_raw_moments` overloads in
[`src/distributions.cpp`](../../src/distributions.cpp) generate raw moments of
analytical distributions, so that a case can start from a distribution rather
than from given moments. A single distribution with absolute number density
$n_0$ (`total_weight`, default one) has $m_k = n_0\,\mathrm{E}[X^k]$, so
$m_0 = n_0$; nothing is normalized. Reconstruction never learns where its
moments came from. The families and their parameters are listed in the
[Reference](../user/7_reference.md#distribution-families).

## Families

| Family | Evaluation of $\mathrm{E}[X^k]$ |
| --- | --- |
| Normal | Recurrence $\mu\,\mathrm{E}[X^{k-1}] + (k-1)\,\sigma^2\,\mathrm{E}[X^{k-2}]$; the integer factor is applied before $\sigma^2$, so an exact zero cannot become $0\cdot\infty$. |
| Log-normal | $\exp(k\mu + k^2\sigma^2/2)$, in log form. |
| Weibull | $\lambda^k\,\Gamma(1 + k/c)$, in log form with one `std::lgamma` per order. |
| Gamma, beta, chi-squared, chi | Products of one factor per order, below. |

The log forms add $\ln n_0$ and check the range before calling `std::exp`.

### Products of simple factors

```math
\begin{aligned}
\text{gamma } (a, \theta):\quad & \mathrm{E}[X^k] = \prod_{i=0}^{k-1} \theta\,(a + i), &
\text{beta } (\alpha, \beta):\quad & \mathrm{E}[X^k] = \prod_{i=0}^{k-1} \frac{\alpha + i}{\alpha + \beta + i}, \\
\text{chi-squared } (\nu):\quad & \mathrm{E}[X^k] = \prod_{i=0}^{k-1} (\nu + 2i), &
\text{chi } (\nu):\quad & \mathrm{E}[X^k] = (\nu + k - 2)\,\mathrm{E}[X^{k-2}].
\end{aligned}
```

`ScaledProduct` holds $n_0$ times the product as a mantissa and a binary
exponent, renormalized with `std::frexp` after every factor. No intermediate
can overflow or underflow, and the final exponent decides exactly whether the
moment fits in a double; subnormal results are returned as they are. Each
factor adds a few roundings, so a moment of order $k$ is accurate to a few
times $k$ units in the last place however large the parameters are, where a
difference of log-gamma values would lose accuracy in proportion to
$|\ln\Gamma|$.

The maximum-order path extends one product order by order
(`compute_from_products`; chi keeps separate even and odd chains), while the
selected path builds each requested order from the start
(`evaluate_product_moment`, `evaluate_chi_moment`). Both multiply the same
factors in the same order, so they return identical values, as do the other
families.

### First moment of the chi distribution

Odd chi moments start from $\mathrm{E}[X] = \sqrt2\,r(\nu/2)$ with
$r(x) = \Gamma(x+\tfrac12)/\Gamma(x)$, which has no finite product form.
`half_step_gamma_ratio` raises $x$ to at least 10 with the exact recurrence
$r(x) = r(x+1)\,x/(x+\tfrac12)$, then sums eight terms of the asymptotic series

```math
\ln\frac{r(x)}{\sqrt{x}} \sim \sum_{n\ \mathrm{odd}} \frac{(2^{-n} - 2)\,B_{n+1}}{n(n+1)\,x^n}
= -\frac{1}{8x} + \frac{1}{192\,x^3} - \frac{1}{640\,x^5} + \cdots,
```

with Bernoulli numbers $B_m$. From $x = 10$ the truncation error is below
$10^{-17}$. The series follows from the expansions of $\ln\Gamma(x+h)$ at
$h=\tfrac12$ and $h=0$ [[DLMF, 5.11.8](https://dlmf.nist.gov/5.11.E8)] with
$B_m(\tfrac12)=-(1-2^{1-m})B_m$ [[DLMF, 24.4.27](https://dlmf.nist.gov/24.4.E27)].

## Mixtures

A `DistributionMixture` has $m_k = \sum_r n_r\,\mathrm{E}_r[X^k]$ with absolute
component densities $n_r$, so $m_0 = \sum_r n_r$. It needs at least one
component, every density finite and non-negative, and at least one positive.
Components of zero density are skipped before their parameters are evaluated.
Each other component is evaluated with $n_r$ as its total weight and added per
order with compensated summation (`add_compensated`): a running correction
per order collects the rounding error of each addition, recovered with the
larger operand first, and is added at the end.

## Independent products

An `IndependentProductDistribution` has one marginal per phase-space
dimension, a distribution or a mixture, and one number density $n_0$:

```math
M_{k_0,\ldots,k_{D-1}} = n_0 \prod_{d} \frac{m^{(d)}_{k_d}}{m^{(d)}_0}.
```

Each marginal is normalized by its own zeroth moment, so mixture densities
inside a marginal set only proportions, and $n_0$ is applied once. The factors
are multiplied directly rather than in log form because normal marginals have
zero or negative odd moments. A zero $n_0$ returns zeros without evaluating
the marginals.

## Full tensors and selected moments

The tensor overload returns every entry up to inclusive maximum orders, with
extents `maximum_orders[d] + 1`.

A `RawMomentRequest` lists phase-space-order multi-indices; the
`SelectedRawMomentResult` keeps request order and duplicates, one value per
index. One-dimensional requests have rank one. For a product, the selected
overload collects the distinct orders per dimension, plus order zero,
evaluates each marginal once through its own selected overload, and assembles
the requested moments. CQMOM initialization from a configuration file requests
exactly the required moments of the layout and converts them with
`create_cqmom_moment_data`; see [CQMOM](cqmom.md#required-moments) and
[Configuration and output](configuration_and_output.md).

## Validation and failure

| `RawMomentStatus` | Cause |
| --- | --- |
| `invalid_input` | A parameter outside its range, a negative or non-finite weight or density, an invalid mixture, a request of the wrong rank, or an order equal to the largest `std::size_t`. |
| `overflow` | A moment is infinite. |
| `underflow` | A moment whose exact value is nonzero rounds to zero; exact zeros, such as odd moments of a zero-mean normal, are not underflows. |
| `non_finite` | Any other non-finite value. |

A single distribution checks its parameters before the zero-weight shortcut,
which then returns zeros. A parameter that is representable but makes moments overflow from some
order on, such as a normal or log-normal variance that overflows, fails only
requests that reach that order.

A failed result carries no moments, only the status, a message and the
location where known: `failed_order` and `failed_component` in
`RawMomentResult`; `failed_component`, `failed_dimension` and the multi-index
`failed_order` in `SelectedRawMomentResult`; `failed_dimension` and
`failed_order` in `TensorMomentResult`. Component and dimension fields that
do not apply hold `invalid_component_index` and `invalid_dimension_index`.

## Adding a family

1. Add a parameter struct to
   [`include/qmomkit/distributions.hpp`](../../include/qmomkit/distributions.hpp)
   and append it to `DistributionSpecification`.
2. Add `compute_<family>` and `compute_<family>_selected` to
   [`src/distributions.cpp`](../../src/distributions.cpp) and route them in
   `RawMomentVisitor` and `SelectedRawMomentVisitor`. Prefer a product of
   simple factors through `ScaledProduct`, otherwise a log form through
   `evaluate_log_expectation`, and give both paths the same arithmetic per
   order.
3. Validate the parameters first and return `invalid_input` with a message.
4. Add the family name to `parse_distribution` and its parameter names to
   `distribution_parameters` in
   [`src/configuration.cpp`](../../src/configuration.cpp), and document them
   in [`docs/user/7_reference.md`](../user/7_reference.md).
5. Add a sampler and a parameter check to the Monte Carlo companion in
   [`tools/monte_carlo/src/sampling.cpp`](../../tools/monte_carlo/src/sampling.cpp),
   which visits every family.
6. Test validation, statuses and selected against full results in
   `tests/unit/distributions.cpp`, parsing in `tests/unit/configuration.cpp`,
   and the moments against independent analytical references in
   `tests/verification/distributions/`.
