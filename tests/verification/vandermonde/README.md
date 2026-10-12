# Vandermonde solver verification

## Purpose

A CQMOM branch recovers the conditional moments at its parent nodes from a
weighted Vandermonde system, solved by the `[vandermonde] solver` the case
requests; `automatic` selects `bjorck_pereyra`
([Linear algebra](../../../docs/developer/linear_algebra.md#conditional-vandermonde-systems)).
These cases verify that the requested solver is used and reported, that every
solver recovers a prescribed quadrature, well-conditioned or clustered, and
that the detailed diagnostics describe the scaled system the solvers decide on
([Diagnostics](../../../docs/developer/linear_algebra.md#diagnostics)).

## Cases

**Product case**, `verification_vandermonde_solver_<solver>`
(`solver_selection.cpp`), for `automatic`, `dense`, `pivoted_qr`, `svd`,
`moore_penrose` and `bjorck_pereyra`. The fixtures `<solver>.qmom` differ only
in the solver (`moore_penrose.qmom` adds `rcond = 1e-14`) and hold the moments
of

```text
xi_0 = {1, 3},    w_0 = {0.4, 0.6}
xi_1 = {-1, 2},   w_1 = {0.25, 0.75},   W = w_0 w_1
```

Each case checks the loaded solver, the terminal nodes and weights against the
rule, the reconstructed against the configured moments, and the requested and
selected solver in the diagnostics.
`verification_vandermonde_solver_dense_matches_bjorck_pereyra` instead compares
the `dense` and `bjorck_pereyra` results with each other.

**Clustered case**, `verification_vandermonde_clustered_<solver>`
(`clustered_solver_policy.cpp`). The prescribed rule, conditioned on $\xi_0$,
has clustered parent nodes and a different child rule at each:

```text
xi_0 = 1.000:  xi_1 = {-1,   2},  W = 0.2 * {0.4,  0.6}
xi_0 = 1.007:  xi_1 = {-0.5, 3},  W = 0.3 * {0.7,  0.3}
xi_0 = 1.014:  xi_1 = {0.25, 4},  W = 0.5 * {0.25, 0.75}
```

The test sums the $6\times4$ raw moments of these six weighted points and
requires `clustered_policy_automatic.qmom` to match them to `1e-12`. For each
solver it compares the fixed $3\times2$ quadrature and all reconstructed
moments with the rule and checks the root diagnostics. The root system, the
Vandermonde matrix of $\xi_{0,i}/2$, has an infinity-norm condition number of
about $5.5\times10^5$; the reported estimate must lie in $[10^5, 10^6]$.

**Diagnostics invariance**, `verification_vandermonde_diagnostics_<case>`
(`diagnostics_invariance.cpp`). Divided by the power of two above the largest
node, the product case's parent nodes are $\{0.25, 0.75\}$ in every
power-of-two unit, so the scaled matrix has condition number 7.
`number_density_does_not_change_diagnostics` scales the number density by
$2^{\pm60}$ and $3.7\times10^{-9}$;
`power_of_two_unit_does_not_change_diagnostics` scales $\xi_0$ by $2^{\pm40}$,
where at $2^{40}$ the unscaled weighted matrix exceeds the default warning
threshold $10^{12}$. The estimate must stay 7, the relative residual at most
`1e-14`, and the warning off.

## Tolerances

- Product case: `1e-12` absolute on nodes, weights, moments and the
  dense/Björck–Pereyra differences; observed errors are `1e-14` or below.
- Clustered case: nodes `5e-4`, weights `2e-5`, moments
  `1e-7 + 1e-8 max(1, |M|)`. Every solver returns nodes of order `1e-5` and
  weights of order `1e-6` from the rule while reproducing the moments to
  about `1e-11`. The error is identical for every solver: it comes from the
  clustered parent nodes, which the moments pin only weakly, not from solver
  rounding. The node bound leaves a factor of about ten and stays below a
  tenth of the parent spacing; the moment bound is a safety bound.
- Diagnostics: the estimate matches 7 to a relative `1e-13`, the rounding of
  the reconstructed nodes, exactly under power-of-two scaling, and to `1e-12`
  at number density $3.7\times10^{-9}$.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_vandermonde_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing_and_validation.md#building-and-running).
