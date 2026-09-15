# Moment-Source Analytical Verification

These cases exercise the public moment-source update,
`qmomkit::advance_qmom_moments_with_source(...)`, which advances raw moments through
a caller-supplied source term. The verification workflow then separately checks
final QMOM/CQMOM inversion and raw-moment reconstruction so the evolved moments
are tested as reconstructable moment states. Initial moment states and
QMOM/CQMOM reconstruction settings are loaded from `.qmom` files. Source
coefficients and SSP-RK3 step plans live in the verifier because the current
configuration format does not define source-model parameters or SSP-RK3
time-integration options.

These are source-operator tests only. They do not include advection node
advection or full phase-space evolution composition.

`moment_source_analytical.cpp` is the canonical verification driver for this
area. It uses named `.qmom` configurations, closed-form source solutions, and
reported error metrics.

Direct source-update and diagnostic regression checks live under
`tests/component_regression/moment_source.cpp`.

## Order-Dependent Relaxation

Configuration: `order_dependent_relaxation_qmom.qmom`

Initial one-dimensional QMOM moments come from:

```text
nodes   = {-0.6, 0.4, 1.3}
weights = {0.25, 0.35, 0.40}
```

The target moments come from:

```text
nodes   = {-0.4, 0.7, 1.5}
weights = {0.25, 0.35, 0.40}
```

For `k > 0`,

```text
dm_k/dt = -0.15 k (m_k - m_k*)
m_k(t0 + tau) = m_k* + (m_k(t0) - m_k*) exp(-0.15 k tau)
```

For `k = 0`, `lambda_0 = 0`, so total mass is conserved exactly by the model.
The verifier uses moments `m0` through `m5`, runs one step from `t=0.2` with
`dt=0.08`, and runs six outer steps from `t=0.2` with `dt=0.04`.

## Common-Rate CQMOM Relaxation

Configuration: `common_rate_relaxation_cqmom.qmom`

The initial and target states are exact two-dimensional terminal quadratures on
the same support with different positive weights. The initial support is:

```text
(x, y) = (-0.5, -1.0), W = 0.20
(x, y) = (-0.5,  0.8), W = 0.30
(x, y) = ( 1.2, -1.0), W = 0.15
(x, y) = ( 1.2,  0.8), W = 0.35
```

The target weights are:

```text
0.25, 0.20, 0.20, 0.35
```

For every nonzero multi-index `alpha`,

```text
dM_alpha/dt = -0.45 (M_alpha - M_alpha*)
M_alpha(t0 + tau) = M_alpha* + (M_alpha(t0) - M_alpha*) exp(-0.45 tau)
```

For `alpha = (0,0)`, the derivative is zero and total mass remains unchanged.
The configuration uses `conditioning_order = {1,0}` with fixed rectangular
storage. The verifier runs one step from `t=0.1` with `dt=0.06`, and five
outer steps from `t=0.1` with `dt=0.06`.

## Linear Coupled Source

Configuration: `linear_coupled_source_qmom.qmom`

Initial one-dimensional QMOM moments come from:

```text
nodes   = {0.5, 1.5}
weights = {0.45, 0.55}
```

The source is lower triangular:

```text
dm0/dt = 0
dm1/dt = -0.2 m1
dm2/dt = 0.15 m1 + 0.1 m2
dm3/dt = 0.6 m3
```

The exact solution is:

```text
m0(t0 + tau) = m0(t0)
m1(t0 + tau) = m1(t0) exp(-0.2 tau)
m2(t0 + tau) = m2(t0) exp(0.1 tau)
             + 0.15 m1(t0) (exp(-0.2 tau) - exp(0.1 tau)) / (-0.2 - 0.1)
m3(t0 + tau) = m3(t0) exp(0.6 tau)
```

The `m2` equation couples to `m1`, so this case guards against accidental
component-wise source assumptions. The verifier runs one step from `t=-0.2`
with `dt=0.08`, and eight outer steps from `t=-0.2` with `dt=0.06`. The chosen
coefficients preserve realizability over the tested interval, so the final
two-node QMOM inversion remains meaningful.

The convergence study retains the same closed-form coupled source problem and
advances to the same final time, `tau = 0.48`, at four timestep refinement
levels. It computes the observed SSP-RK3 order from adjacent final-moment
errors. The two finest observed orders must satisfy `p >= 2.8`, explicitly
verifying third-order temporal behavior of the public moment-source wrapper.

## Invariants And Tolerances

Every case checks the source update and then performs separate reconstruction
checks:

- every evolved raw moment against the closed-form solution;
- conservation of `m0` or `M[0,0]` where prescribed;
- successful SSP-RK3 diagnostics and accepted timestep reporting;
- final QMOM or CQMOM inversion;
- reconstructed raw moments from the final computed quadrature state, with the
  exact reference state also inverted as a realizability cross-check.

Near-zero quantities use absolute-error protection through a unit-floored
relative-error denominator. The pass tolerances are:

```text
evolved moment absolute/relative error <= 2e-6
reconstruction absolute/relative error <= 2e-8
zeroth-moment conservation error <= 2e-8
```

The moment tolerance reflects fixed-step SSP-RK3 truncation error against
advection exact exponentials. It was selected after keeping the time steps
small enough that the analytical discrepancy is several orders below the
source magnitudes.

## Running

Configure and build the standalone test project first:

```bash
cd tests
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../install"
cmake --build build
```

Run all source verification cases:

```bash
ctest --test-dir build -R verification_moment_source_ --output-on-failure
```

Run individual cases:

```bash
ctest --test-dir build -R verification_moment_source_order_dependent_relaxation --output-on-failure
ctest --test-dir build -R verification_moment_source_common_rate_relaxation --output-on-failure
ctest --test-dir build -R verification_moment_source_linear_coupled --output-on-failure
```
