# Verification

Each case checks a mathematical property or public workflow against a
reference built independently of the code under test. Its source header and
the README of its area state the reference, assumptions, error metric,
tolerance and reproduction; registration, checks and build commands are in
[Testing and validation](../../docs/developer/testing_and_validation.md).

| Area | What is verified | Reference |
| --- | --- | --- |
| [QMOM](qmom/README.md) | Wheeler inversion with fixed and adaptive node counts and abscissa bounds; the moment-space boundary; failures | exact discrete quadratures, classical Gaussian rules |
| [CQMOM](cqmom/README.md) | recursive reconstruction in up to six dimensions: conditioning orders, adaptive branch counts, bounds, unit invariance, moment reconstruction | prescribed product, affine and branch-dependent quadratures |
| [Distributions](distributions/README.md) | raw moments of distributions, mixtures and products; configured initialization | closed-form and high-precision moments |
| [Realizability](realizability/README.md) | support-aware classification and diagnostics | Hamburger, Stieltjes and Hausdorff examples |
| [Vandermonde](vandermonde/README.md) | the configured solvers on ordinary and clustered systems; diagnostics under rescaling | prescribed terminal quadratures |
| [Runge–Kutta](runge_kutta/README.md) | SSP-RK3, RK23 and RK45: accuracy, order, step control | closed-form ODE solutions |
| [Node advection](node_advection/README.md) | `advance_nodes(...)` for QMOM and CQMOM | closed-form node trajectories |
| [Moment source](moment_source/README.md) | moment-source updates for QMOM and CQMOM | closed-form source solutions |
| [Phase-space evolution](phase_space_evolution/README.md) | operator splitting: exactness, order, exact split solutions with time-dependent callbacks, adaptive node counts, source-step policies and retries, states from changed moments, untracked source moments, boundary contributions, native sources, failures | closed-form and manufactured solutions |
| [Output](output/README.md) | written files, read back | the values written; an analytical product rule |
| [Configuration](configuration/README.md) | loader errors for faulty `.qmom` files | one fault per fixture |

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_' --output-on-failure
```

See [Building and running](../../docs/developer/testing_and_validation.md#building-and-running).
