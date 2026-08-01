# Vandermonde Solver Selection Verification

These cases verify the user-facing `[vandermonde] solver` configuration used by
recursive CQMOM conditional moment reconstruction.

The selected analytical system is a fixed, well-conditioned 2D independent
product CQMOM case:

```text
x = {1, 3},     w_x = {0.4, 0.6}
y = {-1, 2},    w_y = {0.25, 0.75}
W_terminal = w_x * w_y
```

The parent nodes are separated and all weights are positive, so the
branch-local weighted Vandermonde solves are representative without being a
stress test. The expected terminal quadrature and raw moments are known from
the product formula, independently of the solver implementation.

## Solver Policy Versus Algorithm

`VandermondeSolver` is the user-selected policy loaded from configuration:

- `automatic`
- `dense`
- `pivoted_qr`
- `svd`
- `bjorck_pereyra`
- `bjorck_pereyra_with_dense_fallback`

The algorithm actually used is the numerical path selected after policy
resolution:

- `dense`
- `pivoted_qr`
- `svd`
- `bjorck_pereyra`

For this well-conditioned case, `automatic` is expected to resolve to
`bjorck_pereyra_with_dense_fallback` and use the Björck-Pereyra algorithm
without dense fallback. `pivoted_qr` is an explicit policy only and is expected
to solve the weighted conditional Vandermonde system directly. `svd` is also an
explicit policy only and is expected to solve the weighted conditional
Vandermonde system directly without dense fallback.

The public CQMOM result exposes branch success diagnostics, the requested
Vandermonde solver, the selected solver after policy resolution, and whether
dense fallback occurred. The internal `VandermondeAlgorithm` field identifies
the backend actually used (`dense`, `pivoted_qr`, `svd`, or `bjorck_pereyra`);
that algorithm-used field is not exposed as the same public CQMOM diagnostic
field. This verifier loads public `.qmom` files, verifies public CQMOM results,
and inspects internal recursive diagnostics for algorithm-used checks.

## Checks

The verifier:

- loads one public `.qmom` file for each supported solver policy;
- confirms configuration loading preserved the requested solver policy;
- computes public CQMOM results successfully;
- compares terminal nodes and weights against the independent analytical
  quadrature;
- reconstructs every configured raw moment from the computed quadrature;
- checks that diagnostics report the requested solver policy, selected policy,
  algorithm actually used, and fallback flag;
- compares dense and Björck-Pereyra results as an additional consistency check.

Dense/Björck-Pereyra comparison is not the only reference. Every solver policy
is also checked against the independent analytical quadrature and reconstructed
moments.

## Fallback Behavior

This verification does not force dense fallback. A stable public CQMOM case
that succeeds after deterministic dense fallback is not currently part of the
suite. Repeated-node or singular artificial systems can exercise the fallback
flag internally, but they do not provide a reliable successful user workflow
and are better kept as unit-level failure tests unless a natural public case is
identified.

## Tolerance Policy

The pass tolerance is `1e-12` absolute for terminal nodes, terminal weights,
reconstructed raw moments, and dense/Björck-Pereyra consistency. The system is
small and well-conditioned, so errors should remain near roundoff. Do not
relax this tolerance without first checking whether the configuration,
analytical reference, solver selection, or reconstruction path changed.

## Running

Configure and build the standalone test project first:

```bash
cd tests
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../install"
cmake --build build
```

Run the Vandermonde solver-selection verification:

```bash
ctest --test-dir build -R verification_vandermonde_solver_selection --output-on-failure
```
