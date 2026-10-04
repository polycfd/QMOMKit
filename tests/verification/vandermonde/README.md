# Vandermonde Solver Selection Verification

These cases verify the user-facing `[vandermonde] solver` configuration used by
recursive CQMOM conditional moment reconstruction.

They cover both the ordinary well-conditioned product system below and a
separate clustered, realizable conditional system. The latter is a public CQMOM
verification against a prescribed terminal quadrature, not a comparison between
solver outputs.

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
- `moore_penrose`
- `bjorck_pereyra`

The algorithm actually used is the numerical path selected after policy
resolution:

- `dense`
- `pivoted_qr`
- `svd`
- `moore_penrose`
- `bjorck_pereyra`

For this well-conditioned case, `automatic` is expected to resolve to
`bjorck_pereyra` and use the Björck-Pereyra algorithm.
`pivoted_qr` is an explicit policy only and is expected
to solve the weighted conditional Vandermonde system directly. `svd` is also an
explicit policy only and is expected to solve the weighted conditional
Vandermonde system directly. `moore_penrose` is likewise explicit only; its
fixture sets `rcond = 1.0e-14`.

The public CQMOM result exposes branch success diagnostics, the requested
Vandermonde solver, and the selected solver after policy resolution. The
internal `VandermondeAlgorithm` field identifies
the backend actually used (`dense`, `pivoted_qr`, `svd`, `moore_penrose`, or
`bjorck_pereyra`);
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
- reconstructs the configured required CQMOM moments from the computed
  quadrature;
- checks that diagnostics report the requested solver policy, selected policy,
  and algorithm actually used;
- compares dense and Björck-Pereyra results as an additional consistency check.
  This is its own case, selected by `dense_matches_bjorck_pereyra` and
  registered as `verification_vandermonde_solver_dense_matches_bjorck_pereyra`;
  it loads only those two policies, and the per-policy checks above stay with
  the per-policy cases.

Dense/Björck-Pereyra comparison is not the only reference. Every solver policy
is also checked against the independent analytical quadrature and reconstructed
moments.

## Clustered Conditional-System Verification

`clustered_solver_policy.cpp` prescribes a fixed 2D `3 x 2` terminal rule with
conditioning order `{0,1}`. Its parent rule is

```text
x = {1, 1.007, 1.014},     w_x = {0.2, 0.3, 0.5}
```

and the branch-dependent conditional rules are

```text
x = 1.000 -> y = {-1,   2}, w_y = {0.4,  0.6}
x = 1.007 -> y = {-0.5, 3}, w_y = {0.7,  0.3}
x = 1.014 -> y = {0.25, 4}, w_y = {0.25, 0.75}
```

The test independently sums all `6 x 4` mixed moments from these six terminal
points. The configured `automatic` case in
`clustered_policy_automatic.qmom` is checked against that sum, then the same
public case is reconstructed with every public policy.

At the root conditional solve, the weighted system is

```text
A(p,i) = w_x[i] * x[i]^p,  p,i = 0,1,2.
```

An out-of-tree SVD calculation gives `kappa_2(A) = 2.33e5`; the detailed CQMOM
diagnostic reports `kappa_inf(A)` about `2.79e5`. Thus the parent geometry is
clustered enough to make the conditional solve meaningful without approaching
the effectively unstable `delta = 0.001` geometry.

Every policy must recover the independently prescribed fixed topology, nodes,
weights, and all mixed moments. `automatic` resolves to `bjorck_pereyra`. The
public diagnostics verify requested and selected policy; the internal backend
is intentionally not treated as a public contract.

The case also requires the condition estimate reported by the conditional
solve to lie between `1e5` and `1e6`, which is what makes the system clustered
enough to be a meaningful test.

The node tolerance is `5e-4`, about ten times the observed maximum error
on the selected system and below one tenth of the parent spacing. The weight
tolerance is `2e-5`, and every moment uses an absolute-plus-relative tolerance
of `1e-7 + 1e-8 max(1, |M|)`. These leave cross-platform margin while still
rejecting a materially wrong terminal rule. The independent quadrature, rather
than residual alone, is the primary oracle.

## Tolerance Policy

The pass tolerance is `1e-12` absolute for terminal nodes, terminal weights,
reconstructed raw moments, and dense/Björck-Pereyra consistency. The system is
small and well-conditioned, so errors should remain near roundoff. Do not
relax this tolerance without first checking whether the configuration,
analytical reference, solver selection, or reconstruction path changed.

## Running

Configure and build the test suite from the repository root, in Debug as
described in [Testing](../../../docs/developer/testing.md#debug-build-requirement):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build --parallel
```

Run the Vandermonde solver-selection verification, one test per policy plus
the dense/Björck-Pereyra comparison:

```bash
ctest --test-dir build -R '^verification_vandermonde_solver_' --output-on-failure
```

Run the clustered solver-policy verification:

```bash
ctest --test-dir build -R '^verification_vandermonde_clustered_' --output-on-failure
```
