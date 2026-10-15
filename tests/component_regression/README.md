# Component regression tests

## Purpose

These tests pin behavior a change could break while verification still
passes: operator results and diagnostics, status mapping, agreement between
`CqmomCase` and the default `CqmomPreparedCase`, and recorded CQMOM results
([Testing](../../docs/developer/testing.md#the-test-tree),
[CQMOM](../../docs/developer/cqmom.md),
[Phase-space evolution](../../docs/developer/phase_space_evolution.md)).

## Cases

All names start with `component_regression_`.

- `moment_source_*` (3): `advance_qmom_moments_with_source(...)` with constant,
  linear and coupled sources, against the exact one-step SSP-RK3 result.
- `node_advection_*` (6): `advance_nodes(...)` on QMOM and 2×2 CQMOM
  quadratures, against closed-form trajectories; weights unchanged, step
  counts, and a tighter RK23 tolerance landing closer to the exact node.
- `phase_space_evolution_*` (3): zero rates keep the moments; one node under
  constant velocity and linear source matches the exact translation and the
  SSP-RK3 amplification factor; a Strang-split CQMOM state matches
  `advance_nodes(...)`, so dimension order is kept.
- `cqmom_compute_*` (8): $M_{0,0} < 0$ gives `qmom_failed` through both case
  types; without a realizability policy, a bounded reconstruction rejects only
  the node Wheeler returns outside $[0, 1]$; extraction throws on structural
  mistakes; `CqmomCase` reports `invalid_input` for bad moment data.
- `cqmom_prepared_case_moment_replacement_*` (2): replaced moments match parsed
  ones, and `CqmomMomentData` matches a full tensor, in settings, quadrature,
  diagnostics, state and one evolution step.
- `cqmom_drag_and_degenerate_moments_*` (3): moments of eight prescribed
  three-dimensional nodes, after drag or with two velocity nodes
  $1.5 \times 10^{-6}$ apart, reconstruct and reproduce every required moment,
  and the drag case also recovers its nodes; one branch's conditional sequence
  of the drag case also inverts on its own with 1D QMOM.
- `cqmom_physical_transport`: a fixture of unknown provenance, so it pins
  behavior, not correctness. Adaptive reconstruction collapses each terminal
  branch to one node, and zero-rate evolution keeps those counts.
- `cqmom_reference_*` (20): the matrix of `cqmom_case_matrix.*` over
  dimension, conditioning order, node counts, adaptive branches, bounds, every
  Vandermonde solver and every splitting scheme. Static cases must recover the
  nodes they are built from; all must match their recorded values, with both
  case types agreeing bit for bit. The solver cases agree to about $10^{-15}$
  relative, so solver selection is verified in
  [`vandermonde/`](../verification/vandermonde/README.md) instead.

`cqmom_reference_values.inc` has no shipped generator; changing it needs the
justification of any numerical change, and its source header says how to
record it.

## Tolerances

| Cases | Bound | Observed | Why |
| --- | --- | --- | --- |
| moment source | $10^{-15}$ | $10^{-16}$ | round-off of one exact step |
| SSP-RK3 linear velocity | $2 \times 10^{-7}$, $3 \times 10^{-9}$ | $10^{-7}$, $10^{-9}$ | truncation error of one step |
| RK23 | $2 \times 10^{-8}$, $10^{-7}$ | $5 \times 10^{-9}$, $4 \times 10^{-8}$ | relative tolerance $10^{-8}$ |
| recorded values | $10^{-12}$ relative, $10^{-15}$ floor | $10^{-15}$ | see below |

Every other bound, $10^{-14}$ to $5 \times 10^{-9}$, is a safety bound;
observed errors are of order $10^{-14}$ or less. Recorded values cannot compare
exactly, because another compiler changes the last bits, by up to 27
representable values with GCC; the bound sits far above that and far below any
change of algorithm, ordering or branch selection. Running the suite before and
after a refactor with one compiler pins it bit for bit.

## Reproduce

```bash
ctest --test-dir build-tests -R '^component_regression_' --output-on-failure
```

See [Building and running](../../docs/developer/testing.md#building-and-running).
