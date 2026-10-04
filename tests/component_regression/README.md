# Component Regression Tests

Component regression tests protect important implementation behavior that is
larger than an isolated unit contract but narrower than mathematical
verification. They are useful for diagnostics, status propagation, policy
mapping, direct component composition, and previously fragile behavior that
should remain visible during maintenance.

These tests are separate from mathematical verification. A component regression
may use manufactured inputs or simple analytical expectations, but its primary
purpose is regression protection rather than executable documentation of a
public mathematical workflow.

Component regression tests are not failed verification tests. They protect a
different category of behavior: implementation stability for important
components, overloads, diagnostics, and status propagation.

Component regression tests are also not a replacement for unit tests. Isolated
function, class, option, indexing, storage, and edge-case behavior belongs under
`tests/unit/`. When a regression case becomes a user-facing mathematical
workflow with an independent reference, promote it to `tests/verification/`
instead of growing this category.

When adding a test here, the primary question should be: "Would a future
maintainer accidentally break this important behavior even though the broader
mathematical verification still passes?"

Current component regression coverage includes:

- CQMOM direct compute failure/status mapping across `CqmomPreparedCase` and
  `CqmomCase`;
- direct moment-source update behavior and diagnostics;
- direct node-advection behavior and diagnostics;
- focused phase-space evolution integrator composition and phase-space ordering behavior;
- a historical CQMOM fixture regression for adaptive terminal reduction/ragged
  reinversion and active-count persistence through zero-rate state evolution.
  Its original generation provenance is not recoverable, so it is not
  independent numerical verification.
- CQMOM recorded reference values across a case matrix spanning phase-space
  dimension, conditioning order, node counts, adaptive branch counts, bounded
  coordinate support, every Vandermonde solver policy, and all three splitting
  schemes. Each static case is additionally checked against the quadrature
  nodes its moments were computed from, so its recorded values are analytically
  anchored rather than merely self-consistent. Evolution cases are checked
  against their recorded values and for agreement between the states created
  from the open and the prepared case;
- prepared-case moment replacement: replacing a prepared case's initial
  moments programmatically matches parsing the replacement from configuration,
  and replacing them with `CqmomMomentData` matches replacing them with a full
  tensor, in case settings, quadrature, diagnostics, state, and one evolution
  step;
- CQMOM reconstruction of two moment sets that are harder to invert than a rule
  whose nodes are far apart: one whose velocities have been changed by drag, and
  one with two velocity nodes 1.5e-6 apart.

## The CQMOM case matrix

`cqmom_case_matrix.*` defines the case matrix `cqmom_reference_values.cpp`
runs. That test asks *did the result change?*: status and node count must match
exactly, and every weight and abscissa must agree with its recorded value to
1e-12 relative, with a 1e-15 absolute floor. The budget absorbs the differing last bits another
compiler produces for the same computation, measured at up to 27 representable
values, while staying far below any change of algorithm or ordering. Pinning a
refactor more tightly than that needs a before-and-after run with one compiler,
where results are reproducible bit for bit.

The six Vandermonde solver cases exercise each solver's path through the
recursive reconstruction. The systems they solve are well conditioned, so the
solvers agree to within about 1.5e-15 relative, and some pairs bit for bit. The
1e-12 recorded comparison therefore cannot tell one solver from another: routing
one solver to another implementation would still pass here. Which solver is
selected and which algorithm runs is verified by the Vandermonde solver-policy
verification under `tests/verification/vandermonde/`.

### Regenerating the recorded values

`cqmom_reference_values.inc` is generated, but no generator is shipped: the
values are meant to be stable, and regenerating them should be a deliberate act
with the same justification as any other change to numerical behaviour. To
regenerate, write a short program that links the library, includes
`cqmom_case_matrix.hpp`, walks `static_cases()` and `evolution_cases()` exactly
as the test does, and prints each accepted quadrature at `%.17g`. Keep it
outside the repository.
