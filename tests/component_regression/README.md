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

- CQMOM direct compute failure/status mapping across public case
  representations;
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
  schemes. Each case is additionally checked against the quadrature nodes its
  moments were computed from, so the recorded values are analytically anchored
  rather than merely self-consistent;
- CQMOM reconstruction of two moment sets that are harder to invert than a rule
  whose nodes are far apart: one whose velocities have been changed by drag, and
  one with two velocity nodes 1.5e-6 apart.

## The CQMOM case matrix

`cqmom_case_matrix.*` defines the case matrix `cqmom_reference_values.cpp`
runs. That test asks *did the result change?* and compares against recorded
values to 1e-12 relative. The budget absorbs the differing last bits another
compiler produces for the same computation, measured at up to 27 representable
values, while staying far below any change of algorithm or ordering. Pinning a
refactor more tightly than that needs a before-and-after run with one compiler,
where results are reproducible bit for bit.

The six Vandermonde solver cases can produce identical numbers to one another,
because the systems they solve are well conditioned. They are still six
separate cases: each pins its own dispatch path, so silently routing one solver
to another implementation changes that case's recorded values and no others.

### Regenerating the recorded values

`cqmom_reference_values.inc` is generated, but no generator is shipped: the
values are meant to be stable, and regenerating them should be a deliberate act
with the same justification as any other change to numerical behaviour. To
regenerate, write a short program that links the library, includes
`cqmom_case_matrix.hpp`, walks `static_cases()` and `evolution_cases()` exactly
as the test does, and prints each accepted quadrature at `%.17g`. Keep it
outside the repository.
