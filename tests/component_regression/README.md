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
