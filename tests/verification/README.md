# Verification Framework

Verification tests are executable documentation for QMOMKit's mathematical and
workflow behavior. A verification case should make clear:

- what component or workflow is being exercised;
- why that behavior matters to the implemented design;
- how the expected result is obtained independently of the implementation under
  test;
- why the tolerance is appropriate for the numerical method and time-step or
  reconstruction scale.

A verification case prints nothing when it passes. Everything it would have to
say on success is already in its source header, and CTest hides passing output
anyway. When a case fails it writes the measured values, the bounds they missed,
and any refinement table needed to tell a near miss from a collapse, to standard
error, which `ctest --output-on-failure` surfaces. Files written by the output
workflow cases are the subject under test rather than diagnostics and are
produced on every run.

Verification differs from validation. Verification asks whether the
implementation produces the expected mathematical results for the documented
workflow. Validation asks whether the model is appropriate for external
physical, statistical, or reference data.

Verification is also separate from component regression. Component regression
tests protect implementation stability for important or fragile behavior;
verification tests demonstrate mathematical or public-workflow correctness
against an independent reference.

## Verification Suite Index

Use this table to find the canonical executable documentation for a public
workflow or mathematical component. Local directory READMEs describe the
individual cases, references, and tolerance choices in more detail.

| Area | Purpose | Main executable(s) | Reference type |
| --- | --- | --- | --- |
| QMOM | Verify public one-dimensional `compute_qmom(...)` success and representative direct-compute failures. | `verification_qmom_analytical_*`, `verification_qmom_compute_failure_*`, `verification_wheeler_analytical_*`, `verification_adaptive_qmom_*`, `verification_qmom_moment_space_boundary_*`, `verification_qmom_quadrature_bounds_*` | Manufactured discrete quadratures, classical Gaussian quadrature rules, public failure behavior |
| CQMOM | Verify the default `CqmomPreparedCase` and open `CqmomCase` workflows, adaptive branch counts, and high-dimensional cases. | `verification_cqmom_public_workflows`, `verification_cqmom_*` analytical cases | Analytical product quadratures, manufactured affine/branch-dependent quadratures |
| Recursive quadrature | Verify recursive CQMOM assembly conventions that are not exposed as standalone user workflows. | `verification_recursive_quadrature_manufactured_*`, `verification_cqmom_recursive_*` | Manufactured terminal quadratures and reconstruction identities |
| Phase-space evolution | Verify complete phase-space evolution workflows, stage-order behavior, public diagnostics, and failure propagation. | `verification_phase_space_evolution_*`<br>`verification_qmom_advection_matches_exponential_phase_evolution`<br>`verification_qmom_phase_space_source_evolution_matches_manufactured_solution`<br>`verification_qmom_node_source_matches_exact_translation` | Closed-form commuting/noncommuting evolution, manufactured phase-space evolution compositions |
| Source terms | Verify moment-source raw-moment source updates for QMOM and CQMOM workflows. | `verification_moment_source_*` | Closed-form source ODE solutions |
| Node advection | Verify advection quadrature-node advection for QMOM and CQMOM workflows. | `verification_node_advection_*` | Closed-form node trajectories |
| Realizability | Verify support-aware truncated-moment classification and diagnostics. | `verification_realizability_analytical_*` | Analytical Hamburger, Stieltjes, and Hausdorff moment examples |
| Distributions | Verify analytical distribution families, mixtures, independent products, and configuration-driven initialization, including configured mixture-marginal workflows. | `verification_distributions_analytical_*` | Closed-form distribution moments, mixture sums, independent-product tensors |
| Runge-Kutta | Verify standalone ODE utilities used by phase-space evolution and source components. | `verification_runge_kutta_analytical_*`<br>`verification_ssp_rk3_convergence`<br>`verification_rk23_fixed_step_convergence`<br>`verification_rk45_fixed_step_convergence`<br>`verification_rk23_adaptive_convergence`<br>`verification_rk45_adaptive_convergence` | Closed-form ODE solutions. Formal convergence order is measured only with the step size forced; the adaptive cases measure controller behaviour instead, since the controller chooses the steps and no step size exists to refine |
| Output | Verify public text output workflows for moments and quadratures. | `verification_output_writer_qmom`, `verification_output_writer_cqmom` | Computed public results compared with analytical QMOM/CQMOM references |
| Configuration | Verify representative user-facing configuration failure paths. | `verification_configuration_*` | Invalid `.qmom` fixtures with expected failure messages that name the configuration file and, where the fault sits on a line, a line number |
| Vandermonde solver policy | Verify user-facing CQMOM Vandermonde solver selection on ordinary and clustered conditional systems. | `verification_vandermonde_solver_*`<br>`verification_vandermonde_clustered_*` | Analytical product rule plus independently prescribed branch-dependent clustered CQMOM quadrature |

## Test Styles

Canonical verification cases are the preferred executable documentation. They
normally use named `.qmom` configurations or explicitly named analytical cases,
compare against closed-form or manufactured references, and report the measured
errors only on failure. New workflow coverage should be added in this style.

Component regression checks cover older or narrower component behavior that is
still useful but is not the primary documentation path for a feature. They live
under `tests/component_regression/` and use the `component_regression` CTest
label. When a behavior becomes an important user workflow, promote it to a
canonical named verification case instead of growing the regression test.

Internal manufactured checks may include implementation-only headers when they
verify internal numerical machinery directly. These should remain under
`tests/verification/` only when the mathematical reference is still independent
and the case documents a core numerical convention.

Unit tests are separate from verification. They belong under `tests/unit/` and
are the right place for exhaustive option validation, structural helper
behavior, status helpers, indexing utilities, and narrow failure-path variants.
Do not promote a unit test to verification unless it documents a public
mathematical workflow or an important reference case.

When deciding whether a new test belongs here, ask whether the expected result
comes from an analytical formula, manufactured solution, convergence study, or
independent workflow reference. If the main purpose is local edge-case coverage,
use `tests/unit/`. If the main purpose is preserving implementation behavior or
status propagation, use `tests/component_regression/`.

## Failure Reporting

Verification logic must not use C `assert`, because `assert` can be compiled
out in Release builds when `NDEBUG` is defined. Use
`QMOMKIT_VERIFY_CHECK(...)` from `verification_helpers.hpp` for simple
component checks, or explicit `VerificationResult`/`append_failure(...)`
patterns for larger analytical drivers.

Every verification failure must:

- run in Debug and Release builds;
- print enough context to identify the failed expectation;
- return a non-zero process exit code.

## Remaining Framework Gaps

The suite treats additional failure-path/status variants beyond the
representative phase-space evolution cases as unit or regression coverage
rather than canonical verification.
