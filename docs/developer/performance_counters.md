# Performance counters

Performance counters measure how much numerical work a call did. They are not
evidence of correctness, which is the job of the
[test suite](testing_and_validation.md), and they are not timings: to find
where time goes, profile a `Release` or `RelWithDebInfo` build of a
representative workload, and compare only runs from the same machine,
compiler, build type and load. Counters say why a workload is expensive and
suggest hypotheses; confirm them with a profiler before optimizing.

## The instrumentation option

The option `QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION` (default `OFF`, root
[`CMakeLists.txt`](../../CMakeLists.txt)) adds the compile definition of the
same name to `qmomkit` as a `PUBLIC` definition, so the library and the tests,
which include `src/internal/`, compile the same counter implementation. The
header is not installed, so only the tests read the counters.

[`src/internal/performance.hpp`](../../src/internal/performance.hpp), namespace
`qmomkit::internal::performance`, holds the whole layer:

| Name | Role |
| --- | --- |
| `PerformanceCounters` | one aggregate with a member per counter group, each a struct of `std::size_t` counters |
| `CounterSession` | owns a `PerformanceCounters`; `snapshot()` copies it, `reset()` zeroes it |
| `CounterScope` | RAII: makes a session the active one of the calling thread and restores the previous one on destruction, so scopes nest |
| `record_*()` | event helpers called from `src/`; each increments the active session's counter and does nothing when no session is active |
| `enabled()` | whether the counters are compiled in |

The active session is a `thread_local` pointer, so a scope records work done on
its own thread only.

Each counter is incremented once per event; tests assert exact counts.

| Member | Counter | Event |
| --- | --- | --- |
| `reconstruction` | `qmom_reconstructions` | QMOM reconstruction (`compute_qmom_from_moments`) |
| | `cqmom_reconstructions` | CQMOM reconstruction (`compute_recursive_quadrature_from_validated_moment_data`) |
| | `wheeler_calls` | Wheeler inversion at one coordinate scale (`compute_quadrature_at_scale`) |
| | `moment_layout_builds` | `create_cqmom_moment_layout` or `create_branch_moment_layout` |
| | `cqmom_prepared_metadata_builds` | prepared branch structure built |
| | `cqmom_prepared_metadata_uses`, `cqmom_dynamic_metadata_builds` | CQMOM reconstruction with prepared metadata, or building its own |
| `linear_algebra` | `vandermonde_solves` | conditional Vandermonde solve, any solver |
| | `bjorck_pereyra_solves`, `pivoted_qr_solves`, `svd_solves` | solve with that solver; `moore_penrose` counts as `svd`, and `dense` has no counter of its own |
| | `bjorck_pereyra_transform_updates` | right-hand-side update in the Björck–Pereyra elimination |
| | `dense_linear_solves` | `solve_dense`: a `dense` Vandermonde solve or a condition estimate |
| | `residual_checks`, `condition_estimates` | conditional-solve diagnostic computed |
| `phase_space_evolution` | `post_source_inversions` | post-source reconstruction |
| `cqmom_topology` | `active_count_derivations` | successful adaptive CQMOM reconstruction that returns its active counts |
| | `active_count_validations` | `validate_active_branch_counts` |
| | `active_count_transition_builds` | comparison of active counts around a source step, when either list is non-empty |
| `time_integration` | `source_rk_stages` | SSP-RK3 rate evaluation |
| | `source_rk_retries` | halving of a rejected source step |
| | `node_rhs_evaluations` | node per velocity evaluation in node advection |
| | `embedded_rk_accepted_steps`, `embedded_rk_rejected_steps` | accepted or rejected embedded Runge–Kutta trial step |

With the option off, the header defines the same names with no state: the
`record_*()` helpers are empty inline functions, `CounterSession` holds no
counters and returns zeros, and there is no `thread_local`. Call sites
therefore need no `#if`, and a normal build carries neither storage nor
branches in hot loops.

## Adding a counter

Every counter is asserted by a test; a counter that nothing reads is removed.
To add one:

1. Add the field to its group struct and a row to the table above.
2. Add a `record_*()` helper in both branches of the header, the second one
   empty.
3. Call it where the event happens in `src/`.
4. Assert the expected count in the unit test of the owning subsystem or in
   [`tests/unit/performance_instrumentation.cpp`](../../tests/unit/performance_instrumentation.cpp),
   guarded by `#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)` or
   `performance::enabled()`; the instrumented CI-Debug build runs it.

Counters are internal and change no result, status or output. Diagnostics
tell a caller what happened numerically; counters tell a maintainer how much
work it took.
