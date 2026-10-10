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
same name to `qmomkit` as a `PUBLIC` definition, so every target that links
QMOMKit, installed or not, compiles the same counter implementation.

[`src/internal/performance.hpp`](../../src/internal/performance.hpp), namespace
`qmomkit::internal::performance`, holds the whole layer:

| Name | Role |
| --- | --- |
| `PerformanceCounters` | one struct per counter group, all `std::size_t` |
| `CounterSession` | owns a `PerformanceCounters`; `snapshot()` copies it, `reset()` zeroes it |
| `CounterScope` | RAII: makes a session the active one of the calling thread and restores the previous one on destruction, so scopes nest |
| `record_*()` | event helpers called from `src/`; each increments the active session's counter and does nothing when no session is active |
| `enabled()` | whether the counters are compiled in |

The active session is a `thread_local` pointer, so a scope records work done on
its own thread only.

| Group | Counts |
| --- | --- |
| `ReconstructionCounters` | QMOM and CQMOM reconstructions, Wheeler calls, moment-layout builds, prepared-metadata builds and uses, dynamic-metadata builds |
| `LinearAlgebraCounters` | Vandermonde solves by solver, Björck–Pereyra transform updates, dense solves, residual checks, condition estimates |
| `PhaseSpaceEvolutionCounters` | post-source inversions |
| `CqmomTopologyCounters` | active-count derivations, validations and transition builds |
| `TimeIntegrationCounters` | moment-source Runge–Kutta stages and retries, node right-hand-side evaluations, accepted and rejected embedded Runge–Kutta steps |

With the option off, the header defines the same names with no state: the
`record_*()` helpers are empty inline functions, `CounterSession` holds no
counters and returns zeros, and there is no `thread_local`. Call sites
therefore need no `#if`, and a normal build carries neither storage nor
branches in hot loops.

## Adding a counter

Every counter is asserted by a test; a counter that nothing reads is removed.
To add one:

1. Add the field to its group struct.
2. Add a `record_*()` helper in both branches of the header, the second one
   empty.
3. Call it where the event happens in `src/`.
4. Assert the expected count in
   [`tests/unit/performance_instrumentation.cpp`](../../tests/unit/performance_instrumentation.cpp)
   inside `#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)`; the
   instrumented CI-Debug build runs it.

Counters are internal and change no result, status or output. Diagnostics
tell a caller what happened numerically; counters tell a maintainer how much
work it took.
