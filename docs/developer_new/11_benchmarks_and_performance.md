# Benchmarks and performance

Benchmarks measure time; performance counters measure how much numerical work
a call did. Neither is evidence of correctness, which is the job of the
[test suite](10_testing_and_validation.md).

## The benchmark project

[`tests/benchmarks/`](../../tests/benchmarks) is a standalone CMake project. It
finds an installed QMOMKit with `find_package(QMOMKit CONFIG REQUIRED)`. Google
Benchmark is found with `find_package(benchmark)` or, if absent, fetched with
`FetchContent` from its GitHub repository at tag `v1.8.5`. It is a dependency
of this project only. The benchmarks also see the repository's `src/`, so they
can time internal kernels; they are not installed and are not part of the
package.

Build in `Release` (`RelWithDebInfo` for profiling), against a library
installed with the same configuration:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build-release --target install

cmake -S tests/benchmarks -B build-benchmarks -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build build-benchmarks

./build-benchmarks/cqmom_benchmark
./build-benchmarks/cqmom_benchmark --benchmark_filter=diagnostic_overhead
```

| Executable | Measures |
| --- | --- |
| `wheeler_benchmark` | fixed and adaptive Wheeler inversion |
| `vandermonde_benchmark` | the internal conditional Vandermonde solve, dense against Björck–Pereyra: matrix size, right-hand-side width, clustered, mixed-sign and tiny-weight nodes, and solver selection |
| `cqmom_benchmark` | CQMOM reconstruction: conditioning order, adaptive ragged branches, deep and wide high-dimensional hierarchies, diagnostic level, and the public, internal recursive and internal trusted-moment entry points as dimension and node count grow |
| `adaptive_quadrature_benchmark` | adaptive QMOM and CQMOM node counts: `Reconstruction` and `FullStep` timings for retained, grown, reduced and rejected counts, and `Direct`/`Fallback` pairs that reach the same count directly or by falling back from a larger maximum |
| `allocation_benchmark` | heap allocations and bytes per phase-space step and per reconstruction, through a replaced global `operator new` that counts only around the measured call; steps run with both velocity callback forms (`_packed_velocity` cases pass a `PhaseSpaceVelocity` directly) |
| `performance_instrumentation_benchmark` | representative reconstruction, moment-source, node-advection, splitting-scheme and Vandermonde workloads with the performance counters below, a grouped report, and `performance_instrumentation_overhead_*` cases that run with instrumentation compiled in but no active session |

The instrumentation benchmark reads its counters from the installed library:
build that library with `QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`. The
definition propagates through the installed target, so the benchmark project
needs no flag. Against an uninstrumented library its counter cases skip with a
message.

Interpretation:

- Timings are never correctness evidence. A benchmark may check that its
  fixture is valid and skip itself if not; correctness assertions belong in
  the tests.
- Compare only runs from the same machine, compiler, build type and load.
  Absolute timings do not transfer.
- To measure instrumentation overhead, install two libraries that differ only
  in the option, build one benchmark project against each, and run the same
  filter on both.
- Timings say how expensive a workload is; counters say why. Counters suggest
  hypotheses; confirm them with a profiler before optimizing.
- Each registered case answers a distinct performance question. Fixtures are
  deterministic and built outside the timed loop. Avoid Cartesian products of
  solvers, dimensions and node counts.

## Performance counters

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
| `ReconstructionCounters` | QMOM and CQMOM reconstructions, branch solves, Wheeler calls and adaptive reductions, terminal nodes, moment-layout builds, prepared-metadata builds and uses, dynamic-metadata builds |
| `LinearAlgebraCounters` | Vandermonde solves by solver, Björck–Pereyra transform updates, dense solves, residual checks, condition estimates |
| `PhaseSpaceEvolutionCounters` | macro steps, advection and moment-source operator calls, moment reconstructions, post-source inversions |
| `CqmomTopologyCounters` | active-count derivations, validations and transition builds |
| `TimeIntegrationCounters` | moment-source Runge–Kutta stages and retries, node right-hand-side evaluations, accepted and rejected embedded Runge–Kutta steps |

With the option off, the header defines the same names with no state: the
`record_*()` helpers are empty inline functions, `CounterSession` holds no
counters and returns zeros, and there is no `thread_local`. Call sites
therefore need no `#if`, and a normal build carries neither storage nor
branches in hot loops.

To add a counter:

1. Add the field to its group struct.
2. Add a `record_*()` helper in both branches of the header, the second one
   empty.
3. Call it where the event happens in `src/`.
4. Assert the expected count in
   [`tests/unit/performance_instrumentation.cpp`](../../tests/unit/performance_instrumentation.cpp)
   inside `#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)`; the
   instrumented CI-Debug build runs it.
5. Report it in `performance_instrumentation_benchmark` if it explains a
   workload.

Counters are not diagnostics. They are internal, not public API; they change no
result, status or output, and they are separate from `CqmomDiagnosticLevel` and
`PhaseSpaceEvolutionDiagnostics`. Diagnostics tell a caller what happened
numerically; counters tell a maintainer how much work it took.
