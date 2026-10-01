# Benchmarks

Benchmarks measure computational performance. They do not establish
correctness (see [testing.md](testing.md)) or numerical credibility (see
[validation.md](validation.md)), and benchmark timings are never correctness
evidence.

QMOMKit is embedded in larger simulation codes where phase-space evolution runs
many times per timestep, so reconstruction, traversal, allocation, and
diagnostic costs are tracked deliberately.

## Executables

`tests/benchmarks/` is a standalone CMake project. Each source file builds one
executable:

| Executable | Measures |
|---|---|
| `wheeler_benchmark` | Fixed and adaptive one-dimensional Wheeler inversion. |
| `vandermonde_benchmark` | Vandermonde solver kernels: size and right-hand-side scaling, clustered/mixed-sign/tiny-weight robustness, solver selection, and the public one-level conditional reconstruction. |
| `cqmom_benchmark` | Representative full CQMOM workloads: conditioning-order traversal, adaptive ragged output, high-dimensional recursion, diagnostic overhead, and moment-representation scaling. |
| `adaptive_quadrature_benchmark` | Adaptive QMOM/CQMOM lifecycle: paired `Reconstruction` and `FullStep` timings for retained counts, rejected growth, reductions, and `Direct`/`Fallback` reconstruction-ceiling pairs. |
| `performance_instrumentation_benchmark` | Representative workloads reported with internal workload counters (see below). |
| `allocation_benchmark` | Heap allocations and bytes per phase-space step and per reconstruction, via a counting `operator new`. Phase-space steps are measured for both velocity callback forms, directly written and adapted per node. |

## Build and run

Benchmarks consume an **installed** QMOMKit package, so build and install the
library first, with matching configuration. Use `Release` for timing and
`RelWithDebInfo` when profiling; never `Debug` or an empty `CMAKE_BUILD_TYPE`.

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PWD}/install"
cmake --build build-release --target install

cmake -S tests/benchmarks -B build-benchmarks -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="${PWD}/install"
cmake --build build-benchmarks

./build-benchmarks/cqmom_benchmark
./build-benchmarks/cqmom_benchmark --benchmark_filter=diagnostic_overhead
```

`CMAKE_PREFIX_PATH` must be an absolute path to the prefix containing
`lib/cmake/QMOMKit/QMOMKitConfig.cmake`.

Google Benchmark is a dependency of this project only. It is found via
`find_package`, or fetched at a pinned version if absent. It is never a core
library, test, example, or Monte Carlo dependency.

### Performance instrumentation benchmarks

`performance_instrumentation_benchmark` reports internal counters only when the
installed library was built with
`-DQMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`; the compile definition
propagates from the installed target, so no benchmark-side flag is needed.
Without it the executable still builds and runs, and counter benchmarks skip
with an explanatory message.

To measure instrumentation overhead, install two packages with identical
compiler settings and different prefixes (`...=ON` and `...=OFF`), configure
one benchmark build against each, and run the same
`--benchmark_filter=performance_instrumentation_overhead` on both. Do not link
instrumented and non-instrumented libraries into one executable.

## Reading the output

All executables print the Google Benchmark timing table. Two add counters:

`allocation_benchmark` reports `allocations_per_op` and `bytes_per_op`.
Counting is enabled only around the measured call, so harness allocations are
excluded. Cases named `..._packed_velocity` pass a `PhaseSpaceVelocity` written
directly; the corresponding case without that suffix passes a
`PerNodePhaseSpaceVelocity` adapted through `make_phase_space_velocity`, so the
pair shows what the adapter costs.

`performance_instrumentation_benchmark` additionally prints a grouped report:
workload shape (dimensions, node counts, moment counts, schemes); timing per
reconstruction, source update, node update, or state-step; reconstruction
counters (Wheeler calls, branch solves, terminal nodes, layout builds,
prepared-versus-dynamic metadata); linear-algebra counters (Vandermonde solves
by backend, residual checks, condition estimates); phase-space
evolution counters (operator applications, post-source inversions, SSP-RK
stages, node RHS evaluations, embedded-RK steps); and derived ratios such as
Wheeler calls per reconstruction.

Where a benchmark has an analytical trajectory, it also reports maximum
trajectory error as an **absolute** error in the units of the phase-space
coordinate, taken as the largest component difference over all nodes,
dimensions, and iterations.

Interpret timing and counters together: timing says how expensive a workload
was, counters say why. High Wheeler calls per reconstruction point at recursive
branch solves; high node RHS evaluations per update point at the velocity
callback and node-advection scheme. Counters produce hypotheses; confirm with a
profiler before optimizing. Scheme comparisons (SSP-RK3 versus adaptive RK23)
are cost/accuracy tradeoff studies, not rankings.

Absolute timings and derived percentages are machine-specific local evidence,
not portable guarantees.

## Design principles

Fixtures must be deterministic, and inputs constructed outside the measured
loop when the goal is to time a kernel rather than fixture generation. Each
registered case should answer a distinct performance question; record why a
problem size is meaningful.

Avoid Cartesian products of solver modes, diagnostic levels, dimensions, and
node counts. Prefer one diagnostic-overhead case over every workload at every
level, and keep solver comparisons in `vandermonde_benchmark` rather than
application files. Exhaustive parameter grids are fine for a temporary local
investigation but should not become permanent registered cases.

Correctness assertions belong in unit and verification tests. A benchmark may
perform a setup sanity check and skip itself if the fixture is invalid.

Do not add broad benchmark result tables to developer documentation; keep
detailed results in benchmark tooling or external analysis.

## Public API boundary

Benchmark sources are development artifacts. They are not installed, exported,
or part of the public package. They may use internal headers where a specific
internal kernel must be measured; doing so does not make those headers public
API.
