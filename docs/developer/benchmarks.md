# Benchmarks

This page documents the QMOMKit benchmarking philosophy and benchmark project
architecture for maintainers, contributors, and AI-assisted development
workflows. It answers how QMOMKit measures and tracks computational performance.

For package boundaries and standalone consumer structure, see `architecture.md`.
For correctness testing, see `testing.md`. For numerical credibility and
reference comparisons, see `validation.md`.

## Role Of Benchmarks

QMOMKit separates correctness, credibility, and performance:

```text
testing:
    correctness

validation:
    numerical credibility

benchmarks:
    computational performance
```

Benchmarks measure the cost of representative numerical operations. They do not
establish mathematical correctness, validate physical behavior, define numerical
tolerances, or replace unit and verification coverage.

Performance tracking matters because QMOMKit provides numerical algorithms
intended to be embedded in larger simulation, analysis, and validation
workflows. Changes to
inversion, conditional reconstruction, CQMOM traversal, memory layout, or
diagnostic paths can affect downstream cost even when numerical results remain
unchanged. Benchmarks provide controlled workloads for comparing performance
across implementations, compilers, build settings, and hardware.

Benchmarks also guide Vandermonde solver decisions. Solver benchmarks measure
isolated numerical operations such as dense and Björck-Pereyra Vandermonde
solves. Application benchmarks measure representative CQMOM
workloads where recursive traversal, conditional reconstruction, adaptive
branching, storage, and diagnostics interact. Solver-selection policies should
be informed by both views: solver scaling shows where an algorithm can be
advantageous, and application benchmarks show whether that advantage matters in
realistic traversal patterns.

## Benchmark Architecture

Benchmarks live in:

```text
tests/benchmarks/
```

This directory is a standalone CMake project. It consumes an installed QMOMKit
package with:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(... PRIVATE qmomkit::qmomkit)
```

The benchmark project is separate from:

- the repository root core-library build;
- the normal `tests/` CTest correctness project;
- package-consumer smoke tests;
- examples and Monte Carlo validation tools.

Benchmark executables are not registered as CTest correctness tests. They are
run directly as performance measurement programs.

### Benchmark Build Configuration

Benchmark results require an optimized QMOMKit library and an optimized
benchmark executable. Build and install QMOMKit with `Release` for timing
comparisons:

```sh
cmake -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
cmake --install build-release
```

Then configure the standalone benchmark project with the same configuration:

```sh
cmake -S tests/benchmarks -B build-benchmarks \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=<install>
cmake --build build-benchmarks
```

For the default repository install prefix, use an absolute path. A bare
relative value such as `-DCMAKE_PREFIX_PATH=install` is resolved by CMake from
the benchmark build context and can fail to find the package configuration.
From the repository root, prefer:

```sh
cmake -S tests/benchmarks -B build-benchmarks \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="${PWD}/install"
cmake --build build-benchmarks
```

Run benchmark executables directly:

```sh
./build-benchmarks/wheeler_benchmark
./build-benchmarks/vandermonde_benchmark
./build-benchmarks/cqmom_benchmark
./build-benchmarks/performance_instrumentation_benchmark
```

Use Google Benchmark options to select a subset of registered benchmarks when a
focused comparison is more useful than running the full executable:

```sh
./build-benchmarks/vandermonde_benchmark --benchmark_filter=solver_selection
./build-benchmarks/cqmom_benchmark --benchmark_filter=diagnostic_overhead
./build-benchmarks/performance_instrumentation_benchmark --benchmark_filter=performance_cqmom_transport
./build-benchmarks/performance_instrumentation_benchmark --benchmark_filter=performance_node_advection_scheme
```

Counter-collecting performance-instrumentation benchmarks require QMOMKit to
be built and installed with
`QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`. The no-session overhead
variants can also be built against an instrumentation-disabled package for
timing comparison. Other benchmark executables do not require that option.

### Performance Instrumentation Benchmark Build

`performance_instrumentation_benchmark` uses QMOMKit's optional internal
performance instrumentation layer. Build and install the core library with
instrumentation enabled before configuring the standalone benchmark project:

```bash
INSTALL_PREFIX="${PWD}/install-performance"

cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DQMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON

cmake --build build
cmake --install build
```

Then configure and build the benchmark project against that installed package:

```bash
cmake -S tests/benchmarks -B build-benchmarks \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="${INSTALL_PREFIX}"

cmake --build build-benchmarks --target performance_instrumentation_benchmark
```

The benchmark project uses `find_package(QMOMKit CONFIG REQUIRED)`, so
`CMAKE_PREFIX_PATH` must name the installation prefix that contains
`lib/cmake/QMOMKit/QMOMKitConfig.cmake`. If a local shell variable is
inconvenient, pass the absolute path directly:

```bash
cmake -S tests/benchmarks -B build-benchmarks \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=/path/to/install-performance
```

Equivalently, point CMake at the package configuration directory:

```bash
cmake -S tests/benchmarks -B build-benchmarks \
    -DCMAKE_BUILD_TYPE=Release \
    -DQMOMKit_DIR=/path/to/install-performance/lib/cmake/QMOMKit
```

Run the full instrumentation benchmark executable directly:

```bash
./build-benchmarks/performance_instrumentation_benchmark
```

The default command prints the normal Google Benchmark timing table followed
by a human-readable QMOMKit performance instrumentation report. The report is
the preferred starting point for understanding why a representative workload is
expensive.

Use Google Benchmark filters for focused runs:

```bash
./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter=performance_cqmom_transport

./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter=performance_vandermonde

./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter=performance_source
```

If QMOMKit is built without
`QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`, the executable builds. The
enabled-session counter benchmarks skip with a message explaining the missing
configuration, while timing-only and no-session overhead variants still run.
This keeps the benchmark target discoverable while preserving the default
zero-overhead core-library build.

No additional benchmark-side instrumentation flag is required. The installed
`qmomkit::qmomkit` target propagates the
`QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION` compile definition when the core
library was built with that option enabled.

Use `RelWithDebInfo` instead when profiling benchmark workloads so optimization
remains enabled while debug symbols are retained:

```sh
cmake -B build-profile -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build-profile
cmake --install build-profile

cmake -S tests/benchmarks -B tests/benchmarks/build-profile \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_PREFIX_PATH="${PWD}/install"
cmake --build tests/benchmarks/build-profile
```

Do not use Debug or an empty `CMAKE_BUILD_TYPE` for benchmark timing results.
Debug builds keep assertions enabled and avoid aggressive optimization. Empty
single-configuration builds are generator-default builds and are not a stable
performance baseline. The installed QMOMKit package and benchmark executable
should be built with the same intended configuration before comparing results.

### Google Benchmark

The benchmark project uses Google Benchmark as an external dependency. This
dependency is local to `tests/benchmarks/`; it is not a QMOMKit runtime
dependency and is not required to configure, build, install, test, or consume
the core QMOMKit package.

The benchmark CMake project first looks for an installed Google Benchmark
package that exports `benchmark::benchmark_main`. If that target is not
available, the standalone benchmark project may fetch the pinned Google
Benchmark dependency for local benchmark builds. That fallback is scoped to the
benchmark project and should not leak into normal QMOMKit builds.

## Benchmark Organization

Current benchmark coverage is split by the performance question being asked.

### Application Benchmarks

Application-level CQMOM benchmarks live in:

```text
tests/benchmarks/cqmom_benchmark.cpp
```

These benchmarks measure representative QMOMKit CQMOM workloads. They are
intended to characterize full recursive reconstruction behavior, including
conditioning-order traversal, branch-local conditional reconstruction, adaptive
branch counts, terminal quadrature storage, diagnostics, and associated memory
movement.

Application benchmarks should represent meaningful workflows rather than solver
microbenchmarks. They normally use the automatic internal Vandermonde solver
selection path, because that is the application-facing behavior maintainers
need to track. Dense/Björck-Pereyra Vandermonde solver comparisons
belong in dedicated solver benchmarks unless a CQMOM benchmark is explicitly
studying an application-level solver-policy effect.

The current CQMOM application benchmark categories are:

- 2D conditioning-order baseline: a small fixed reconstruction with a
  non-default conditioning order, used to track public CQMOM overhead and
  dimension-order handling.
- Adaptive ragged 3D baseline: a branch-specific adaptive reconstruction that
  produces ragged terminal output, used to track adaptive traversal and
  terminal-path bookkeeping.
- Moderate fixed-dimensional reconstruction: a deterministic fixed CQMOM
  workload with moderate phase-space dimension and branch-local reconstruction
  work.
- High-dimensional low-node-count recursion: a deep tree with small node
  counts in each dimension, used to isolate recursive traversal and branch
  proliferation rather than large Vandermonde solves.
- Wide conditional reconstruction: a high-dimensional workload with a wider
  root conditional reconstruction, used to track the cost of broad
  branch-local right-hand sides inside full CQMOM traversal.
- Diagnostic overhead: one representative fixed workload run at each CQMOM
  diagnostic level, used to measure the incremental cost of diagnostic
  collection without multiplying every application case.
- Dense-vs-lean representation scaling: explicit dense `RawMomentTensor` and
  lean `CqmomMomentData` CQMOM workflows over selected dimensions and node
  counts. These benchmarks report runtime plus dense and lean moment counts so
  maintainers can quantify the memory/runtime tradeoff without introducing
  automatic representation selection.

### Performance Instrumentation Benchmarks

Instrumentation-based workload benchmarks live in:

```text
tests/benchmarks/performance_instrumentation_benchmark.cpp
```

The `performance_instrumentation_benchmark` executable provides representative
QMOMKit workload measurements using the optional performance instrumentation
layer. It reports Google Benchmark wall-clock/CPU timing together with
internal workload counters, such as reconstruction calls, branch solves,
Wheeler calls, Vandermonde solves, fallback counts, terminal node assembly,
transport operator calls, source SSP-RK stages, and node RHS evaluations.

This benchmark is intended for:

- performance analysis;
- workload characterization;
- optimization studies.

It is not:

- a correctness test;
- a verification case;
- a replacement for external profiling tools.

Use this benchmark to understand how much computational work a representative
workflow performs. Use profiling tools when the benchmark and counters suggest
that elapsed time, allocation behavior, cache effects, or callback cost need
closer inspection. Do not treat a counter value by itself as an optimization
result.

The current instrumentation benchmark categories are:

- CQMOM transport workload: constructs one `CqmomCase`, creates many
  independent `CqmomState` objects, and repeatedly advances them through
  `advance_transport_cqmom(...)` with deterministic velocity, source, and
  boundary callbacks. Strang, continuous-discontinuous, and
  discontinuous-continuous splitting variants make the expected operator
  sequence visible in normalized counters.
- CQMOM reconstruction and Vandermonde policy comparisons: isolate
  `compute_cqmom(...)` for 2D 2 x 2, 2D 4 x 4, and 3D 2 x 2 x 2 fixed-node
  cases with automatic, explicit Björck-Pereyra, explicit dense, and detailed
  automatic policies. These variants make residual checks, condition
  estimates, dense solves, and fallback paths visible at the application
  reconstruction level.
- CQMOM source integration and layout comparisons: repeatedly call
  `advance_cqmom_moments_with_source(...)` with a deterministic source model
  to expose source SSP-RK stage counts, retry behavior, and layout construction.
  Companion variants measure source callback evaluation alone,
  `CqmomMomentData` copying, and SSP-RK3 source integration around an already
  available layout.
- CQMOM transport component comparisons: isolate continuous node advection,
  discontinuous/source update, raw-moment reconstruction from an existing
  quadrature, and CQMOM inversion for a common 2D 2 x 2 workload. Use these
  with the complete transport variants to check whether integrated transport
  cost is consistent with the isolated kernels.
- Node-advection scheme comparisons: call `advance_nodes(...)` directly for
  small 1D QMOM and 2D 2 x 2 CQMOM quadratures. Constant-velocity cases show
  the fixed-stage work of SSP-RK3 and the adaptive bookkeeping of RK23 on a
  simple trajectory. Affine-velocity cases add an analytical non-trivial
  trajectory so the benchmark can report both workload counters and maximum
  trajectory error. These variants compare `NodeAdvectionScheme::ssp_rk3` with
  `NodeAdvectionScheme::adaptive_rk23` without changing transport splitting,
  source integration, or node-advection methods.
- Instrumentation overhead comparisons: run equivalent no-session workloads
  against instrumentation-enabled and instrumentation-disabled QMOMKit builds.
  Enabled-session variants then show the additional cost of active counter
  accumulation for the same workload family.

Useful focused filters are:

```bash
./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter='performance_source|performance_ssp_source|performance_cqmom_moment_data'

./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter='performance_vandermonde_(auto_basic|bp_basic|dense_basic|auto_detailed)'

./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter='performance_transport_component|performance_cqmom_transport_(strang|cd|dc)'

./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter=performance_node_advection_scheme

./build-benchmarks/performance_instrumentation_benchmark \
    --benchmark_filter=performance_instrumentation_overhead
```

### Instrumentation Overhead Workflow

Do not link instrumented and non-instrumented QMOMKit libraries into one
benchmark executable. Build two installed packages with identical compiler
settings and different install prefixes:

```bash
cmake -S . -B build-release-no-instrumentation \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PWD}/install-release-no-instrumentation" \
    -DQMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=OFF

cmake --build build-release-no-instrumentation
cmake --install build-release-no-instrumentation

cmake -S . -B build-release-instrumented \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PWD}/install-release-instrumented" \
    -DQMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON

cmake --build build-release-instrumented
cmake --install build-release-instrumented
```

Then configure two benchmark executables against those installations:

```bash
cmake -S tests/benchmarks -B build-benchmarks-no-instrumentation \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="${PWD}/install-release-no-instrumentation"

cmake --build build-benchmarks-no-instrumentation --target performance_instrumentation_benchmark

cmake -S tests/benchmarks -B build-benchmarks-instrumented \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="${PWD}/install-release-instrumented"

cmake --build build-benchmarks-instrumented --target performance_instrumentation_benchmark
```

Run the same no-session filter against both executables:

```bash
./build-benchmarks-no-instrumentation/performance_instrumentation_benchmark \
    --benchmark_filter=performance_instrumentation_overhead

./build-benchmarks-instrumented/performance_instrumentation_benchmark \
    --benchmark_filter=performance_instrumentation_overhead
```

The non-instrumented executable measures the compile-time-disabled path. The
instrumented no-session executable measures the enabled-build path without an
active counter session. To measure active counter accumulation, compare with
the corresponding enabled-session families, such as:

```bash
./build-benchmarks-instrumented/performance_instrumentation_benchmark \
    --benchmark_filter='performance_cqmom_transport_strang|performance_vandermonde_auto_basic_2d_2x2|performance_source_normal'
```

Compare runtime per reconstruction, runtime per source update, and runtime per
transport state-step. These comparisons quantify instrumentation cost for the
selected workload only; they do not prove anything about unrelated application
workloads.

### Reading The Instrumentation Report

`performance_instrumentation_benchmark` emits two complementary views:

- the Google Benchmark table reports elapsed time, CPU time, iteration count,
  throughput, and raw counter columns;
- the QMOMKit performance instrumentation report groups the same workload
  counters by numerical subsystem and adds derived metrics.

The report sections are:

- Workload: benchmark name, workload type, phase-space dimensions, node counts,
  requested terminal nodes, moment counts, transport scheme, node-advection
  scheme, state count, and step count when applicable.
- Timing: total measured benchmark-body runtime and average runtime per
  reconstruction, source update, node update, or state-step. This answers how
  expensive the measured workload was.
- Reconstruction counters: CQMOM/QMOM reconstruction counts, Wheeler calls,
  branch solves, terminal node assembly, and moment-layout construction.
- Linear algebra counters: Vandermonde solves, solver-path counts, dense
  fallback counts, residual checks, and condition estimates.
- Transport/time-integration counters: continuous/discontinuous operator
  applications, moment reconstruction calls, post-source inversions, source
  SSP-RK stages and retries, node RHS evaluations, and adaptive RK23 step
  counts when adaptive RK23 node advection is selected.
- Node-advection diagnostics: scheme-neutral node-advection steps, rejected
  steps, RHS evaluations, and maximum trajectory error for benchmark cases with
  an analytical trajectory.
- Trajectory error metric: the reported maximum trajectory error is an
  absolute error. For each benchmark update, the benchmark compares every
  advanced quadrature coordinate with the analytical trajectory value for that
  node and phase-space dimension, takes the largest absolute component
  difference over all nodes and dimensions, and reports the maximum seen across
  benchmark iterations. Values such as `1e-5` or `1e-9` therefore have the same
  units as the phase-space coordinate being advanced; they are not relative
  errors.
- Derived metrics: normalized ratios such as Wheeler calls per reconstruction,
  Vandermonde solves per reconstruction, dense fallback percentage, and
  transport operations per state update.

Interpret timing and counters together. Timing tells how expensive a workload
was; counters help explain why.

Example interpretations:

- High Wheeler calls per reconstruction usually indicate many recursive CQMOM
  branch solves. Compare this with dimensions, node counts, and terminal node
  assembly before optimizing Wheeler itself.
- High Vandermonde solves per reconstruction indicate conditional
  reconstruction is prominent in the workload. Dense fallback percentage helps
  distinguish ordinary branch traversal from difficult conditioning.
- High node RHS evaluations per state update indicate continuous node
  advection work. Profile the velocity callback and node-advection scheme
  settings before changing reconstruction code.
- Node-advection scheme-comparison benchmarks should be interpreted as a
  cost/accuracy tradeoff study, not a ranking. SSP-RK3 has predictable cost
  because it evaluates a fixed three-stage trajectory step over the supplied
  timestep. Adaptive RK23 may perform substantially more RHS evaluations
  because it estimates error and adjusts internal substeps. The resulting
  trajectory error depends on the velocity field, tolerance settings, and
  supplied transport timestep, so the useful scheme choice is application
  dependent.
- High source SSP-RK stages with low reconstruction and Vandermonde counts
  indicate source-update cost, callback cost, or CQMOM moment-context
  construction may dominate that benchmark.

Counters identify hypotheses. Use external profiling tools to confirm where
wall-clock time is spent before making optimization changes.

### Solver Benchmarks

Solver-level benchmarks live in:

```text
tests/benchmarks/vandermonde_benchmark.cpp
```

These benchmarks compare numerical algorithms, measure scaling behavior, and
evaluate dense and Björck-Pereyra Vandermonde solvers. They cover
matrix size scaling, right-hand-side width scaling, representative numerical
robustness cases, wrapper behavior, automatic selection, and public
conditional-reconstruction kernel cost.

Solver benchmarks are the right place for dense/Björck-Pereyra comparisons,
Vandermonde-specific solver stress cases, and numerical conditioning studies.

Other benchmark categories remain focused on their own numerical operations:

- Wheeler inversion: fixed one-dimensional inversion and adaptive Wheeler
  reduction behavior.
- Vandermonde and reconstruction operations: conditional moment reconstruction
  through the weighted Vandermonde solve.

These categories are intended to measure representative operations, not to
publish universal timing claims. Individual benchmark results depend on the
workload, compiler, build type, machine, dependency source, and runtime
environment.

## Vandermonde Solver Benchmarks

The Vandermonde benchmark suite evaluates Vandermonde solver choices
independently from CQMOM application workloads. These benchmarks
measure solver kernels and the public one-level reconstruction kernel without
the recursive CQMOM tree traversal, adaptive terminal storage, or CQMOM
diagnostic machinery measured by `cqmom_benchmark.cpp`.

### Solver Scaling

Solver-scaling benchmarks compare dense and Björck-Pereyra behavior as the
Vandermonde matrix size changes. They also compare behavior as the
number of right-hand sides changes for a fixed matrix size. These cases use a
small representative set of sizes and report counters such as `N`, `RHS`,
`N3`, and `N2_RHS` so results can be interpreted against the expected scaling
terms.

### Numerical Robustness

Numerical-robustness benchmarks characterize solver behavior for difficult but
meaningful systems. The permanent suite keeps representative cases for:

- clustered nodes;
- mixed-sign nodes;
- tiny weights.

These cases are not exhaustive numerical verification. Correctness,
near-singular rejection, and edge-case status behavior belong in unit or
verification tests. The benchmark role is to track the performance cost of
stable representative difficult systems.

### Solver Selection

Solver-selection benchmarks exercise the Vandermonde solver wrapper with
`dense`, `bjorck_pereyra`, `bjorck_pereyra_with_dense_fallback`, and
`automatic` modes. They report counters for the selected solver, algorithm
actually used, and dense fallback occurrence. These counters make policy
behavior visible in benchmark output and help verify the public CQMOM
configuration setting.

### Public Reconstruction Kernel

The public reconstruction-kernel benchmark measures the end-to-end cost of
`reconstruct_conditional_moments(...)` for a deterministic weighted
Vandermonde-backed conditional reconstruction. This is still a solver-level
benchmark: it includes public request/result overhead for one conditional
reconstruction level, but not full recursive CQMOM traversal.

## Benchmark Design Principles

### Reproducibility

Benchmarks should use deterministic fixtures. Inputs should be constructed
before the measured loop whenever the goal is to time a computational kernel rather
than fixture generation.

A benchmark fixture should define:

- the operation being measured;
- the moment sequence, tensor, hierarchy, or request;
- the node counts and adaptive settings;
- the problem size and dimensionality;
- whether setup and validation checks are outside the measured loop;
- the public or internal boundary intentionally exercised.

Each benchmark should answer a clear performance question. Representative
workloads are more useful than overly tiny smoke paths or artificially large
cases with unclear interpretation. When adding a larger case, record why the
problem size is meaningful.

Avoid Cartesian multiplication of solver modes, diagnostic modes, dimensions,
and node counts in application benchmark files. Add a matrix only when every
combination answers a distinct question. Prefer a single diagnostic-overhead
benchmark over running every CQMOM workload at every diagnostic level. Prefer
dedicated solver benchmarks for algorithm comparisons and solver stress cases.

CQMOM application benchmarks should stay interpretable: choose phase-space
dimension, node counts, adaptive settings, and moment construction so the
resulting workload maps to recursive traversal, conditional reconstruction,
ragged output, diagnostic overhead, or another recognizable CQMOM workflow.

Solver benchmark names should describe the question being asked, not only the
parameter values being swept. Prefer names such as matrix-size scaling,
right-hand-side scaling, clustered-node robustness, solver selection, or public
conditional reconstruction. Exhaustive parameter grids are appropriate for
temporary local investigations, but they should not become the permanent
benchmark suite unless each registered case has a distinct long-term
interpretation.

### Interpretation

Benchmark results are timing measurements for a specific workload. They are not
portable performance guarantees.

Interpret results in the context of:

- algorithm choice;
- phase-space dimension and moment count;
- node-count hierarchy and adaptive behavior;
- memory allocation and data movement;
- compiler, build type, and optimization flags;
- hardware, CPU scheduling, and thermal state;
- Google Benchmark dependency source and version.

Use benchmarks to compare controlled alternatives or track performance trends.
Do not treat one benchmark number as a universal metric for QMOMKit performance.

## Adding Benchmarks

When adding a benchmark, identify the numerical operation being measured before
choosing a fixture. The benchmark should make the cost boundary obvious: for
example, Wheeler inversion, conditional moment reconstruction, full recursive
CQMOM, transport composition, or output formatting.

Prefer measuring the kernel of interest. If setup is not part of the intended
measurement, construct inputs outside the measured loop and verify that the
fixture is numerically valid before timing begins.

Use representative deterministic inputs:

- raw moments from a known quadrature or analytical construction;
- phase-space hierarchies with meaningful node counts;
- adaptive settings that exercise the intended branch behavior;
- tensor extents that match the operation being measured.

Benchmark code should remain independent from correctness tests. It may reuse
stable public APIs and simple fixture construction patterns, but correctness
assertions belong in unit or verification tests. A benchmark may perform a setup
sanity check and skip itself if the fixture is invalid; that is not a substitute
for correctness coverage.

Do not add broad benchmark result tables to developer documentation. Store and
compare detailed results through the benchmark tooling or external analysis
appropriate to the performance investigation. Compact regression baselines for
the performance instrumentation suite may live in
[`performance_instrumentation.md`](performance_instrumentation.md#performance-baselines)
when they document representative workload counters and clearly identify the
environment behind any timing numbers.

## Relationship With Profiling

Benchmarks and profiling answer different questions:

- Benchmarks provide controlled performance comparisons for selected workloads.
- Profiling identifies where time is spent in a running program or workload.

Profiling can explain a benchmark result or guide optimization work, but
profiling output is not itself a stable developer-documentation artifact. Use
profiling after a benchmark identifies a workload worth investigating.

## Public/API Boundary

Benchmark sources are repository development artifacts. They are not installed,
exported, or part of the QMOMKit public package.

Benchmarks may use public QMOMKit APIs to represent downstream-visible costs.
They may also use repository-local fixture helpers or internal access if a
specific internal kernel must be measured, but doing so does not make those
helpers or internal headers public API.

Downstream users should consume the installed QMOMKit package and public headers.
Benchmark code should not be treated as a supported user-facing interface.
