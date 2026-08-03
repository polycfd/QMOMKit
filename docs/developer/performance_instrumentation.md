# Performance Instrumentation

## Purpose and philosophy

QMOMKit is designed to be embedded in larger simulation, validation, and
analysis codes. In large CFD-style applications, QMOM and CQMOM reconstruction,
moment transport, source integration, and linear solves may be called many
millions or billions of times. A small per-call cost can therefore become a
significant application-level cost.

Performance work in QMOMKit should follow a simple rule:

```text
Measure before optimizing.
```

Performance instrumentation exists to help maintainers answer workload
questions before changing algorithms or implementation structure:

- How much computational work was performed?
- Which numerical operations dominate a workload?
- How does the workload scale with dimensions, node counts, branch counts,
  moment counts, transport steps, or solver choices?

The instrumentation layer does not directly measure wall-clock performance. It
records counters that describe the shape of the computational workload. Timing,
cache behavior, vectorization, branch prediction, memory bandwidth, and
allocator behavior remain profiling concerns.

## Relationship with diagnostics and profiling

QMOMKit separates numerical diagnostics, performance instrumentation, and
external profiling because they answer different questions.

Numerical diagnostics explain what happened numerically. They report solver
choices, reconstruction status, failure context, branch diagnostics, residuals,
condition estimates, and related information needed to interpret the validity
and numerical path of a result.

Performance instrumentation counts computational work. It records events such
as reconstruction calls, Wheeler calls, Vandermonde solver usage, fallback
usage, transport operator applications, source SSP-RK stages, source retries,
node RHS evaluations, and adaptive embedded-RK node-advection steps. These counters
describe workload composition, not correctness.

Profiling tools measure elapsed time and hardware-level behavior. Use external
profilers to identify where time is spent, whether memory allocation or cache
behavior matters, and whether compiler optimization is effective for a
representative workload.

All three are useful, but they should not be combined into one reporting path.
Numerical diagnostics answer "what happened numerically?" Performance
instrumentation answers "how much work was performed?" Profilers answer "where
did time and hardware resources go?"

## Performance analysis workflow

Use performance instrumentation as part of a measurement-driven workflow:

1. Create a representative benchmark or application workload.
2. Collect performance counters for that workload.
3. Identify dominant operation counts and scaling behavior.
4. Use external profiling tools if the counters indicate a meaningful hotspot.
5. Optimize only measured bottlenecks.
6. Re-run the benchmark or workload to confirm the improvement.

Avoid speculative optimization. A frequently mentioned operation is not
necessarily a bottleneck, and a locally expensive operation may not matter if it
is rarely called. Counters help identify plausible dominant operations; they do
not replace benchmark timing or profiling evidence.

The standalone
[`performance_instrumentation_benchmark`](benchmarks.md#performance-instrumentation-benchmarks)
provides initial representative QMOMKit workloads that combine Google Benchmark
timing with the internal counter layer. Build it against an installed QMOMKit
package configured with `QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`; the
benchmark project discovers that package through `find_package(QMOMKit)`.
Its default output includes a grouped human-readable report so maintainers can
connect timings to reconstruction, linear algebra, transport, and
time-integration workload counters. The benchmark also includes comparison
variants for source-layout work, Vandermonde solver policy, transport
components, and instrumentation overhead; those comparisons are designed to
test performance hypotheses before any optimization work begins.

## Instrumentation design

Performance instrumentation is optional and compile-time controlled by
`QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION`. It is disabled by default.

Disabled builds compile instrumentation event calls to no-ops and allocate no
counter storage. This preserves normal library behavior and avoids adding
instrumentation branches or storage to hot loops in production builds that do
not request instrumentation.

Enabled builds record counters through internal event helpers when an internal
thread-local counter session is active. The event helpers are intent-based, so
call sites describe stable computational events rather than exposing detailed
implementation mechanics.

The first instrumentation layer records counters instead of timers. Counters
are less intrusive, easier to aggregate across repeated calls, and useful for
scaling analysis. Timing initially belongs in external profilers,
which can measure wall-clock time, CPU behavior, allocation cost, and hardware
effects without making timing policy part of the core numerical API.

The instrumentation layer is internal. It is not a public API, does not affect
output formats, and is not connected to `CqmomDiagnosticLevel`,
`TransportDiagnostics`, or numerical status reporting.

## Available counter categories

Current counters are grouped by workload area.

Reconstruction counters summarize QMOM and CQMOM reconstruction work, including
QMOM/CQMOM reconstruction calls, Wheeler calls, recursive branch solves,
adaptive Wheeler reductions, terminal node assembly, and moment-layout
construction.

Linear algebra counters summarize Vandermonde and dense solve work, including
Vandermonde solves, solver paths, dense fallback usage, residual checks, and
condition estimates.

Transport counters summarize operator-split transport work, including macro
steps, continuous and discontinuous operator applications, moment
reconstruction calls, and post-source inversions.

Time-integration counters summarize RK-related work, including source SSP-RK
stages, source retries, node RHS evaluations, and adaptive embedded-RK accepted
or rejected steps when RK23 or RK45 node advection is selected.

These groups intentionally avoid per-branch timing, per-dimension breakdowns,
allocation tracking, callback timing, and other finer-grained details until
measurements justify them.

## Interpreting counters

Counters describe workload, not performance by themselves. Interpret them
together with benchmark timings, profiler output, and the numerical problem
being solved.

A high Wheeler-call count may indicate frequent reconstruction, many CQMOM
branches, many transport macro steps, or a modeling choice that requires
repeated inversion. It does not by itself prove that Wheeler inversion is the
runtime bottleneck.

A high dense Vandermonde fallback count may indicate difficult conditioning,
an intentionally conservative solver policy, or a workload where fallback
robustness dominates the conditional reconstruction path. It may suggest an
optimization opportunity, but only profiling and benchmark comparisons can
show whether that opportunity is meaningful.

A high moment-layout-build count indicates metadata preparation work: required
CQMOM moment-layout generation, branch-local layout construction, or validation
paths that rebuild expected layouts. Interpret this with the benchmark type.
For public input paths, one expected-layout build may be intentional validation
work. For trusted `CqmomCase` or `CqmomState` paths, repeated global-layout
builds are usually more interesting than branch-local layouts.

A high node-RHS-evaluation count may reflect many terminal nodes, long
transport intervals, difficult velocity dynamics, or tight embedded-RK
tolerances when RK23 or RK45 is selected. The next step is to compare counters with
node-advection timings, including the
`performance_node_advection_scheme` benchmark cases, and profile the velocity
callback path.

Use counters to form hypotheses. Use benchmarks and profiling to test them.

## Performance Baselines

Performance baselines are regression references for future development. They
help maintainers answer where QMOMKit spends computational effort for
representative workloads, and they make counter changes visible when an
optimization or refactor claims to preserve numerical work.

Absolute timings depend on hardware, compiler, build type, dependency versions,
CPU scheduling, and thermal state. Workload counters are often more
transferable than timings because they describe the numerical work performed.
Interpret benchmark timings together with counters: timing says how expensive a
workload was in one environment, while counters explain the workload shape.

The following baseline was captured from
`performance_instrumentation_benchmark` on July 24, 2026 after the CQMOM source
and trusted reconstruction optimizations. Timings are normalized by the
operation named in the table and should not be treated as universal performance
claims.

Reproducibility information to record with future baselines:

- QMOMKit build type, normally `Release`.
- Compiler and version.
- Whether `QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`.
- Installed QMOMKit package used by the standalone benchmark project.
- Benchmark executable, filter, arguments, and Google Benchmark options.
- Hardware, operating system, and relevant runtime conditions.

Current local baseline environment:

- QMOMKit core: `Release`,
  `QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`.
- Benchmark project: `Release`, linked against the instrumentation-enabled
  installed QMOMKit package.
- Compiler: Apple clang 17.0.0, `/usr/bin/c++`, arm64 target.
- OS: macOS 26.5.2.
- Benchmark option used for the representative slices:
  `--benchmark_min_time=0.1s`.

### CQMOM Reconstruction

These rows use fixed-node public lean reconstruction workloads with automatic
Vandermonde selection and basic diagnostics. Counters are normalized per
successful reconstruction.

| Benchmark | Dimensions | Node Counts | Type | Time / Reconstruction | CQMOM Reconstructions | Wheeler Calls | Vandermonde Solves | Layout Builds | Source RK Stages | Node RHS Evaluations |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `performance_transport_component_inversion_2d_2x2` | 2 | `2 x 2` | CQMOM reconstruction | 8.953 us | 1 | 3 | 1 | 3 | 0 | 0 |
| `performance_vandermonde_auto_basic_2d_4x4` | 2 | `4 x 4` | CQMOM reconstruction | 28.670 us | 1 | 5 | 1 | 3 | 0 | 0 |
| `performance_vandermonde_auto_basic_3d_2x2x2` | 3 | `2 x 2 x 2` | CQMOM reconstruction | 22.903 us | 1 | 7 | 3 | 4 | 0 | 0 |

### CQMOM Transport

These rows use fixed-node `2 x 2` `CqmomState` transport workloads with 32
states and 4 macro steps per state. Counters are normalized per state-step.

| Benchmark | Dimensions | Node Counts | Type | Time / State-Step | CQMOM Reconstructions | Wheeler Calls | Vandermonde Solves | Layout Builds | Source RK Stages | Node RHS Evaluations |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `performance_cqmom_transport_strang_2d_2x2/32/4` | 2 | `2 x 2` | Strang CQMOM transport | 25.006 us | 2 | 6 | 2 | 5 | 3 | 32 |
| `performance_cqmom_transport_cd_2d_2x2/32/4` | 2 | `2 x 2` | CD CQMOM transport | 23.052 us | 2 | 6 | 2 | 5 | 3 | 16 |
| `performance_cqmom_transport_dc_2d_2x2/32/4` | 2 | `2 x 2` | DC CQMOM transport | 22.840 us | 2 | 6 | 2 | 5 | 3 | 16 |

### Source Integration

This row uses the fixed-node `3 x 3` CQMOM source integration workload with 64
states and 8 source updates per state. Counters are normalized per source
update.

| Benchmark | Dimensions | Node Counts | Type | Time / Source Update | CQMOM Reconstructions | Wheeler Calls | Vandermonde Solves | Layout Builds | Source SSP-RK Stages | Node RHS Evaluations |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `performance_source_normal_2d_3x3/64/8` | 2 | `3 x 3` | CQMOM source integration | 4.085 us | 0 | 0 | 0 | 1 | 3 | 0 |

The lower-level `performance_ssp_source_existing_layout_2d_3x3/64/8` benchmark
measured the same SSP-RK3 source-update kernel with existing layout metadata at
512.726 ns per source update and 3 source SSP-RK stages per update. It is useful
for separating callback/RK cost from public source-update validation and
metadata preparation.

### Interpretation Guidance

Use these baselines as workload fingerprints, not just timing rows:

- High Wheeler counts indicate recursive reconstruction complexity, usually
  from more branches or repeated transport inversions.
- Vandermonde solves indicate conditional reconstruction work; dense fallback
  counts, when nonzero, point to difficult conditional solves or conservative
  solver policy.
- Layout builds indicate metadata preparation work. Public paths may retain
  validation layout builds by design; trusted `CqmomCase` and `CqmomState`
  workflows should avoid repeated immutable global-layout validation.
- Source SSP-RK stages indicate discontinuous/source integration work. The normal
  SSP-RK3 source-update baseline performs 3 stages per update.
- Node RHS evaluations indicate continuous transport workload. Strang performs
  two continuous half-steps, which is why its node RHS count is higher than CD
  and DC for the same state and macro-step counts.

Recent optimization effects visible in this baseline:

- CQMOM source integration no longer performs redundant layout preparation in
  every SSP-RK3 stage; the public source-update path has one layout build per
  update, while the existing-layout SSP-RK3 kernel has none.
- Trusted `CqmomCase` and `CqmomState` reconstruction avoids repeated
  immutable metadata validation and global required-layout rebuilding. The
  `CqmomState` transport rows retain branch-local layout work but avoid
  extra global-layout validation work.
- Numerical workload counters such as Wheeler calls and Vandermonde solves are
  unchanged by these metadata optimizations.

## Future extensions

Future instrumentation may add optional timing, allocation analysis,
finer-grained counter groups, or workload histograms if measurements show that
they would answer useful performance questions.

These are possible extensions, not committed implementation plans. New
instrumentation should remain separate from numerical diagnostics, preserve
zero overhead when disabled, and avoid exposing implementation details unless a
stable public API is deliberately designed.
