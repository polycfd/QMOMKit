# Performance Instrumentation

## Contents

- [Purpose and philosophy](#purpose-and-philosophy)
- [Relationship with diagnostics and profiling](#relationship-with-diagnostics-and-profiling)
- [Performance analysis workflow](#performance-analysis-workflow)
- [Instrumentation design](#instrumentation-design)
- [Available counter categories](#available-counter-categories)
- [Interpreting counters](#interpreting-counters)
- [Performance baselines](#performance-baselines)
- [Future extensions](#future-extensions)

## Purpose and philosophy

QMOMKit is designed to be embedded in larger simulation, validation, and
analysis codes. In large CFD-style applications, QMOM and CQMOM reconstruction,
phase-space evolution, moment-source integration, and linear solves may be
called many millions or billions of times. A small per-call cost can therefore
become a significant application-level cost.

Performance work in QMOMKit should follow a simple rule:

```text
Measure before optimizing.
```

Performance instrumentation exists to help maintainers answer workload
questions before changing algorithms or implementation structure:

- How much computational work was performed?
- Which numerical operations dominate a workload?
- How does the workload scale with dimensions, node counts, branch counts,
  moment counts, phase-space evolution steps, or solver choices?

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
as reconstruction calls, Wheeler calls, Vandermonde solver usage, phase-space
evolution operator applications, adaptive CQMOM active-count derivation,
validation, and transition construction, moment-source SSP-RK stages,
moment-source retries, node RHS evaluations, and adaptive embedded-RK
node-advection steps. These counters describe workload composition, not
correctness.

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
necessarily a bottleneck, and a locally expensive operation may not matter if
it is rarely called. Counters help identify plausible dominant operations; they
do not replace benchmark timing or profiling evidence.

The standalone
[`performance_instrumentation_benchmark`](benchmarks.md#performance-instrumentation-benchmarks)
provides initial representative QMOMKit workloads that combine Google Benchmark
timing with the internal counter layer. Build it against an installed QMOMKit
package configured with `QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`; the
benchmark project discovers that package through `find_package(QMOMKit)`. Its
default output includes a grouped human-readable report so maintainers can
connect timings to reconstruction, linear algebra, phase-space evolution, and
time-integration workload counters. The benchmark also includes comparison
variants for source-layout work, Vandermonde solver policy, phase-space
evolution components, and instrumentation overhead; those comparisons are
designed to test performance hypotheses before any optimization work begins.

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
`PhaseSpaceEvolutionDiagnostics`, or numerical status reporting.

## Available counter categories

Current counters are grouped by workload area.

Reconstruction counters summarize QMOM and CQMOM reconstruction work, including
QMOM/CQMOM reconstruction calls, Wheeler calls, recursive branch solves,
adaptive Wheeler reductions, terminal node assembly, and moment-layout
construction. Fixed normal CQMOM workflows also report prepared-metadata
builds, prepared-metadata uses, and dynamic-metadata builds.

Linear algebra counters summarize Vandermonde and dense solve work, including
Vandermonde solves, solver paths, BP transform updates,
residual checks, and condition estimates. BP transform updates count the
innermost-loop RHS updates in the two BP transforms; they are deterministic
kernel-work counts, not timings.

Phase-space evolution counters summarize operator-split work, including macro
steps, advection and moment-source operator applications, moment
reconstruction calls, and post-source inversions.

CQMOM topology counters summarize adaptive active-count work. A derivation is
recorded only after a complete active-count collection has been derived from an
accepted adaptive quadrature; fixed-topology early returns are not counted. A
validation is recorded when recursive reconstruction validates a non-empty
active-count constraint collection. A transition build is recorded when
non-empty before/after active counts are processed for transition diagnostics.

Time-integration counters summarize RK-related work, including moment-source
SSP-RK stages, moment-source retries, node RHS evaluations, and adaptive
embedded-RK accepted or rejected steps when RK23 or RK45 node advection is
selected.

These groups intentionally avoid per-branch timing, per-dimension breakdowns,
allocation tracking, callback timing, and other finer-grained details until
measurements justify them.

## Interpreting counters

Counters describe workload, not performance by themselves. Interpret them
together with benchmark timings, profiler output, and the numerical problem
being solved.

A high Wheeler-call count may indicate frequent reconstruction, many CQMOM
branches, many phase-space evolution macro steps, or a modeling choice that
requires repeated inversion. It does not by itself prove that Wheeler inversion
is the runtime bottleneck.

A high moment-layout-build count indicates metadata preparation work: required
CQMOM moment-layout generation, branch-local layout construction, or validation
paths that rebuild expected layouts. Interpret this with the benchmark type.
For public input paths, one expected-layout build may be intentional validation
work. For trusted `CqmomCase` or `CqmomState` paths, repeated global-layout
builds are usually more interesting than branch-local layouts.

`cqmom_prepared_metadata_builds` increments once when a fixed normal case or a
direct fixed normal state starts a new lineage and eagerly constructs its
immutable required-layout and branch-layout graph. It does not increment for
copies, transformations that preserve structure, or phase-space evolution
successors. `cqmom_prepared_metadata_uses` increments once when a fixed trusted
reconstruction consumes compatible metadata. `cqmom_dynamic_metadata_builds`
increments once per reconstruction that creates root/branch metadata
dynamically, such as an adaptive normal case or an explicit dense/lean case. A
valid fixed normal case/state never contributes a dynamic-metadata build;
missing metadata is an internal contract violation. These coarse counters do
not replace `moment_layout_builds`: the latter still counts individual
required-layout and branch-layout objects, including the objects created during
one eager prepared build.

A high node-RHS-evaluation count may reflect many terminal nodes, long
phase-space evolution intervals, difficult velocity dynamics, or tight
embedded-RK tolerances when RK23 or RK45 is selected. The next step is to
compare counters with node-advection timings, including the
`performance_node_advection_scheme` benchmark cases, and profile the velocity
callback path.

Use counters to form hypotheses. Use benchmarks and profiling to test them.

## Performance baselines

Counter values are deterministic integers and reproduce exactly across runs, so
they are the regression gate for optimization work. Wall-clock timings vary by
machine, build, and load; treat them as corroboration only and do not record
them here.

Baseline recorded at commit `7276a58`, Release build with
`QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON`.

Heap allocations per operation, from `allocation_benchmark`:

| Operation | allocations | bytes |
| --- | --- | --- |
| CQMOM 2D 2x2 Strang phase-space step | 564 | 23,528 |
| CQMOM 3D 2x2x2 Strang phase-space step | 927 | 56,928 |
| QMOM n=3 Strang phase-space step | 267 | 5,792 |
| CQMOM 2D 2x2 reconstruction | 203 | 8,392 |
| QMOM n=3 inversion | 13 | 920 |

### Current record

The baseline above is a fixed historical marker and is not updated. Record
current values here after work that changes them, so the two can be compared
directly. Values are the same whether performance instrumentation is compiled
in or out.

Recorded on the fill-style `OdeRhs` change, Release build, from
`allocation_benchmark`. Cite a commit here only once it is on `main`: a hash
taken from a branch stops resolving when that branch is squash-merged.

The phase-space steps are listed for both velocity callback forms: a
`PhaseSpaceVelocity` written directly, and a `PerNodePhaseSpaceVelocity`
adapted through `make_phase_space_velocity`, which costs one call and one
returned vector per node per right-hand-side evaluation.

| Operation | allocations, adapted | allocations, packed | bytes, adapted | bytes, packed |
| --- | --- | --- | --- | --- |
| CQMOM 2D 2x2 Strang phase-space step | 356 | 252 | 15,680 | 14,016 |
| CQMOM 3D 2x2x2 Strang phase-space step | 601 | 393 | 36,272 | 31,280 |
| QMOM n=3 Strang phase-space step | 168 | 90 | 3,760 | 3,136 |

| Operation | allocations | bytes |
| --- | --- | --- |
| CQMOM 2D 2x2 reconstruction | 153 | 7,152 |
| QMOM n=3 inversion | 11 | 776 |

Two changes to the Runge-Kutta layer produced these numbers. The controller
first began reusing its stage buffers across the attempts of an integration,
which took exactly 30 allocations off every phase-space step regardless of
dimension or callback form, because the number of attempts per step follows the
step size and the tolerances rather than the size of the node state. `OdeRhs`
then became a fill-style callback, which removed the per-stage derivative
vector: a CQMOM step sheds a further 20 allocations and a QMOM step a further
17. Reconstruction and inversion are unchanged, as they run no time integration.

The second change also made the cost structure qualitatively different. An
embedded-RK integration now allocates a fixed number of vectors — the state,
the result, four working buffers, and one per tableau stage — no matter how
many steps it takes. Integrating twenty RK23 steps of a 24-component state
costs 10 allocations where it once cost 162.

Workload counters per state-step are unchanged from the baseline, so the
reduction comes from removing redundant allocation rather than from doing less
numerical work.

Workload counters per state-step, CQMOM 2D 2x2 Strang
(`performance_cqmom_phase_space_evolution_strang_2d_2x2`):

| Counter | Per step |
| --- | --- |
| CQMOM reconstructions | 1 |
| Wheeler calls | 3 |
| CQMOM branch solves | 3 |
| Vandermonde solves | 1 |
| moment layout builds | 0 |
| prepared metadata uses | 1 |
| advection operators | 2 |
| moment-source operators | 1 |
| moment reconstructions | 3 |
| moment-source SSP-RK stages | 3 |
| node RHS evaluations | 32 |

Per-operation counters for fixed reconstruction, showing that the prepared and
dynamic metadata paths perform identical numerical work
(`wheeler_per_op` and `vandermonde_per_op` agree for each case) while differing
in layout construction:

| Case | Wheeler | Vandermonde | layout builds (dynamic) | layout builds (prepared) |
| --- | --- | --- | --- | --- |
| 2D 2x2 | 3 | 1 | 3 | 0 |
| 2D 3x3 | 4 | 1 | 3 | 0 |
| 3D 2x2x2 | 7 | 3 | 4 | 0 |
| 4D 2^4 | 15 | 7 | 5 | 0 |
| 5D 2^5 | 31 | 15 | 6 | 0 |

## Future extensions

Future instrumentation may add optional timing, allocation analysis,
finer-grained counter groups, or workload histograms if measurements show that
they would answer useful performance questions.

These are possible extensions, not committed implementation plans. New
instrumentation should remain separate from numerical diagnostics, preserve
zero overhead when disabled, and avoid exposing implementation details unless a
stable public API is deliberately designed.