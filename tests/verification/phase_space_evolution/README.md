# Phase-Space Evolution Analytical Verification

These cases verify the public phase-space evolution facility as a complete path:
configuration loading, QMOM/CQMOM inversion, advection operator calls,
raw-moment reconstruction, moment-source operator calls, repeated moment/quadrature
representation changes, diagnostics, and phase-space evolution scheme behavior for
`strang`, `advsrc`, and `srcadv`. They are automated analytical tests that complement
the user-facing phase-space evolution examples under `examples/qmom_1d_breakage_advection/`
and `examples/cqmom_2d_minimal/`.

All initial QMOMKit cases are loaded from `.qmom` files through the public
configuration API. Configuration may define the initial numerical setup and the
operator-splitting scheme (`strang`, `advsrc`, or `srcadv`). Operator
coefficients, exact references, timestep plans, RK23/SSP-RK3 options, and
tolerances live in `phase_space_evolution_analytical.cpp` because the current configuration
format does not define physics callbacks, model coefficients, timestep plans, or
detailed integration controls.

## Commuting QMOM

Configuration: `commuting_qmom.qmom`

The one-dimensional advection operator is:

```text
dxi/dt = a xi,  a = 0.2
```

so each raw moment satisfies `dm_k/dt = k a m_k` during the advection part.
The moment-source operator is diagonal:

```text
dm_k/dt = -lambda_k m_k
lambda_k = 0.05 k
lambda_0 = 0
```

Both moment-space actions are diagonal, so they commute. The exact combined
solution is:

```text
m_k(t) = exp((k a - lambda_k) t) m_k(0)
```

The test runs one complete step, a multi-step run to `T = 0.24`, and a
refinement sequence with `dt = 0.12, 0.06, 0.03`. The refinement isolates
ordinary RK and moment-source integration error because the phase-space evolution composition itself
has no commuting-operator error. The residual refinement behavior is expected
to be approximately third order because the remaining visible error is
dominated by the SSP-RK3 moment-source update rather than by phase-space evolution composition.

## QMOM Phase-Space Evolution Failure Propagation

Configuration: `qmom_phase_space_evolution_failures.qmom`

This case checks how QMOM phase-space evolution reports velocity-callback failures and
source-step failures. The configured moment set is realizable; the advection
failures are caused only by user callback mistakes:

- returning a velocity vector with the wrong dimension;
- returning a non-finite velocity value.

Both failures should stop during the first advection half-step. The verifier
checks that `PhaseSpaceEvolutionStatus::advection_failed`,
`PhaseSpaceEvolutionStage::first_advection_half_step`, and
`NodeAdvectionStatus::velocity_failed` are reported, that the velocity-callback
message is preserved in phase-space evolution diagnostics, and that no source or later
reconstruction stage is reported as completed. The returned moments are checked
against the input moments to document that a failed advection substep is not a
valid final phase-space evolution state.

The same case also verifies adaptive embedded-RK minimum-step propagation. It uses a
finite, dimensionally valid, deliberately stiff velocity field and RK23 options
with `initial_step = max_step = 0.04` and `min_step = 0.02`. The first Strang
advection half-step rejects the requested step, then advances finite updates
at the effective minimum step instead of failing solely because the controller
would prefer a smaller value. The verifier checks successful phase-space evolution,
successful embedded-RK diagnostics, rejected-step counts, and exact endpoint arrival.

The QMOM case also verifies post-source reconstruction failure. A finite,
correctly sized source callback drives `m0` negative during the moment-source
stage. The moment-source update succeeds, but the evolved moments are no longer
realizable as a QMOM quadrature. The verifier checks
`PhaseSpaceEvolutionStatus::inversion_failed`, `PhaseSpaceEvolutionStage::post_source_inversion`,
the non-success post-source QMOM status, the successful moment-source diagnostics, a
realizability-oriented message, and that the second advection half-step is not
reported as completed.

Configuration: `qmom_phase_space_evolution_inversion_failures.qmom`

This companion case verifies initial inversion failure propagation through the
same public QMOM phase-space evolution workflow. The input moments are finite and loaded
successfully, but `m0 < 0`, so the moment sequence is not realizable before
phase-space evolution can begin. The verifier checks that phase-space evolution reports
`PhaseSpaceEvolutionStatus::inversion_failed`, identifies
`PhaseSpaceEvolutionStage::initial_inversion`, preserves the QMOM realizability message,
does not advance time, and does not report source or reconstruction stages as
completed.

## Noncommuting QMOM

Configuration: `noncommuting_qmom.qmom`

The advection operator is constant translation:

```text
dxi/dt = b,  b = 0.45
```

For moments `u = [m0, m1, m2, m3]^T`, this gives the lower-bidiagonal action:

```text
dm0/dt = 0
dm1/dt = b m0
dm2/dt = 2 b m1
dm3/dt = 3 b m2
```

The moment-source operator is:

```text
dm_k/dt = -gamma k m_k,  gamma = 0.17
```

With matrices `A` for translation and `B = diag(0, -gamma, -2 gamma,
-3 gamma)`, the commutator is nonzero because the diagonal source rates differ
across the subdiagonal couplings in `A`.

The exact reference is the closed-form exponential of the lower-bidiagonal
matrix `A + B`. The verifier represents each component as a finite sum of
exponentials:

```text
m_k(t) = sum_{r=0}^k c_{k,r} exp(-gamma r t)
```

and computes the coefficients recursively from the subdiagonal entries
`k b`. This reference does not call QMOMKit RK solvers, inversion,
phase-space evolution, or an external numerical library.

The convergence study uses `T = 1.2` and:

```text
dt = 0.3, 0.15, 0.075, 0.0375, 0.01875
```

The error norm is the scaled infinity norm over all configured moments:

```text
max_k |computed_k - exact_k| / max(1, |exact_k|)
```

The finest refinement pairs are required to approach second order. RK23 is run
with tight tolerances and small maximum step, and the moment-source operator is a
third-order SSP-RK3 update, so the tested range is chosen to keep the
second-order composition error visible. The observed order should approach
two on the finest refinement pairs.

The same noncommuting manufactured case also verifies `advsrc` and `srcadv` as
first-order Lie splitting schemes. Those schemes are expected to approach first
order under refinement and to differ from each other at finite timestep,
documenting that the public scheme option is active rather than ignored.

## Noncommuting CQMOM 2D

Configuration: `noncommuting_cqmom_2d.qmom`

This fixed (2 x 2) public CQMOM workflow uses conditioning order {0,1} and
the positive discrete quadrature x={1,3}, wx={0.4,0.6}, y={-1,2},
wy={0.25,0.75}. Its operators are:

    advection: dx/dt = x, dy/dt = 0
    source:     dM_ij/dt = j M_(i+1,j-1)

The exact terminal trajectories are x(t)=exp(t)x(0) and
y(t)=y(0)+(exp(t)-1)x(0). The verifier advances a public CqmomState through
advance_cqmom_phase_space with Strang, `advsrc`, and `srcadv` splitting, compares
reconstructed raw moments with this discrete reference, and checks
second-order Strang versus first-order Lie convergence. `advsrc` and `srcadv` are also
required to differ at the coarse timestep.

The fixed final time is `T = 1.2` with macro steps
`0.3, 0.15, 0.075, 0.0375, 0.01875`. The reference is constructed directly
from the listed four terminal nodes after applying the exact trajectories; it
does not use a phase-space-evolved result. The error is the scaled infinity norm over
all configured raw moments, `max |computed - exact| / max(1, |exact|)`, and the
phase-space evolution moment-layout moments are checked separately. Phase-space advection RK23 tolerances
are `1e-12` relative and `1e-14` absolute with maximum step `0.01`, keeping
sub-integrator error below the measured splitting error. The observed-order
acceptance window is 0.2 around order two for Strang and order one for each Lie
ordering.

## Commuting CQMOM 2D

Configuration: `commuting_cqmom_2d.qmom`

This compact optional case verifies the same complete phase-space evolution path for a
two-dimensional CQMOM state with non-identity conditioning order. The
advection operator is diagonal:

```text
dx/dt = 0.16 x
dy/dt = -0.11 y
```

and the moment-source operator is:

```text
dM_ij/dt = -0.025 (i + j) M_ij
```

The exact solution is:

```text
M_ij(t) = exp((0.16 i - 0.11 j - 0.025 (i + j)) t) M_ij(0)
```

## Adaptive CQMOM State 2D

Configuration: `commuting_cqmom_state_adaptive_2d.qmom`

This case uses the same diagonal advection operator with a zero CQMOM source
callback through the recommended configured CQMOM workflow: `load_cqmom(...)`,
`CqmomCase::initial_state()`, and `advance_cqmom_phase_space(CqmomState, ...)`.
The hierarchy has
`node_counts = 2, 2`, which is the maximum branch capacity. With
`adaptive_branch_counts = true`, CQMOM automatically uses adaptive Wheeler for
marginal and conditional inversions. Reconstruction determines the active
branch counts and may return ragged terminal storage with fewer than four
terminal nodes. Phase-space advection advances the existing ragged quadrature,
preserves terminal paths, and reconstructs moments by summation from the moved
nodes; active counts cannot grow in this stage. The verifier checks that the
reconstructed full moments follow the advection analytical reference above
with the source term omitted.

## Automatic Scaling Through Adaptive CQMOM Evolution

The `automatic_scaling_cqmom_state_adaptive_{reference,rescaled}.qmom`
fixtures use the normal configured public workflow with
`automatic_scaling = true`, adaptive `{2,2}` capacity, and non-identity
conditioning order `{1,0}`. They represent the same ragged physical rule
`(-1,-1,0.10)`, `(1,-1,0.30)`, `(0.5,2,0.60)` in units `(1,1)` and
`(1e6,1e-6)`.

One Strang macro step applies constant physical node velocities
`(0.4,-0.25)` and the CQMOM-native source `dM_alpha/dt = 0.1 M_alpha`.
Constant advection is exact under the selected SSP-RK3 node step; the source
oracle is the independently evaluated SSP-RK3 polynomial
`P(z)=1+z+z^2/2+z^3/6`. The verifier checks final ragged paths, weights,
coordinates, and all `4 x 4` reconstructed moments against that manufactured
rule. Its source callback uses only `CqmomSourceContext::for_each_terminal_node`
and verifies that source-stage coordinates are in the configured physical
units after the first advection half-step, not normalized internal units.

The `2e-9` scaled tolerance covers reconstruction and source-stage roundoff;
the observed normalized errors are near machine precision. This is
scaling-enabled evolution interaction evidence: it verifies physical semantics
through evolution, not that scaling is numerically necessary for this
particular representation. The separate CQMOM feature-sensitive scaling
verification establishes that necessity on a deliberately difficult raw-moment
representation.

## CQMOM-Native Source Callback

Configuration: `cqmom_native_source.qmom`

This case verifies the public CQMOM physics callback workflow:

```text
load_cqmom(...)
    -> CqmomCase
    -> compute_cqmom(...)
    -> CqmomCase::initial_state()
    -> CqmomState
    -> advance_cqmom_phase_space(...)
    -> CqmomMomentSource
```

The configuration is a fixed 2D independent-product quadrature with
`conditioning_order = {0,1}`. The callback uses only the public
`CqmomSourceContext` API. On the first invocation it checks:

- `context.size()` is the expected number of phase-space-evolved CQMOM moments;
- `context.index(i)` follows the explicit phase-space multi-index order
  `{0,0}`, `{1,0}`, `{2,0}`, `{3,0}`, `{0,1}`, `{0,2}`, `{0,3}`, `{1,1}`,
  `{1,2}`, `{1,3}`;
- `context.value(i)` matches the corresponding initial phase-space-evolved moment;
- `context.hierarchy()` reports the configured hierarchy.

The returned rate is ordered exactly like the context:

```text
dM_alpha/dt = gamma M_alpha + B_alpha
```

where `gamma = -0.12` and `B_alpha` is constructed analytically from source
rates on the known terminal weights. Because this source changes only terminal
weights on the same support, the final state remains representable by the same
2D quadrature support. The verifier compares the reconstructed final moments
against the exact one-step SSP-RK3 update for this affine linear source.

## Combined CQMOM Phase-Space Advection and Native Moment Source

Configuration: `cqmom_combined_native_source.qmom`

This case verifies complete public CQMOM phase-space evolution through `CqmomState` with
both physics callbacks active in the same macro step:

```text
load_cqmom(...)
    -> CqmomCase
    -> compute_cqmom(...)
    -> CqmomCase::initial_state()
    -> CqmomState
    -> advance_cqmom_phase_space(...)
```

The advection velocity is the same diagonal 2D field used by the commuting
CQMOM case:

```text
dx/dt = 0.16 x
dy/dt = -0.11 y
```

so a moment `M_ij` is multiplied by `exp((0.16 i - 0.11 j) tau)` over a
advection substep of length `tau`. The native source callback applies a
non-zero diagonal moment source:

```text
dM_ij/dt = gamma M_ij
gamma = -0.09
```

The default scheme is Strang splitting, so the analytical one-step reference for
one macro step is:

```text
M_ij(dt) = exp(L_ij dt/2) P(gamma dt) exp(L_ij dt/2) M_ij(0)
L_ij = 0.16 i - 0.11 j
P(z) = 1 + z + z^2/2 + z^3/6
```

`P(z)` is the exact one-step polynomial for the SSP-RK3 source update applied
by the phase-space evolution framework. The source callback checks that its first context
contains the moments after the initial advection half-step, with the same
explicit multi-index ordering documented in the source-only case. A failure
there indicates that the public callback received stale values, reordered
indices, or rates that would be mapped back to the wrong phase-space-evolved moments.

## CQMOM-Native Source Callback Failures

Configuration: `cqmom_native_source_failures.qmom`

This case uses the same public `CqmomState` phase-space evolution workflow as the
successful CQMOM-native source case, but supplies invalid user callbacks:

- a callback returning the wrong number of rate values;
- a callback returning `NaN`;
- a callback returning infinity.

These mistakes simulate common physics-callback failures: forgetting that one
rate is required per `CqmomSourceContext` entry, or producing an invalid
floating-point derivative. Each case must fail during the moment-source
stage, report `PhaseSpaceEvolutionStatus::moment_source_failed`, preserve a moment-source diagnostic
message, and avoid silently proceeding to post-source reconstruction.

`CqmomMomentSource` has no explicit "reject this candidate" return type.
Source rejection in this public callback workflow is therefore represented by
invalid returned data or by phase-space evolution failure. Generic SSP-RK3 accept/reject
predicates exist in moment-source update options, but they are not returned by the
`CqmomMomentSource` callback and are not duplicated here.

## CQMOM Phase-Space Evolution Failure Propagation

Configuration: `cqmom_phase_space_evolution_failures.qmom`

This case checks how CQMOM phase-space evolution through `CqmomState` reports failures
through:

```text
load_cqmom(...)
    -> CqmomCase
    -> CqmomCase::initial_state()
    -> CqmomState
    -> advance_cqmom_phase_space(...)
```

The advection-callback failures mirror the QMOM case, but through the
recommended `CqmomState` workflow and the CQMOM-native source callback
signature. The verifier checks that wrong-dimension and non-finite velocity
returns are reported as phase-space evolution failures during the first
advection half-step, with meaningful diagnostic messages and no later-stage
success diagnostics.

The CQMOM case also repeats the RK23 minimum-step acceptance path through
`CqmomState`. The velocity callback is valid and finite, so the verifier checks
that phase-space evolution succeeds at the effective minimum step while preserving
rejected-step diagnostics.

The post-source reconstruction failure uses a finite, correctly sized
`CqmomMomentSource` callback that drives `M[0,0]` negative. The moment-source update
therefore succeeds, but the following CQMOM reconstruction is not physically
realizable. The verifier checks that phase-space evolution reports
`PhaseSpaceEvolutionStatus::inversion_failed`, identifies
`PhaseSpaceEvolutionStage::post_source_inversion`, preserves the successful source
diagnostics, records CQMOM diagnostics for the failed inversion, and stops
before the second advection half-step.

The same valid configuration also drives two further failure checks, which
differ in the entry point a caller uses rather than in how moments are stored:

- `verification_phase_space_evolution_cqmom_dense_phase_space_evolution_failures` advances
  the configured `CqmomMomentData` directly;
- `verification_phase_space_evolution_cqmom_lean_phase_space_evolution_failures` advances a
  `CqmomState` built from the same configured case.

For equivalent deterministic velocity, source, and post-source reconstruction
failures, both entry points should report the same high-level status, failed
stage, preserved diagnostics, and useful messages as the default `CqmomState`
workflow. The tests do not require identical internal diagnostic details.

Configuration: `cqmom_phase_space_evolution_inversion_failures.qmom`

This companion case verifies initial CQMOM inversion failure propagation
through the normal `CqmomState` workflow. The configuration loads finite
moments with `M[0,0] < 0`, so the initial recursive CQMOM reconstruction is
not realizable. The verifier checks that phase-space evolution distinguishes this from
advection or source failure, reports `PhaseSpaceEvolutionStage::initial_inversion`,
preserves CQMOM diagnostics and a realizability message, does not advance
time, and does not report later stages as completed.

The companion cases,
`verification_phase_space_evolution_cqmom_dense_phase_space_evolution_inversion_failures` and
`verification_phase_space_evolution_cqmom_lean_phase_space_evolution_inversion_failures`, reuse the
same invalid configured moments through the two entry points above. They document that CQMOM phase-space evolution through
`CqmomState` rejects an initially
unreconstructable moment state before invoking advection or source callbacks.

## Running

Configure and build the standalone test project first:

```bash
cd tests
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../install"
cmake --build build
```

Run individual cases:

```bash
ctest --test-dir build -R verification_phase_space_evolution_commuting_qmom --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_qmom_phase_space_evolution_failures --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_qmom_phase_space_evolution_inversion_failures --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_noncommuting_convergence --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_commuting_cqmom_2d --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_commuting_cqmom_state_adaptive_2d --output-on-failure
ctest --test-dir build -R 'verification_phase_space_evolution_cqmom_automatic_scaling_is_unit_invariant' --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_cqmom_native_source --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_cqmom_combined_native_source --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_cqmom_native_source_failures --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_cqmom_phase_space_evolution_failures --output-on-failure
ctest --test-dir build -R verification_phase_space_evolution_cqmom_phase_space_evolution_inversion_failures --output-on-failure
```

Run all complete phase-space evolution integrator analytical cases:

```bash
ctest --test-dir build -R verification_phase_space_evolution_ --output-on-failure
```

The verifier reports final exact and computed moments, maximum absolute moment
error, scaled moment error, zeroth-moment error, final reconstruction error,
and a refinement table with observed orders for convergence cases.

## Machine-Readable Output

The verifier can optionally write a whitespace-delimited text table:

```bash
mkdir -p tests/verification/phase_space_evolution/output

./build/phase_space_evolution_analytical_verification \
  tests/verification/phase_space_evolution/commuting_qmom.qmom \
  commuting_qmom \
  --output tests/verification/phase_space_evolution/output/commuting_qmom_refinement.txt

./build/phase_space_evolution_analytical_verification \
  tests/verification/phase_space_evolution/noncommuting_qmom.qmom \
  noncommuting_qmom \
  --output tests/verification/phase_space_evolution/output/noncommuting_qmom_refinement.txt
```

The output directory under `tests/verification/phase_space_evolution/output/` is ignored by
Git. Generated files are manual-analysis artifacts, not authoritative fixtures,
and are not committed by default.

The same files can be generated into the build tree with:

```bash
cmake --build build --target generate_phase_space_evolution_verification_outputs
```

The file format uses `std::setprecision(17)` and one row per refinement:

```text
case dt step_count final_time error_norm observed_order zeroth_moment_error reconstruction_error retry_count success
```

`error_norm` is the scaled infinity norm defined above. `observed_order` is:

```text
p = log(e_h / e_{h/2}) / log(2)
```

The first refinement row has `observed_order = nan` because there is no coarser
row in that table. `reconstruction_error` is the scaled final quadrature
reconstruction error. The generated rows are written from the same `RunSummary`
values used by the verifier's assertions; there is no separate hard-coded
results table.
