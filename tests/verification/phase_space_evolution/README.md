# Phase-Space Evolution Analytical Verification

These cases verify the public phase-space evolution facility as a complete path:
configuration loading, QMOM/CQMOM inversion, advection operator calls,
raw-moment reconstruction, moment-source operator calls, repeated moment/quadrature
representation changes, diagnostics, and phase-space evolution scheme behavior for
`strang`, `advsrc`, and `srcadv`. They are automated analytical tests that complement
the user-facing phase-space evolution examples under `examples/qmom_1d_breakage_advection/`
and `examples/cqmom_2d_minimal/`.

All initial QMOMKit cases are loaded from `.qmom` files through the public
configuration API, and the fixtures define only the initial case. The splitting
scheme (`strang`, `advsrc`, or `srcadv`), the RK23/SSP-RK3 options, operator
coefficients, exact references, timestep plans, and tolerances are set in
`phase_space_evolution_analytical.cpp`, so each case controls them explicitly.
The configuration format does not define physics callbacks, model coefficients,
or timestep plans.

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

This companion case covers a QMOM moment sequence that is not realizable before
phase-space evolution can begin. The input moments are finite and load
successfully, but `m0 < 0`, so `QmomCase::initial_state()` throws and no
accepted state can enter phase-space evolution. The verifier checks that it
throws.

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
`0.3, 0.15, 0.075, 0.0375, 0.01875`. The reference applies the exact
trajectories to the initial CQMOM quadrature that `compute_cqmom` recovers from
the configured moments; it does not use a phase-space-evolved result. That
initial quadrature is not compared with the four nodes listed above, so the
reference depends on the initial reconstruction. The error is the larger of the
scaled infinity norms, `max |computed - exact| / max(1, |exact|)`, over the
state's layout moments and over the reconstructed raw-moment tensor. Phase-space advection RK23 tolerances
are `1e-12` relative and `1e-14` absolute with maximum step `0.01`, keeping
sub-integrator error below the measured splitting error. The penultimate and
finest observed orders must be at least 1.8 for Strang and at least 0.8 for
each Lie ordering.

## Commuting CQMOM 2D

Configuration: `commuting_cqmom_2d.qmom`

This compact case verifies the same complete phase-space evolution path for a
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
callback through the recommended configured CQMOM workflow: `load_cqmom_prepared(...)`,
`CqmomPreparedCase::initial_state()`, and `advance_cqmom_phase_space(CqmomState, ...)`.
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

One Strang macro step applies constant node velocities `(0.4,-0.25)` in the
reference units, scaled by the fixture units, and the CQMOM-native source `dM_alpha/dt = 0.1 M_alpha`.
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
load_cqmom_prepared(...)
    -> CqmomPreparedCase
    -> compute_cqmom(...)
    -> CqmomPreparedCase::initial_state()
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
load_cqmom_prepared(...)
    -> CqmomPreparedCase
    -> compute_cqmom(...)
    -> CqmomPreparedCase::initial_state()
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
load_cqmom_prepared(...)
    -> CqmomPreparedCase
    -> CqmomPreparedCase::initial_state()
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

The same valid configuration also drives
`verification_phase_space_evolution_cqmom_case_phase_space_evolution_failures`.
It advances a `CqmomState` built directly from the configured `CqmomCase`
rather than created by a prepared case. That state has no accepted quadrature,
so phase-space evolution performs the initial CQMOM inversion itself. For the
same velocity-callback and post-source reconstruction failures, plus a source
callback that returns a wrong-sized rate vector, it checks the same high-level
status, failed stage, preserved diagnostics, and useful messages, and that the
moments are unchanged after a failed advection. It does not require identical
internal diagnostic details. It also checks that diagnostic-level,
Vandermonde-solver, and adaptive-branch-count settings on the case reach the
CQMOM reconstructions during phase-space evolution; for this probe a copy of the
case is set to `basic` diagnostics, the dense solver, and adaptive branch
counts.

Configuration: `cqmom_phase_space_evolution_inversion_failures.qmom`

This companion case verifies initial CQMOM inversion failure propagation
through `CqmomState`. The configuration loads finite moments with
`M[0,0] < 0`, so the initial recursive CQMOM reconstruction is not realizable.
`CqmomPreparedCase::initial_state()` would throw for these moments, so the
state is built directly from the configured `CqmomCase`, with no accepted
quadrature, and the failure surfaces in phase-space evolution. The verifier
checks that phase-space evolution distinguishes this from advection or source
failure, reports `PhaseSpaceEvolutionStage::initial_inversion`, preserves
CQMOM diagnostics and a realizability message, does not advance time, and does
not report later stages as completed.

## Boundary Contributions

`qmom_boundary_contribution` reuses `commuting_qmom.qmom` with zero velocity and
one `srcadv` step of `0.04`. A boundary callback injects mass at `0.25` with
rate `0.08` and removes it at `1.1` with rate `-0.03`. Because the rates are
constant, the exact moments are the initial moments plus `dt` times the rate;
the verifier compares the state's moments and the moments reconstructed from its
final quadrature with that reference. A second run from the same initial state
adds a constant moment source to the same boundary rates and compares against
the corresponding sum.

`cqmom_boundary_contribution` reuses `commuting_cqmom_2d.qmom` with zero
velocity and one Strang step of `0.02`. The boundary rate of each stored moment
is the moment of prescribed weight rates at the nodes of the initial quadrature,
so the exact moments are again the initial moments plus `dt` times that rate.
The callback also checks that its source context has the expected size,
hierarchy, conditioning order `{1,0}`, and index order.

## Linear Advection Against a Closed Form

`qmom_linear_advection_manufactured.cpp` uses the velocity `u(x) = x`, for which
every node follows `x(t) = x(0) exp(t)` and advection leaves the weights
unchanged, so the raw moments follow `M_k(t) = M_k(0) exp(k t)`. The initial
moments come from three chosen nodes and weights, not from a configuration file.

- `verification_qmom_advection_matches_exponential_phase_evolution`: after
  advection the positions must be within `2e-12` of `x(0) exp(0.1)`, the weights
  must be unchanged bit for bit, and the rebuilt moments must agree to `2e-11`
  relative.
- `verification_qmom_phase_space_source_evolution_matches_manufactured_solution`:
  with a source that adds `0.1` to every moment per unit time, the exact
  solution is `M_0(t) = M_0(0) + 0.1 t` and
  `M_k(t) = M_k(0) exp(k t) + 0.1 (exp(k t) - 1) / k` for `k >= 1`. The error
  must fall at every refinement, and over the last two refinements the observed
  order must exceed `1.7` for Strang and lie between `0.7` and `1.3` for the
  two first-order orderings.

## Running

Configure and build the test suite from the repository root, in Debug as
described in [Testing](../../../docs/developer/testing.md#debug-build-requirement):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build --parallel
```

Run individual cases:

```bash
ctest --test-dir build -R '^verification_phase_space_evolution_commuting_qmom$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_qmom_boundary_contribution$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_qmom_phase_space_evolution_failures$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_qmom_phase_space_evolution_inversion_failures$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_noncommuting_convergence$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_noncommuting_cqmom_convergence$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_commuting_cqmom_2d$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_boundary_contribution$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_commuting_cqmom_state_adaptive_2d$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_automatic_scaling_is_unit_invariant$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_native_source$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_combined_native_source$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_native_source_failures$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_phase_space_evolution_failures$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_phase_space_evolution_inversion_failures$' --output-on-failure
ctest --test-dir build -R '^verification_phase_space_evolution_cqmom_case_phase_space_evolution_failures$' --output-on-failure
ctest --test-dir build -R '^verification_qmom_advection_matches_exponential_phase_evolution$' --output-on-failure
ctest --test-dir build -R '^verification_qmom_phase_space_source_evolution_matches_manufactured_solution$' --output-on-failure
```

Run all phase-space evolution verification cases in this directory:

```bash
ctest --test-dir build -R '^verification_(phase_space_evolution_|qmom_advection_matches|qmom_phase_space_source)' --output-on-failure
```

The verifier is silent on success. On failure it writes a per-run error table
(scheme, timestep, steps, final time, error norm, observed order, absolute and
scaled moment errors, reconstruction errors, zeroth-moment error, retry count,
and success) and the failure messages to standard error.
