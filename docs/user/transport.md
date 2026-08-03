# Transport And Time Evolution

QMOM and CQMOM turn raw moments into quadrature nodes and weights. Transport
uses that quadrature representation to evolve a distribution in time, then
returns updated moments and quadrature for the next step or for output.

QMOMKit treats transport as a set of numerical utilities. Your application
supplies the phase-space velocity, moment sources, boundary exchange
contributions, timestep loop, and physical interpretation of coordinates that
are present in your model. QMOMKit supplies the public routines that advance
nodes, advance moment sources, reconstruct moments, and compose these pieces
into a transport step.

This division is useful for models with both smooth and source-like moment
evolution. Smooth effects such as growth, drift, acceleration, or relaxation can often be represented by moving
quadrature nodes through phase space. Internal or volumetric effects such as
breakage, aggregation, jumps, or externally computed rates are often more
naturally represented as raw-moment source terms. Phase-space boundary exchange
such as nucleation, injection, removal, or evaporation/condensation boundary
contributions can be supplied as a `MomentBoundaryContribution`. QMOMKit
combines the transport steps while leaving the physics in caller-supplied
callbacks.

QMOMKit separates transport physics into independent callbacks. Provide only
the callbacks required by your model: `PhaseSpaceVelocity` for continuous
phase-space motion, `MomentSource` or `CqmomMomentSource` for moment-space
evolution, and `MomentBoundaryContribution` or
`CqmomMomentBoundaryContribution` for phase-space boundary exchange. Examples
that show several callbacks together are illustrative combinations, not a
requirement for every simulation.

## Transport Workflow

A typical transport workflow alternates between the two representations:
quadrature is convenient for continuous node motion, while raw moments are
convenient for source updates and the next reconstruction.

```text
initial moments/quadrature
        ↓
node evolution, if present
        ↓
source or boundary-rate evolution, if present
        ↓
time integration
        ↓
updated quadrature/moments
```

The two main operations are:

- node advection: continuous evolution of quadrature abscissas in phase space;
- moment source update: evolution of raw moments through user-supplied source
  and optional boundary contribution rates.

For a full transport step, QMOMKit can compose these operations with a selected
transport scheme. Advanced applications can instead call the lower-level
composable routines directly when they manage individual transport operations
themselves.

If a model has only continuous phase-space motion, the advanced composable
operation is `advance_nodes(...)`, followed by moment reconstruction as needed.
If it has only moment-space evolution, the advanced composable operation is
`advance_moments_with_source(...)`. Use the higher-level
`advance_transport_qmom(...)` or `advance_transport_cqmom(...)` helpers when
you want QMOMKit to compose continuous motion and source-like moment evolution
in one macro step. Where a convenience helper signature includes a callback
slot for physics that is absent from your model, use a simple zero-rate
callback for that slot.

## Node Advection

Node advection moves quadrature node coordinates through a continuous
phase-space velocity field. Use it when the model says that existing
representative states move through phase space. In one-dimensional QMOM, this
advances the abscissas of the QMOM nodes. In multidimensional CQMOM, this
advances the terminal node coordinates in phase-space dimension order.

Users provide a `PhaseSpaceVelocity` callback when continuous phase-space
motion is part of the model:

```cpp
qmomkit::PhaseSpaceVelocity velocity =
    [](double time, const std::vector<double>& coordinates)
    {
        const double x = coordinates[0];
        return std::vector<double>{1.0 - x};
    };
```

The callback receives the current time and one node coordinate vector. It must
return `dxi/dt` with the same number of components. QMOMKit does not assign a
physical meaning to those coordinates; `x`, `diameter`, `velocity`, or any
other interpretation comes from your application and must match the
phase-space order used elsewhere. If a returned component is non-finite, the
diagnostic includes the evaluation time, zero-based phase-space dimension,
offending value, input coordinate vector, and complete returned velocity
vector.

Use the advanced `advance_nodes(...)` operation when you already have a
quadrature and want to evolve only the node positions:

```cpp
qmomkit::NodeAdvectionOptions node_options;
node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

const qmomkit::QmomNodeAdvectionResult advanced =
    qmomkit::advance_nodes(qmom_result, time, time_step, velocity, node_options);

if (!advanced.ok())
{
    std::cerr << "Node advection failed: " << advanced.message << '\n';
    return 1;
}
```

During node advection, quadrature weights are unchanged. If you need updated
raw moments after moving nodes, reconstruct them from the updated quadrature:

```cpp
const std::vector<double> moments =
    qmomkit::reconstruct_raw_moments(advanced.quadrature, moment_count);
```

For CQMOM, use the overload that accepts `CqmomQuadrature`. This is a
lower-level operation for workflows that only need to move node coordinates:

```cpp
const qmomkit::CqmomNodeAdvectionResult advanced =
    qmomkit::advance_nodes(cqmom_result.quadrature, time, time_step, velocity, node_options);
```

The full QMOM and CQMOM transport helpers described below keep accepted state
together with the reconstruction settings needed for the next step.

`NodeAdvectionOptions::scheme` selects how node trajectories are integrated:

- `NodeAdvectionScheme::rk23` is the default adaptive Bogacki-Shampine 3(2)
  method.
- `NodeAdvectionScheme::rk45` is the adaptive Dormand-Prince RK5(4)7M method
  from Table 2 of Dormand and Prince (1980).
- `NodeAdvectionScheme::ssp_rk3` takes one separate fixed three-stage step over
  the supplied timestep.

All three choices integrate the same continuous quadrature-node trajectory equation.
SSP-RK3 node advection does not retry, reduce the timestep, or perform
realizability checks. Its accuracy is controlled by the timestep chosen by the
caller. Realizability remains a property of reconstructed moments and later
quadrature inversions, not of the coordinate trajectory step itself.

Configuration-driven drivers can load the complete node-advection controls
with `qmomkit::config::load_node_advection_options(...)`. For RK23 and RK45,
`NodeAdvectionOptions::adaptive_options` carries scalar relative and absolute
tolerances, an initial step that defaults to `1.0e-6 s`, maximum step size,
minimum step size, and maximum internal accepted steps. The initial step may
be overridden but must be finite and strictly positive. QMOMKit does not
estimate it from the state, derivatives, or tolerances; accepted states use
the higher-order embedded solution and the methods reuse FSAL derivatives
internally.

## Moment Source Terms

Moment sources represent internal or volumetric changes that act directly on
raw moments. Use them when the model gives rates for moments, or when an event
changes weights or moments without being well described as smooth coordinate
motion. Examples include aggregation, fragmentation, relaxation, chemistry, or
any other model that can be written as:

```text
dm/dt = source(t, m)
```

Users provide a `MomentSource` callback when moment-space evolution is part of
the model:

```cpp
qmomkit::MomentSource source =
    [](double time, const std::vector<double>& moments)
    {
        std::vector<double> derivative(moments.size(), 0.0);
        derivative[0] = 0.0;
        derivative[1] = -0.1 * moments[1];
        return derivative;
    };
```

The returned vector must have the same size as the input moment vector. QMOMKit
does not know what each component means; the ordering must match the moment
ordering used by your QMOM setup. CQMOM transport uses the context-specific
callback below so multidimensional source terms can follow the configured
CQMOM moment layout.

Use the advanced `advance_moments_with_source(...)` operation when you want to
advance only the moment-source part of a workflow:

```cpp
const qmomkit::MomentSourceUpdateResult update =
    qmomkit::advance_moments_with_source(moments, time, time_step, source);

if (!update.ok())
{
    std::cerr << "Moment source update failed: " << update.message << '\n';
    return 1;
}
```

After a source update, reconstruct quadrature before applying another
continuous node-advection step:

```cpp
const qmomkit::QmomResult post_source =
    qmomkit::compute_qmom(qmomkit::MomentSet(update.moments), node_count);
```

## Boundary Contributions

Boundary contributions are optional additive raw-moment rates caused by
phase-space boundary exchange. Use `MomentBoundaryContribution` only when the
model describes nucleation, injection, removal, evaporation/condensation
exchange, or similar production/destruction through a phase-space boundary:

```text
dm_k/dt = source_k(t, m) + boundary_k(t, m)
```

For one-dimensional QMOM:

```cpp
qmomkit::MomentBoundaryContribution boundary =
    [](double time, const std::vector<double>& moments)
    {
        constexpr double rate = 0.02;
        constexpr double boundary_coordinate = 0.1;

        std::vector<double> values(moments.size(), 0.0);
        double power = 1.0;
        for (double& value : values)
        {
            value = rate * power;
            power *= boundary_coordinate;
        }

        (void)time;
        return qmomkit::MomentBoundaryRate{values};
    };
```

For an application-owned spatial operator or time integrator, the
four-argument QMOM evaluator combines both callbacks at one stage:

```cpp
const qmomkit::QmomDiscontinuousMomentRateResult rates =
    qmomkit::evaluate_qmom_discontinuous_moment_rates(time, moments, source, boundary);
```

It evaluates `source(time, moments) + boundary(time, moments)` without
advancing time or reconstructing a quadrature; the source callback runs before
the boundary callback. Use `advance_moments_with_source(...)` when QMOMKit
should perform SSP-RK3 source integration, or the transport helpers when it
should also own the configured operator-splitting sequence.

For CQMOM, use `CqmomMomentBoundaryContribution`. It receives the same
`CqmomSourceContext` as `CqmomMomentSource` and returns one additive rate per
context entry:

```cpp
qmomkit::CqmomMomentBoundaryContribution boundary =
    [](double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> values(context.size(), 0.0);

        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);
            values[entry] = /* boundary exchange rate for M_order */;
        }

        (void)time;
        return qmomkit::CqmomMomentRate{values};
    };
```

Boundary contributions are not a generic PDE boundary-condition framework.
They do not define mesh faces, normals, reflection, clipping, or boundary
geometry. They are runtime physics callbacks, just like velocities and moment
sources, and are not stored in `QmomCase`, `CqmomCase`, or `CqmomState`.

When a boundary callback is supplied to a transport helper, it is evaluated
with the moment source during each SSP-RK3 source stage:

```text
effective_rate(t, m) = source(t, m) + boundary(t, m)
```

The boundary-aware transport overloads take an explicit options argument to
keep overload resolution unambiguous.

## Physics Callbacks

Physics callbacks are the main place where an application extends QMOMKit.
QMOMKit handles moment storage, moment ordering, quadrature reconstruction,
state updates, and transport integration. Your application supplies the
physical model: continuous phase-space motion and discontinuous or source-like
moment rates.

The callback extension points are independent. A model may use only
`PhaseSpaceVelocity`, only `MomentSource`, source plus boundary contribution,
velocity plus source, or all three categories. Choose the public routine that
matches the needed composition, and keep absent physics explicit as either no
callback in lower-level workflows or a zero-rate callback in full transport
helpers whose signature requires one.

Most applications change three things:

- configuration: phase-space names, node counts, initial condition, output,
  and time-integration settings;
- physics callbacks: whichever velocity fields, moment source terms, or
  boundary exchange contributions the model actually needs;
- transport options: timestep controls, tolerances, and splitting scheme.

Most applications can treat internal moment storage, reconstruction logic,
quadrature layout, and CQMOM representation details as QMOMKit-owned parts of
the workflow.

### Continuous Sources

Continuous physics moves quadrature nodes through phase space. Provide a
`PhaseSpaceVelocity` callback only when the physical model can be written as a
coordinate velocity:

```text
dx/dt = velocity(t, x)
```

The callback receives the current time and one terminal node coordinate vector
in phase-space dimension order. It returns one velocity component per
coordinate. This is the right place for drift, growth, relaxation, acceleration,
or any effect that is naturally represented by moving nodes.

### Discontinuous Sources

Discontinuous physics updates raw moments directly. Provide a `MomentSource`
only for internal or volumetric moment evolution and, when needed, a
`MomentBoundaryContribution` for phase-space boundary exchange. For CQMOM use
`CqmomMomentSource` and `CqmomMomentBoundaryContribution`. The combined model
is:

```text
dM_alpha/dt = source_alpha(t, M) + boundary_alpha(t, M)
```

This is the right place for effects that change moment values rather than
continuously moving existing nodes. The transport helpers compose continuous
and discontinuous parts with the selected scheme. Boundary contributions enter
the existing source update; they do not add a separate transport operator.

### CqmomSourceContext

`CqmomSourceContext` describes the CQMOM moments that need source or boundary
contribution rates during a moment-source update. It is passed to
`CqmomMomentSource` and `CqmomMomentBoundaryContribution` callbacks so each
callback can compute rates for the exact moment entries required by the
configured hierarchy.

Use:

- `context.time()` for the current SSP-RK source-stage time;
- `context.size()` for the number of moment entries that need derivatives;
- `context.moment_count()` as a synonym for that source moment count;
- `context.dimension_count()` for the hierarchy dimension count;
- `context.indices()` and `context.values()` for the context-ordered
  multi-index and value arrays;
- `context.index(i)` for the `MultiIndex` of entry `i`;
- `context.value(i)` for the current raw moment value of entry `i`;
- `context.position(alpha)` or `context.value(alpha)` when a source term needs
  a particular retained moment by multi-index;
- `context.hierarchy()` for the CQMOM hierarchy associated with the state.

The returned source values correspond to the context ordering. For entry `i`,
compute the derivative for `context.index(i)` from the current state and write
that derivative to return entry `i`. The multi-index uses phase-space
dimension order.

QMOMKit may request only the moments required by the configured hierarchy, so
CQMOM source and boundary callbacks are most robust when they iterate over the
context entries instead of assuming a dense rectangular tensor.

CQMOM source callbacks usually look like this:

```cpp
qmomkit::CqmomMomentSource source =
    [](double time, const qmomkit::CqmomSourceContext& context)
    {
        const qmomkit::PhaseSpaceHierarchy& hierarchy = context.hierarchy();
        std::vector<double> derivatives(context.size(), 0.0);

        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);
            const double moment = context.value(entry);

            double total_order = 0.0;
            for (std::size_t dimension = 0;
                 dimension < hierarchy.number_of_dimensions();
                 ++dimension)
            {
                total_order += static_cast<double>(order[dimension]);
            }

            const double rate = -0.05 * (1.0 + total_order);
            derivatives[entry] = rate * moment;
        }

        (void)time;
        return qmomkit::CqmomMomentRate{derivatives};
    };
```

This pattern is dimension-independent: it uses the hierarchy dimension count
instead of assuming a two-dimensional phase space.

## Bounded CQMOM Support During Transport

When a `CqmomCase` or `CqmomState` carries coordinate support from
`phase_space.support_lower` and `phase_space.support_upper`, every CQMOM
reconstruction for that coordinate uses the same support-aware acceptance check.
This includes initial inversion, post-source inversion, current active-count
preservation, and adaptive growth attempts.

Fixed branch counts are strict: if the configured count is not realizable on the
declared interval, the stage fails and the returned state rolls back to the
input state. Adaptive branch counts first try to preserve the current count,
then reduce to a feasible lower count when needed. Growth trials are allowed
only toward the configured maximum; if a growth candidate violates support, the
current feasible count is retained.

After bounded realizability succeeds, QMOMKit verifies the returned nodes using
the production realizability tolerance. Boundary nodes and nodes inside the
interval within that tolerance are accepted unchanged. Nodes outside the
interval by more than the tolerance are rejected. QMOMKit does not clip,
project, or otherwise move nodes to fit the interval.

Successful post-source stages store canonical moments reconstructed from the
accepted quadrature and active counts. A source candidate whose higher-order
moments would require unsupported nodes is not preserved as an infeasible
sequence; the returned state contains the accepted bounded quadrature closure.

## Time Integration

QMOMKit uses two user-visible integration roles:

- node-advection integration advances quadrature coordinates through
  continuous phase-space velocities using default adaptive RK23, adaptive RK45,
  or selected fixed-stage SSP-RK3;
- source integration advances raw moments through discontinuous or source-like
  moment updates with SSP-RK3.

Most users interact with these through `NodeAdvectionOptions`,
`MomentSourceUpdateOptions`, or the higher-level transport options:

```cpp
qmomkit::QmomTransportOptions options;
options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
options.node_options.adaptive_options.relative_tolerance = 1.0e-7;
options.source_options.ssp_rk3_options.max_retries = 4;
```

The full transport helpers support three schemes:

```text
strang
cd
dc
```

`strang` is the default. It applies a continuous half-step, a source full-step,
and a second continuous half-step. This is second-order Strang splitting and is
the usual choice when both continuous node motion and source-like moment
evolution matter.

`cd` applies continuous evolution first, then the source update. `dc` applies
the source update first, then continuous evolution. These first-order Lie
splitting schemes can be useful for comparisons, simple workflows, or cases
where your application deliberately wants a fixed splitting order.

## Full Transport Step

Use `QmomCase::initial_state()` to construct an accepted one-dimensional state,
then use `advance_transport_qmom(...)` to advance that state by one transport
step:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
    const std::size_t node_count = 2;
    std::vector<double> moments{1.0, 0.8, 2.8, 4.4};

    const qmomkit::PhaseSpaceVelocity velocity =
        [](double, const std::vector<double>& coordinates)
        {
            return std::vector<double>{1.0 - coordinates[0]};
        };

    const qmomkit::MomentSource source =
        [](double, const std::vector<double>& current_moments)
        {
            return std::vector<double>(current_moments.size(), 0.0);
        };

    qmomkit::QmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::strang;

    qmomkit::QmomState state =
        qmomkit::QmomCase{qmomkit::MomentSet(moments), node_count, options.qmom_options}
            .initial_state();

    double time = 0.0;
    const double final_time = 0.1;
    const double time_step = 0.01;

    while (time < final_time)
    {
        const double step = std::min(time_step, final_time - time);
        const qmomkit::QmomTransportResult result =
            qmomkit::advance_transport_qmom(
                state, node_count, time, step, velocity, source, options);

        if (!result.ok())
        {
            std::cerr << "Transport failed: " << result.message << '\n';
            return 1;
        }

        state = result.state;
        moments = state.moments();
        time = result.final_time;
    }
}
```

Use `advance_transport_qmom(...)` when QMOMKit should own one configured
Strang or Lie splitting step. If the surrounding application owns stage
orchestration, coupling, or acceptance, call the one-dimensional state stages
directly:

```cpp
const qmomkit::QmomTransportResult first =
    qmomkit::advance_qmom_continuous(
        state, time, 0.5 * time_step, velocity, options);

if (!first.ok())
{
    // Handle failure.
}

const qmomkit::QmomTransportResult source_stage =
    qmomkit::advance_qmom_discontinuous(
        first.state, node_count, time + 0.5 * time_step, time_step, source, options);

if (!source_stage.ok())
{
    // Handle failure.
}

const qmomkit::QmomTransportResult second =
    qmomkit::advance_qmom_continuous(
        source_stage.state, time + 0.5 * time_step, 0.5 * time_step, velocity, options);
```

Continuous stages move the accepted quadrature topology and reconstruct the
stored moment vector without another QMOM inversion. Discontinuous stages own
the source update and the same fixed or adaptive post-source reconstruction
policy used by the full transport helper. Application-managed workflows should
pass `QmomState` objects between stages instead of reconstructing quadrature or
mutating moments manually.

If the model also has phase-space boundary exchange, define a
`MomentBoundaryContribution` and call the boundary-aware overload with explicit
transport options:

```cpp
const qmomkit::MomentBoundaryContribution boundary =
    [](double, const std::vector<double>& current_moments)
    {
        constexpr double injection_rate = 0.02;
        constexpr double injection_coordinate = 0.1;

        std::vector<double> rate(current_moments.size(), 0.0);
        double power = 1.0;
        for (double& value : rate)
        {
            value = injection_rate * power;
            power *= injection_coordinate;
        }
        return qmomkit::MomentBoundaryRate{rate};
    };

const qmomkit::QmomTransportResult result =
    qmomkit::advance_transport_qmom(
        state, node_count, time, step, velocity, source, boundary, options);
```

`QmomState` exposes read-only `moments()` and `quadrature()` accessors. Its
stored moments are reconstructed from its accepted quadrature, and
`QmomState::active_node_count()` is derived from the accepted quadrature size.
Users may copy and move accepted states for retry or in-process checkpointing,
but direct construction from independent moments and quadrature is intentionally
unavailable; no public restart/import factory is currently provided.
In adaptive QMOM, `node_count` remains the configured maximum. Continuous node
advection preserves the active count and reconstructs the stored moment vector
from the advected quadrature without another Wheeler inversion. After an
accepted source update, transport first reconstructs at the current active count
and then tries at most one additional node when the configured maximum permits
it. Fixed-count QMOM remains strict at `node_count`.

The one-dimensional `QuadratureWriter` writes one row per active QMOM node:

```text
time node_index weight x0
```

Create it with the configured maximum node count. Fixed and adaptive QMOM use
the same writer, and adaptive output writes only active nodes without padding
inactive entries.

For CQMOM, use `advance_transport_cqmom(...)` with the `CqmomState` obtained
from the loaded CQMOM case when QMOMKit should own one complete configured
splitting step:

```cpp
const auto configured = qmomkit::config::load_cqmom("configuration.qmom");
if (!configured)
{
    // Report configured.error(), then stop.
}

const qmomkit::CqmomCase cqmom_case = configured.value();

const qmomkit::CqmomResult initial_quadrature = qmomkit::compute_cqmom(cqmom_case);
if (!initial_quadrature.ok())
{
    // Report initial_quadrature.message, then stop.
}

qmomkit::CqmomState state = cqmom_case.initial_state();

qmomkit::CqmomTransportOptions options;
options.scheme = qmomkit::TransportScheme::strang;

const qmomkit::CqmomMomentSource cqmom_source =
    [](double, const qmomkit::CqmomSourceContext& context)
    {
        return qmomkit::CqmomMomentRate{
            std::vector<double>(context.size(), 0.0)};
    };

const qmomkit::CqmomTransportResult result =
    qmomkit::advance_transport_cqmom(
        state, time, time_step, velocity, cqmom_source, options);

if (result.ok())
{
    state = result.state;
}
```

The CQMOM boundary-aware overload has the same shape, with an explicit options
argument:

```cpp
const qmomkit::CqmomMomentBoundaryContribution cqmom_boundary =
    [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);
            rate[entry] = /* phase-space boundary exchange rate for M_order */;
        }
        return qmomkit::CqmomMomentRate{rate};
    };

const qmomkit::CqmomTransportResult result_with_boundary =
    qmomkit::advance_transport_cqmom(
        state, time, time_step, velocity, cqmom_source, cqmom_boundary, options);
```

Applications that own C/D stage orchestration can use the public CQMOM stage
helpers directly:

```cpp
const qmomkit::CqmomTransportResult first =
    qmomkit::advance_cqmom_continuous(
        state, time, 0.5 * time_step, velocity, options);

if (!first.ok())
{
    // Handle failure.
}

const qmomkit::CqmomTransportResult source_stage =
    qmomkit::advance_cqmom_discontinuous(
        first.state, time + 0.5 * time_step, time_step, cqmom_source, options);

if (!source_stage.ok())
{
    // Handle failure.
}

const qmomkit::CqmomTransportResult second =
    qmomkit::advance_cqmom_continuous(
        source_stage.state, time + 0.5 * time_step,
        0.5 * time_step, velocity, options);
```

The boundary-aware discontinuous stage has the same shape as the high-level
boundary overload:

```cpp
const qmomkit::CqmomTransportResult source_with_boundary =
    qmomkit::advance_cqmom_discontinuous(
        first.state, time + 0.5 * time_step, time_step,
        cqmom_source, cqmom_boundary, options);
```

Continuous CQMOM stages preserve terminal paths and active branch counts, move
the accepted terminal quadrature nodes, and reconstruct moments from that
advected quadrature without a post-stage CQMOM inversion. Discontinuous stages
advance the CQMOM transport moment layout, then apply the same recursive fixed
or adaptive reconstruction policy used by `advance_transport_cqmom(...)`.
Source stages ignore node-advection options, continuous stages ignore source
integration options, and standalone stage helpers do not use the splitting
scheme option. On failure, the returned state is the stage input state; the
returned quadrature, when present, is diagnostic unless `result.ok()` is true.

Here `compute_cqmom(...)` computes the initial quadrature, while
`CqmomCase::initial_state()` creates the runtime state that transport advances.

Adaptive branch counts configured on the `CqmomCase` are preserved by
`CqmomState` transport. In adaptive mode, the hierarchy node counts remain the
maximum branch counts; the state also retains the current branch-specific
active counts associated with its ragged quadrature. Continuous node advection
preserves terminal paths and therefore preserves active counts. Moments are
reconstructed by summation from the advected quadrature; this stage does not
grow active branch counts.

At reconstruction boundaries, existing branches use their current active count
as the ceiling, so active counts may stay the same or decrease but cannot grow.
After an accepted discontinuous source update, post-source inversion first
tests the current active count; if that count remains valid and the configured
maximum is larger, it tries exactly one additional node. Growth is accepted only
when that trial returns the requested `Nactive + 1` nodes. If a parent topology
change creates a genuinely new descendant branch, the descendant begins from
its own configured maximum node count, and adaptive Wheeler determines its
initial active count.

Adaptive `CqmomState` construction seeds active counts from the initial CQMOM
inversion; if that inversion fails, state construction fails instead of creating
a state with inconsistent adaptive metadata. Failed CQMOM transport through
`CqmomState` returns the previous state metadata so callers can retry from the
last accepted step.

`CqmomTransportOptions` contains only transport-stage policies: node advection,
source integration, and the C/D sequencing scheme.
CQMOM inversion policies such as Wheeler options, adaptive branch counts,
diagnostic level, and Vandermonde solver belong to `CqmomCase` and are preserved
by `CqmomState`. Successful CQMOM transport through `CqmomState` persists active
adaptive topology for the next call. Explicit dense and lean case types remain
available for advanced static inversion and specialized initialization, but
normal CQMOM transport uses `CqmomCase::initial_state()` and
`advance_transport_cqmom(...)` with the storage-independent `CqmomState` API.

For storage-independent CQMOM physics, use `CqmomMomentSource` and
`CqmomMomentBoundaryContribution`. `CqmomSourceContext` exposes the CQMOM moment
layout selected for transport, with the same indices, values, ordering, hierarchy,
and lookup behavior for `CqmomState` with dense or lean moment storage.

`result.ok()` identifies whether returned moments or quadrature describe an
accepted transport step. On failure, `result.message` and `result.diagnostics`
indicate whether the inversion, node advection, moment source update, or
moment reconstruction failed.

## Configuration Workflow

Configuration files can select initialization, quadrature settings, transport
scheme, and output options for drivers that read those sections. The transport
scheme is set in:

```text
[time_integration]
scheme = strang
```

Accepted values are:

```text
strang
cd
dc
```

Node-advection integration is configured separately because it controls the
continuous trajectory integrator, not the operator-splitting scheme:

```text
[time_integration.node_advection]
scheme = rk23
```

Accepted values are:

```text
rk23
rk45
ssp_rk3
```

`rk23` is the default adaptive Bogacki-Shampine 3(2) method. `rk45` is the
adaptive Dormand-Prince RK5(4)7M pair from Table 2 of Dormand and Prince
(1980). `ssp_rk3` remains a separate fixed-stage method. The embedded methods
use their higher-order state for accepted steps and reuse FSAL derivatives.

Adaptive controls are configured in a nested section:

```text
[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-7
absolute_tolerance = 1.0e-10
initial_step = 1.0e-4
max_step = 1.0e-2
min_step = 1.0e-12
max_steps = 200000
```

Omitted adaptive-RK keys keep the `AdaptiveRkOptions` defaults:
`relative_tolerance = 1.0e-6`, `absolute_tolerance = 1.0e-9`,
`initial_step = 1.0e-6 s`, `max_step = full integration interval`,
`min_step = machine-dependent automatic minimum`, and
`max_steps = 1000000000`. `initial_step` may be overridden and must be finite
and strictly positive whenever `rk23` or `rk45` is selected; QMOMKit does not
estimate it from the state, derivatives, or tolerances. Startup summaries in the transport examples report the effective
node-advection scheme and controls. For fixed-step schemes, the summary reports
only the selected scheme. Discontinuous source integration is unchanged.

Output can be requested with:

```text
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
```

Configuration files do not define the velocity callback, source callback,
boundary contribution callback, or outer time loop. Those remain application
code. Drivers can read the scheme with the public configuration API and assign
it to `QmomTransportOptions` or `CqmomTransportOptions`. See
`configuration.md` for the configuration syntax and the example files under:

```text
examples/qmom_1d_continuous_breakage/
examples/cqmom_transport/
```

## Practical Guidance

Coordinate order links the transport callback to the configured model. A
`PhaseSpaceVelocity` callback receives coordinates in phase-space dimension
order, the same order used by CQMOM moments, node counts, and output.

Each moment callback must return the same number of rates as the moments it
receives. A `MomentSource` callback returns the same number of moment
derivatives it receives. A `MomentBoundaryContribution` callback returns the
same number of additive boundary rates. A `CqmomMomentSource` or
`CqmomMomentBoundaryContribution` callback returns one rate per
`CqmomSourceContext` entry.

For CQMOM source and boundary rates, the context order defines the return
order. For entry `i`, compute the rate for `context.index(i)` from the current
value `context.value(i)`, then store it in return entry `i`.

Source terms update raw moments. Quadrature nodes are rebuilt from those
moments before another node-advection step.

Boundary contributions are added to the moment-source rate and integrated by
the same SSP-RK3 source-update path. Source retries, acceptability predicates,
and post-source inversion behavior are shared.

The evolved moment set still needs to support the requested quadrature.
Transport can fail if updated moments no longer support the requested QMOM or
CQMOM reconstruction.

Callbacks are application code. QMOMKit does not provide meshes, boundary
conditions, boundary geometry, physical source models, or units. The library
advances the moment and quadrature data provided by the surrounding workflow.
