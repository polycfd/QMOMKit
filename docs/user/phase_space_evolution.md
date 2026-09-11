# Phase-Space Evolution and Time Integration

QMOM and CQMOM turn raw moments into quadrature nodes and weights. Phase-space evolution
uses that quadrature representation to evolve a distribution in time, then
returns updated moments and quadrature for the next step or for output.

QMOMKit treats phase-space evolution as a set of numerical utilities. Your application
supplies the phase-space velocity, moment sources, boundary exchange
contributions, timestep loop, and physical interpretation of coordinates that
are present in your model. QMOMKit supplies the public routines that advance
nodes, advance moment sources, reconstruct moments, and compose these pieces
into a phase-space evolution step.

This division is useful for models with both smooth and source-like moment
evolution. Smooth effects such as growth, drift, acceleration, or relaxation can often be represented by moving
quadrature nodes through phase space. Internal or volumetric effects such as
breakage, aggregation, jumps, or externally computed rates are often more
naturally represented as raw-moment source terms. Phase-space boundary exchange
such as nucleation, injection, removal, or evaporation/condensation boundary
contributions can be supplied as a `MomentBoundaryContribution`. QMOMKit
combines the phase-space evolution steps while leaving the physics in caller-supplied
callbacks.

QMOMKit separates phase-space evolution physics into independent callbacks. Provide only
the callbacks required by your model: `PhaseSpaceVelocity` for phase-space
advection, `MomentSource` or `CqmomMomentSource` for moment-source evolution,
and `MomentBoundaryContribution` or
`CqmomMomentBoundaryContribution` for phase-space boundary exchange. Examples
that show several callbacks together are illustrative combinations, not a
requirement for every simulation.

## Phase-Space Evolution Workflow

A typical phase-space evolution workflow alternates between the two representations:
quadrature is convenient for phase-space advection, while raw moments are
convenient for moment-source updates and the next reconstruction.

```text
initial moments/quadrature
        ↓
node evolution, if present
        ↓
moment-source or boundary-rate evolution, if present
        ↓
time integration
        ↓
updated quadrature/moments
```

The two main operations are:

- node advection: phase-space advection of quadrature abscissas in phase space;
- moment-source update: evolution of raw moments through user-supplied source
  and optional boundary contribution rates.

For a full phase-space evolution step, QMOMKit can compose these operations with a selected
phase-space evolution scheme. Advanced applications can instead call the lower-level
composable routines directly when they manage individual phase-space evolution operations
themselves.

If a model has only phase-space advection, the advanced composable
operation is `advance_nodes(...)`, followed by moment reconstruction as needed.
If it has only moment-source evolution, the advanced composable operation is
`advance_moments_with_source(...)`. Use the higher-level
`advance_qmom_phase_space(...)` or `advance_cqmom_phase_space(...)` helpers when
you want QMOMKit to compose phase-space advection and source-like moment evolution
in one macro step. Where a convenience helper signature includes a callback
slot for physics that is absent from your model, use a simple zero-rate
callback for that slot.

## Node Advection

Node advection moves quadrature node coordinates through a phase-space
velocity field. Use it when the model says that existing
representative states move through phase space. In one-dimensional QMOM, this
advances the abscissas of the QMOM nodes. In multidimensional CQMOM, this
advances the terminal node coordinates in phase-space dimension order.

Users provide a `PhaseSpaceVelocity` callback when phase-space advection
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
return $\mathrm{d}\xi/\mathrm{d}t$ with the same number of components. QMOMKit does not assign a
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

Configured abscissa bounds do not alter node trajectories: node advection does
not clip, reflect, project, or reject a node solely because it crosses a
configured interval. Bounds are instead checked when QMOM or CQMOM reconstructs
quadrature from moments, including source-stage reconstructions and accepted
post-source inversions.

For CQMOM, use the overload that accepts `CqmomQuadrature`. This is a
lower-level operation for workflows that only need to move node coordinates:

```cpp
const qmomkit::CqmomNodeAdvectionResult advanced =
    qmomkit::advance_nodes(cqmom_result.quadrature, time, time_step, velocity, node_options);
```

The full QMOM and CQMOM phase-space evolution helpers described below keep accepted state
together with the reconstruction settings needed for the next step.

`NodeAdvectionOptions::scheme` selects how node trajectories are integrated:

- `rk23`: the adaptive Bogacki-Shampine 3(2) scheme [[Bogacki & Shampine, 1989](https://www.doi.org/10.1016/0893-9659(89)90079-7)];
- `rk45`: the adaptive Dormand-Prince RK5(4)7M scheme [[Dormand & Prince, 1980](https://www.doi.org/10.1016/0771-050X(80)90013-3)];
- `ssp_rk3`: the fixed-stage SSP-RK3 scheme [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)].

All three choices integrate the same advection quadrature-node trajectory equation.
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
raw moments. They are integrated using the  explicit third-order strong-stability-preserving 
Runge-Kutta step [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)].
Use moment sources when the model gives rates for moments, or when an event
changes moments without being well described as smooth coordinate
motion. Examples include aggregation, fragmentation, relaxation, chemistry, or
any other model that can be written as:
$$\frac{\mathrm{d}m}{\mathrm{d}t} = S(t, m)$$


Users provide a `MomentSource` callback when moment-source evolution is part of
the model:

```cpp
qmomkit::MomentSource source =
    [](double time, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> derivative(context.size(), 0.0);
        derivative[0] = 0.0;
        derivative[1] = -0.1 * context.value(1);
        return derivative;
    };
```

The returned vector must have the same size as the current stage vector.
QMOMKit does not know what each component means; the ordering must match the
moment ordering used by your QMOM setup. CQMOM phase-space evolution uses the
context-specific callback below so multidimensional source terms can follow the
configured CQMOM moment layout.

### QmomSourceContext

`QmomSourceContext` is a non-owning, stage-local view of the one-dimensional
QMOM moment state. Use:

- `context.value(order)` for a tracked current-stage value;
- `QmomSourceContext::moment(order)` for ordinary integer-order access. Tracked
  orders return their stage value directly, while untracked orders are evaluated from
  a quadrature reconstructed from that same SSP-RK stage when QMOM
  reconstruction metadata is present. Generic source-only updates reject
  untracked requests;
- `context.generalized_moment(order)` for a finite, non-negative scalar order
  from that same stage-local quadrature. Integer and fractional orders use the
  same validation and negative-abscissa domain rules as
  `evaluate_generalized_moment(...)`; this derived query is not added to the
  transported state.

Use the advanced `advance_moments_with_source(...)` operation when you want to
advance only the moment-source part of a workflow:

```cpp
const qmomkit::MomentSourceUpdateResult update =
    qmomkit::advance_moments_with_source(moments, time, time_step, source);

if (!update.ok())
{
    std::cerr << "Moment-source update failed: " << update.message << '\n';
    return 1;
}
```

Generic source-only updates do not carry reconstruction metadata. An untracked
`context.moment(...)` request is therefore rejected.

After a moment-source update, reconstruct quadrature before applying another
phase-space advection step:

```cpp
const qmomkit::QmomResult post_source =
    qmomkit::compute_qmom(qmomkit::MomentSet(update.moments), node_count);
```

## Boundary Contributions

Boundary contributions are optional additive raw-moment rates caused by
phase-space boundary exchange. Use `MomentBoundaryContribution` only when the
model describes nucleation, injection, removal, evaporation/condensation
exchange, or similar production/destruction $B$ through a phase-space boundary:
$$\frac{\mathrm{d}m}{\mathrm{d}t} = S(t, m) + B(t, m)$$



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

Use `advance_moments_with_source(...)` when QMOMKit should perform SSP-RK3
moment-source integration, or use the phase-space evolution helpers when it
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

When a boundary callback is supplied to a phase-space evolution helper, it is evaluated
with the moment source during each SSP-RK3 moment-source stage:

```text
effective_rate(t, m) = source(t, m) + boundary(t, m)
```

The boundary-aware phase-space evolution overloads take an explicit options argument to
keep overload resolution unambiguous.

## Physics Callbacks

Physics callbacks are the main place where an application extends QMOMKit.
QMOMKit handles moment storage, moment ordering, quadrature reconstruction,
state updates, and operator-split phase-space evolution. Your application supplies the
physical model: phase-space advection and moment-source or source-like
moment rates.

The callback extension points are independent. A model may use only
`PhaseSpaceVelocity`, only `MomentSource`, moment source plus boundary contribution,
velocity plus moment source, or all three categories. Choose the public routine that
matches the needed composition, and keep absent physics explicit as either no
callback in lower-level workflows or a zero-rate callback in full phase-space evolution
helpers whose signature requires one.

Most applications change three things:

- configuration: phase-space names, node counts, initial condition, output,
  and time-integration settings;
- physics callbacks: whichever velocity fields, moment-source terms, or
  boundary exchange contributions the model actually needs;
- phase-space evolution options: timestep controls, tolerances, and splitting scheme.

Most applications can treat internal moment storage, reconstruction logic,
quadrature layout, and CQMOM representation details as QMOMKit-owned parts of
the workflow.

### Phase-Space Advection

Advection physics moves quadrature nodes through phase space. Provide a
`PhaseSpaceVelocity` callback only when the physical model can be written as a
coordinate velocity:
$$\frac{\mathrm{d}x}{\mathrm{d}t} = \dot{\xi}(t, x)$$


The callback receives the current time and one terminal node coordinate vector
in phase-space dimension order. It returns one velocity component per
coordinate. This is the right place for drift, growth, relaxation, acceleration,
or any effect that is naturally represented by moving nodes.

### Moment Sources

Moment-source physics updates raw moments directly. Provide a `MomentSource`
only for internal or volumetric moment evolution and, when needed, a
`MomentBoundaryContribution` for phase-space boundary exchange. For CQMOM use
`CqmomMomentSource` and `CqmomMomentBoundaryContribution`. The combined model
is:
$$\frac{\mathrm{d}M_\alpha}{\mathrm{d}t} = S_\alpha(t, M) + B_\alpha(t, M)$$


This is the right place for effects that change moment values rather than
by phase-space advection moving existing nodes. The phase-space evolution helpers compose phase-space advection
and moment-source parts with the selected scheme. Boundary contributions enter
the existing moment-source update; they do not add a separate phase-space evolution operator.

### CqmomSourceContext

`CqmomSourceContext` describes the CQMOM moments that need moment-source or boundary
contribution rates during a moment-source update. It is passed to
`CqmomMomentSource` and `CqmomMomentBoundaryContribution` callbacks so each
callback can compute rates for the exact moment entries required by the
configured hierarchy.

Use:

- `context.time()` for the explicitly supplied current SSP-RK moment-source stage
  time; it is the same value as the callback's first argument. Callers that
  represent the initial instant must provide `0.0` deliberately;
- `context.size()` for the number of moment entries that need derivatives;
- `context.moment_count()` as a synonym for that source moment count;
- `context.dimension_count()` for the hierarchy dimension count;
- `context.indices()` and `context.values()` for the context-ordered
  multi-index and value arrays;
- `context.index(i)` for the `MultiIndex` of entry `i`;
- `context.value(i)` for the current raw moment value of entry `i`;
- `context.position(alpha)` or `context.value(alpha)` when a source term needs
  a particular retained moment by multi-index;
- `context.moment(alpha)` to obtain an ordinary integer-order raw moment by
  `MultiIndex`. Retained entries return their current SSP-RK stage value
  directly; an untracked entry is evaluated from a quadrature built from that
  same stage's moments. During intermediate SSP-RK stages, this is not the
  accepted pre-source quadrature;
- `context.generalized_moment(orders)` to obtain a derived non-negative
  real-order raw moment from that same source-stage quadrature. This uses the
  same lazy, context-local quadrature cache as untracked `moment(alpha)`
  requests and does not add a generalized moment to the transported state.
  The order vector has one `double` per phase-space dimension. A negative
  abscissa is valid only when its requested order is integral; a fractional
  order at a negative abscissa has no real-valued result and is rejected;
- `context.for_each_terminal_node(...)` for a quadrature-level source closure.
  It visits the global weight and phase-space-order coordinates of each
  terminal node from that same stage quadrature. The coordinate reference is
  read-only and valid only during its visitor callback;
- `context.hierarchy()` for the CQMOM hierarchy associated with the state.

Use `context.moment(...)` for moment-based source closures and
`context.for_each_terminal_node(...)` for quadrature-level source closures.
Source callbacks do not reconstruct CQMOM quadrature manually during
phase-space source evolution; both operations use the same context-local
stage reconstruction when needed. Generic source-only APIs do not carry
reconstruction metadata and reject untracked
`context.moment(...)` and terminal-node visitor requests.

The returned source values correspond to the context ordering. For entry `i`,
compute the derivative for `context.index(i)` from the current state and write
that derivative to return entry `i`. The multi-index uses phase-space
dimension order.

QMOMKit may request only the moments required by the configured hierarchy, so
CQMOM moment-source and boundary callbacks are most robust when they iterate over the
context entries instead of assuming a dense rectangular tensor.

For example, in a volume-based phase space `(V,T,u,v,w)`, a diameter factor
is proportional to `V^(1/3)`. A breakup closure for transported order `p` can
request its untracked fractional factor directly from the stage context:

```cpp
const double m = context.generalized_moment(
    {static_cast<double>(p) + 1.0 / 3.0,
     static_cast<double>(q), static_cast<double>(r),
     static_cast<double>(s), static_cast<double>(t)});
```

This is a derived quadrature query only; it does not change CQMOM moment
indices, hierarchy layout, or source-state storage.

CQMOM moment-source callbacks usually look like this:

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

## Bounded CQMOM Support During Phase-Space Evolution

When a `CqmomCase` or `CqmomState` carries coordinate support from
`phase_space.support_lower` and `phase_space.support_upper`, every CQMOM
reconstruction for that coordinate checks the returned abscissas against the
same interval. This includes initial inversion, post-source inversion, current
active-count preservation, and adaptive growth attempts.

Fixed branch counts are strict: if reconstruction returns nodes outside the
declared interval, the stage fails and the returned state rolls back to the
input state. Adaptive branch counts first try to preserve the current count,
then reduce to a lower count when needed. Growth trials are allowed only toward
the configured maximum; if a growth candidate returns out-of-interval nodes,
the current accepted count is retained.

QMOMKit verifies the returned nodes using the quadrature-bound tolerance.
Boundary nodes and nodes inside the interval within that tolerance are accepted
unchanged. Nodes outside the interval by more than the tolerance are rejected.
QMOMKit does not clip, project, or otherwise move nodes to fit the interval.

Successful post-source stages store canonical moments reconstructed from the
accepted quadrature and active counts. The returned state contains the accepted
bounded quadrature closure.

## Time Integration

QMOMKit uses two user-visible integration roles:

- node-advection integration advances quadrature coordinates through
  phase-space advection velocities using default adaptive RK23, adaptive RK45,
  or selected fixed-stage SSP-RK3;
- moment-source integration advances raw moments through moment-source or source-like
  moment updates with SSP-RK3.

Most users interact with these through `NodeAdvectionOptions`,
`MomentSourceUpdateOptions`, or the higher-level phase-space evolution options:

```cpp
qmomkit::QmomPhaseSpaceEvolutionOptions options;
options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
options.node_options.adaptive_options.relative_tolerance = 1.0e-7;
options.source_options.ssp_rk3_options.max_retries = 4;
```

The full phase-space evolution helpers support three schemes:

```text
strang
advsrc
srcadv
```

`strang` is the default. It applies a phase-space advection half-step, a moment-source full stage,
and a second phase-space advection half-step. This is second-order Strang splitting 
[[Strang, 1968)](https://www.doi.org/10.1137/0705041)] and is
the usual choice when both phase-space advection and source-like moment
evolution matter.

`advsrc` applies phase-space advection first, then the moment-source update. `srcadv` applies
the moment-source update first, then phase-space advection. These first-order Lie
splitting schemes [[Descombes et al., 2007](https://www.doi.org/10.1080/00207160701458716)] 
can be useful for comparisons, simple workflows, or cases
where your application deliberately wants a fixed splitting order.

## Full Phase-Space Evolution Step

Use `QmomCase::initial_state()` to construct an accepted one-dimensional state,
then use `advance_qmom_phase_space(...)` to advance that state by one phase-space evolution
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
        [](double, const qmomkit::QmomSourceContext& context)
        {
            return std::vector<double>(context.size(), 0.0);
        };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

    qmomkit::QmomState state =
        qmomkit::QmomCase{qmomkit::MomentSet(moments), node_count, {}}
            .initial_state();

    double time = 0.0;
    const double final_time = 0.1;
    const double time_step = 0.01;

    while (time < final_time)
    {
        const double step = std::min(time_step, final_time - time);
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_phase_space(
                state, time, step, velocity, source, options);

        if (!result.ok())
        {
            std::cerr << "Phase-space evolution failed: " << result.message << '\n';
            return 1;
        }

        state = result.state;
        moments = state.moments();
        time = result.final_time;
    }
}
```

Use `advance_qmom_phase_space(...)` when QMOMKit should own one configured
Strang or Lie splitting step. If the surrounding application owns stage
orchestration, coupling, or acceptance, call the one-dimensional state stages
directly:

```cpp
const qmomkit::QmomPhaseSpaceEvolutionResult first =
    qmomkit::advance_qmom_advection(
        state, time, 0.5 * time_step, velocity, options);

if (!first.ok())
{
    // Handle failure.
}

const qmomkit::QmomPhaseSpaceEvolutionResult source_stage =
    qmomkit::advance_qmom_moment_source(
        first.state, time + 0.5 * time_step, time_step, source, options);

if (!source_stage.ok())
{
    // Handle failure.
}

const qmomkit::QmomPhaseSpaceEvolutionResult second =
    qmomkit::advance_qmom_advection(
        source_stage.state, time + 0.5 * time_step, 0.5 * time_step, velocity, options);
```

Phase-space advection stages move the accepted quadrature topology and reconstruct the
stored moment vector without another QMOM inversion. Moment-source stages own
the moment-source update and the same fixed or adaptive post-source reconstruction
policy used by the full phase-space evolution helper. Application-managed workflows should
pass `QmomState` objects between stages instead of reconstructing quadrature or
mutating moments manually.

If the model also has phase-space boundary exchange, define a
`MomentBoundaryContribution` and call the boundary-aware overload with explicit
phase-space evolution options:

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

const qmomkit::QmomPhaseSpaceEvolutionResult result =
    qmomkit::advance_qmom_phase_space(
        state, time, step, velocity, source, boundary, options);
```

`QmomState` exposes read-only `moments()`, `quadrature()`,
`maximum_node_count()`, and `qmom_options()` accessors. Its
stored moments are reconstructed from its accepted quadrature, and
`QmomState::active_node_count()` is derived from the accepted quadrature size.
Users may copy and move accepted states for retry or in-process checkpointing,
but direct construction from independent moments and quadrature is intentionally
unavailable; no public restart/import factory is currently provided.
In adaptive QMOM, `QmomState::maximum_node_count()` remains the configured maximum. Node advection
preserves the active count and reconstructs the stored moment vector
from the advected quadrature without another Wheeler inversion. After an
accepted moment-source update, phase-space evolution first reconstructs at the current active count
and then tries at most one additional node when the configured maximum permits
it. Fixed-count QMOM remains strict at `node_count`.

The one-dimensional `QuadratureWriter` writes one row per active QMOM node:

```text
time node_index weight x0
```

Create it with the configured maximum node count. Fixed and adaptive QMOM use
the same writer, and adaptive output writes only active nodes without padding
inactive entries.

For CQMOM, use `advance_cqmom_phase_space(...)` with the `CqmomState` obtained
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

qmomkit::CqmomPhaseSpaceEvolutionOptions options;
options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

const qmomkit::CqmomMomentSource cqmom_source =
    [](double, const qmomkit::CqmomSourceContext& context)
    {
        return qmomkit::CqmomMomentRate{
            std::vector<double>(context.size(), 0.0)};
    };

const qmomkit::CqmomPhaseSpaceEvolutionResult result =
    qmomkit::advance_cqmom_phase_space(
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

const qmomkit::CqmomPhaseSpaceEvolutionResult result_with_boundary =
    qmomkit::advance_cqmom_phase_space(
        state, time, time_step, velocity, cqmom_source, cqmom_boundary, options);
```

Applications that own advection and moment-source stage orchestration can use the public CQMOM stage
helpers directly:

```cpp
const qmomkit::CqmomPhaseSpaceEvolutionResult first =
    qmomkit::advance_cqmom_advection(
        state, time, 0.5 * time_step, velocity, options);

if (!first.ok())
{
    // Handle failure.
}

const qmomkit::CqmomPhaseSpaceEvolutionResult source_stage =
    qmomkit::advance_cqmom_moment_source(
        first.state, time + 0.5 * time_step, time_step, cqmom_source, options);

if (!source_stage.ok())
{
    // Handle failure.
}

const qmomkit::CqmomPhaseSpaceEvolutionResult second =
    qmomkit::advance_cqmom_advection(
        source_stage.state, time + 0.5 * time_step,
        0.5 * time_step, velocity, options);
```

The boundary-aware moment-source evolution stage has the same shape as the high-level
boundary overload:

```cpp
const qmomkit::CqmomPhaseSpaceEvolutionResult source_with_boundary =
    qmomkit::advance_cqmom_moment_source(
        first.state, time + 0.5 * time_step, time_step,
        cqmom_source, cqmom_boundary, options);
```

Phase-space advection CQMOM stages preserve terminal paths and active branch counts, move
the accepted terminal quadrature nodes, and reconstruct moments from that
advected quadrature without a post-stage CQMOM inversion. Moment-source stages
advance the CQMOM phase-space evolution moment layout, then apply the same recursive fixed
or adaptive reconstruction policy used by `advance_cqmom_phase_space(...)`.
Moment-source stages ignore node-advection options, phase-space advection stages ignore moment-source
integration options, and standalone stage helpers do not use the splitting
scheme option. On failure, the returned state is the stage input state; the
returned quadrature, when present, is diagnostic unless `result.ok()` is true.

Here `compute_cqmom(...)` computes the initial quadrature, while
`CqmomCase::initial_state()` creates the runtime state that phase-space evolution advances.

Adaptive branch counts configured on the `CqmomCase` are preserved by
`CqmomState` phase-space evolution. In adaptive mode, the hierarchy node counts remain the
maximum branch counts; the state also retains the current branch-specific
active counts associated with its ragged quadrature. Node advection preserves
terminal paths and therefore preserves active counts. Moments are
reconstructed by summation from the advected quadrature; this stage does not
grow active branch counts.

At reconstruction boundaries, existing branches use their current active count
as the ceiling, so active counts may stay the same or decrease but cannot grow.
After an accepted moment-source update, post-source inversion first
tests the current active count; if that count remains valid and the configured
maximum is larger, it tries exactly one additional node. Growth is accepted only
when that trial returns the requested `Nactive + 1` nodes. If a parent topology
change creates a genuinely new descendant branch, the descendant begins from
its own configured maximum node count, and adaptive Wheeler determines its
initial active count.

Adaptive `CqmomState` construction seeds active counts from the initial CQMOM
inversion; if that inversion fails, state construction fails instead of creating
a state with inconsistent adaptive metadata. Failed CQMOM phase-space evolution through
`CqmomState` returns the previous state metadata so callers can retry from the
last accepted step.

`CqmomPhaseSpaceEvolutionOptions` contains only phase-space evolution stage policies: node advection,
moment-source integration, and the phase-space evolution scheme.
CQMOM inversion policies such as Wheeler options, adaptive branch counts,
diagnostic level, and Vandermonde solver belong to `CqmomCase` and are preserved
by `CqmomState`. Successful CQMOM phase-space evolution through `CqmomState` persists active
adaptive topology for the next call. Explicit dense and lean case types remain
available for advanced static inversion and specialized initialization, but
normal CQMOM phase-space evolution uses `CqmomCase::initial_state()` and
`advance_cqmom_phase_space(...)` with the storage-independent `CqmomState` API.

For storage-independent CQMOM physics, use `CqmomMomentSource` and
`CqmomMomentBoundaryContribution`. `CqmomSourceContext` exposes the CQMOM moment
layout selected for phase-space evolution, with the same indices, values, ordering, hierarchy,
and lookup behavior for `CqmomState` with dense or lean moment storage.

`result.ok()` identifies whether returned moments or quadrature describe an
accepted phase-space evolution step. On failure, `result.message` and `result.diagnostics`
indicate whether the inversion, node advection, moment-source update, or
moment reconstruction failed.

## Configuration Workflow

Configuration files can select initialization, quadrature settings, phase-space evolution
scheme, and output options for drivers that read those sections. The phase-space evolution
scheme is set in:

```ini
[time_integration]
scheme = strang
```

Accepted values are:

```ini
strang
advsrc
srcadv
```

Node-advection integration is configured separately because it controls the
phase-space advection trajectory integrator, not the operator-splitting scheme:

```ini
[time_integration.node_advection]
scheme = rk23
```

Accepted values are:

```ini
rk23
rk45
ssp_rk3
```

`rk23` is the default adaptive Bogacki-Shampine 3(2) method. `rk45` is the
adaptive Dormand-Prince RK5(4)7M pair from Table 2 of Dormand and Prince
(1980). `ssp_rk3` remains a separate fixed-stage method. The embedded methods
use their higher-order state for accepted steps and reuse FSAL derivatives.

Adaptive controls are configured in a nested section:

```ini
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
estimate it from the state, derivatives, or tolerances. Startup summaries in the phase-space evolution examples report the effective
node-advection scheme and controls. For fixed-step schemes, the summary reports
only the selected scheme. Moment-source integration is unchanged.

Output can be requested with:

```ini
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
```

Configuration files do not define the velocity callback, source callback,
boundary contribution callback, or outer time loop. Those remain application
code. Drivers can read the scheme with the public configuration API and assign
it to `QmomPhaseSpaceEvolutionOptions` or `CqmomPhaseSpaceEvolutionOptions`. See
`configuration.md` for the configuration syntax and the example files under:

```text
examples/qmom_1d_breakage_advection/
examples/cqmom_2d_minimal/
```

## Batch Convenience Wrappers

When an application owns several independent accepted states with shared
phase-space evolution controls, it can call the `*_batch(...)` accepted-state phase-space evolution
wrappers. QMOM and CQMOM provide matching phase-space advection stage, moment-source evolution stage,
and complete phase-space evolution batch families, including the corresponding
boundary-aware overloads.

A batch call is a serial convenience wrapper around the corresponding scalar
call. It processes input positions in increasing order and returns one result
per input in the same order. Each entry retains its own adaptive substeps,
retries, active node or branch topology, diagnostics, and phase-space evolution status. An
empty input returns an empty result without invoking callbacks or validating the
shared arguments.

Use ordinary batch functions such as `advance_qmom_phase_space_batch(...)` when
one scalar callback applies to every entry. When external data vary by input
position, use the batch-index-aware `*_batch_with_index(...)` functions with
the `Batch*` callback types. Their first `batch_index` argument is only the
zero-based position in the supplied state collection; it is not a cell, mesh,
or spatial identifier. The distinct function names select the callback family,
so generic and dual-signature callables remain unambiguous.

```cpp
const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> results =
    qmomkit::advance_qmom_phase_space_batch(
        states, time, time_step, velocity, source, options);

const qmomkit::BatchPhaseSpaceVelocity indexed_velocity =
    [](std::size_t, double, const std::vector<double>& coordinates)
    {
        return std::vector<double>(coordinates.size(), 0.0);
    };

const qmomkit::BatchMomentSource indexed_source =
    [&local_rates](std::size_t batch_index, double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rates[i] = local_rates[batch_index] * context.value(i);
        }
        return rates;
    };

const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> indexed_results =
    qmomkit::advance_qmom_phase_space_batch_with_index(
        states, time, time_step, indexed_velocity,
        indexed_source, options);
```

Time, time step, and phase-space evolution options are shared by a batch; each
`QmomState` supplies its own reconstruction metadata. Applications needing heterogeneous controls should make explicit
scalar calls. An unsuccessful scalar result is retained at its input position
and does not stop later entries. An exception that escapes a scalar call also
escapes the batch call unchanged and stops processing.

Batch wrappers do not provide threading, MPI, scheduling, or spatial
organization. Calling scalar APIs inside an application-managed parallel loop
remains fully supported.

## Practical Guidance

Coordinate order links the phase-space evolution callback to the configured model. A
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

Moment sources update raw moments. Quadrature nodes are rebuilt from those
moments before another node-advection step.

Boundary contributions are added to the moment-source rate and integrated by
the same SSP-RK3 moment-source update path. Moment-source retries, acceptability predicates,
and post-source inversion behavior are shared.

The evolved moment set still needs to support the requested quadrature.
Phase-space evolution can fail if updated moments no longer support the requested QMOM or
CQMOM reconstruction.

Callbacks are application code. QMOMKit does not provide meshes, boundary
conditions, boundary geometry, physical source models, or units. The library
advances the moment and quadrature data provided by the surrounding workflow.
