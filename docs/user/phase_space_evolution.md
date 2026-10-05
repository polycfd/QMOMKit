# Phase-space evolution and time integration

QMOM and CQMOM turn raw moments into quadrature nodes and weights. This guide
covers what to do with that quadrature next: advancing it in time and getting
updated moments back for the following step or for output.

QMOMKit supplies the numerical machinery — advancing nodes, integrating
sources, reconstructing moments, and composing those into one step. Your
application supplies the physics and the surrounding timestep loop, and decides
what the phase-space coordinates mean.

The split between the two kinds of physics is the useful part. Smooth effects
such as growth, drift, acceleration, or relaxation move quadrature nodes
through phase space. Effects that create or destroy population — breakage,
aggregation, chemistry, or rates computed elsewhere — change the moments
directly instead. Exchange across a phase-space boundary, such as nucleation,
injection, removal, or evaporation, is a third contribution added to the source
rates.

Each kind is a separate callback, and you provide only the ones your model
needs. Examples showing several together are illustrative combinations, not a
requirement. [Callbacks](callbacks.md) gives their signatures and rules.

## What is on this page

| Section | Covers |
| --- | --- |
| [Workflow](#phase-space-evolution-workflow) | How a step is composed, and which routine to call |
| [Physics callbacks](#physics-callbacks) | Which physics belongs in which callback |
| [Node advection](#node-advection) | Moving quadrature nodes through a velocity field |
| [Moment source terms](#moment-source-terms) | Rates applied directly to raw moments |
| [Boundary contributions](#boundary-contributions) | Exchange across a phase-space boundary |
| [Time integration](#time-integration) | Splitting schemes and node-advection schemes |
| [Full step](#full-phase-space-evolution-step) | Complete QMOM and CQMOM steps, and driving the stages yourself |
| [Configuration](#configuration-workflow) | Loading these settings from a `.qmom` file |
| [Batch wrappers](#batch-convenience-wrappers) | Advancing many states with shared settings |
| [Practical guidance](#practical-guidance) | Ordering, sizes, and what to check |

For callback signatures and examples see [Callbacks](callbacks.md); when a step
fails, see [Troubleshooting](troubleshooting.md).

## Phase-space evolution workflow

A step alternates between the two representations. Quadrature is the convenient
form for moving nodes; raw moments are the convenient form for source terms and
for the next reconstruction.

```text
initial moments/quadrature
        ↓
node advection, if present
        ↓
moment sources and boundary rates, if present
        ↓
updated quadrature/moments
```

There are two operations to compose:

- **node advection** moves quadrature abscissas through a velocity field;
- **the moment-source update** applies source and boundary rates to the raw
  moments.

Most applications want both composed for them, which is what
`advance_qmom_phase_space(...)` and `advance_cqmom_phase_space(...)` do. Where
one of those helpers asks for physics your model does not have, pass a callback
that returns zero rates.

If you would rather drive the pieces yourself, call `advance_nodes(...)` for
advection alone — reconstructing moments afterwards as needed — or
`advance_qmom_moments_with_source(...)` for sources alone. That is the right choice
when your application needs to interleave its own work between the two.

## Physics callbacks

This section covers *which* physics belongs in which callback. For the
signatures, worked examples, and the rules every callback must satisfy, see
[Callbacks](callbacks.md).

Physics callbacks are the main place where an application extends QMOMKit.
QMOMKit handles moment storage, moment ordering, quadrature reconstruction,
state updates, and the operator splitting. Your application supplies the
physical model: phase-space advection and moment-source or source-like
moment rates.

Provide a `PhaseSpaceVelocity` when the model can be written as a velocity in
the phase-space coordinates:

$$\frac{\mathrm{d}\xi}{\mathrm{d}t} = \dot{\xi}(t, x)$$

This is the right place for drift, growth, relaxation, acceleration, or any
effect naturally represented by moving existing nodes.

Provide a `QmomMomentSource` — `CqmomMomentSource` in several dimensions — for
effects that change moment values instead, together with an optional boundary
contribution for exchange across a phase-space boundary:

$$\frac{\mathrm{d}M_\alpha}{\mathrm{d}t} = S_\alpha(t, M) + B_\alpha(t, M)$$

Boundary contributions are added to the source rates rather than forming a
separate operator.

The extension points are independent: a model may use only a velocity, only a
source, a source with a boundary contribution, a velocity with a source, or all
of them. Choose the routine matching the composition you need, and where a
composed helper asks for physics your model does not have, pass a callback
returning zero rates.

Most applications change three things:

- configuration: phase-space names, node counts, initial condition, output,
  and time-integration settings;
- physics callbacks: whichever velocity fields, moment-source terms, or
  boundary exchange contributions the model actually needs;
- step options: timestep controls, tolerances, and splitting scheme.

Most applications can treat internal moment storage, reconstruction logic,
quadrature layout, and CQMOM representation details as QMOMKit-owned parts of
the workflow.

## Node advection

Node advection moves quadrature node coordinates through a phase-space
velocity field. Use it when the model says that existing
representative states move through phase space. In one-dimensional QMOM, this
advances the abscissas of the QMOM nodes. In multidimensional CQMOM, this
advances the terminal node coordinates in phase-space dimension order.

Users provide a `PhaseSpaceVelocity` callback when phase-space advection
motion is part of the model. The callback is evaluated once per right-hand-side
evaluation for all quadrature nodes together:

```cpp
// One-dimensional phase space, so `dimensions` is 1 and each node owns
// exactly one entry.
qmomkit::PhaseSpaceVelocity velocity =
    [](double time, const std::vector<double>& coordinates,
       std::vector<double>& velocities, std::size_t)
    {
        for (std::size_t node = 0; node < coordinates.size(); ++node)
        {
            const double x = coordinates[node];
            velocities[node] = 1.0 - x;
        }
    };
```

In more than one dimension every phase-space dimension of every node must be
written, for example:

```cpp
qmomkit::PhaseSpaceVelocity velocity =
    [](double time, const std::vector<double>& coordinates,
       std::vector<double>& velocities, std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            for (std::size_t d = 0; d < dimensions; ++d)
            {
                velocities[node * dimensions + d] = -coordinates[node * dimensions + d];
            }
        }
    };
```

`coordinates` is packed node-major: node `i` occupies
`[i * dimensions, (i + 1) * dimensions)`, and within a node the values are
ordered by phase-space dimension. `velocities` arrives with the same size and
must be filled with $\mathrm{d}\xi/\mathrm{d}t$ in the same layout. Do not
resize it; a size mismatch is reported as a velocity failure.

Fill every entry. Leaving one unwritten is a mistake rather than a request for
zero velocity, and the affected node would simply not move. A debug build
catches it and names the node and phase-space dimension; an optimised build
cannot, so it is worth developing callbacks against a debug build.

QMOMKit does not assign a physical meaning to those coordinates; `x`,
`diameter`, `velocity`, or any other interpretation comes from your application
and must match the phase-space order used elsewhere. If a component is
non-finite, the diagnostic includes the evaluation time, zero-based phase-space
dimension, offending value, the coordinate vector of the offending node, and
that node's complete velocity vector.

Evaluating all nodes together keeps the cost of one right-hand-side evaluation
independent of the terminal node count, which is the product of the
per-dimension node counts and therefore grows quickly with phase-space
dimension.

### Adapting a per-node callback

When the physical model is naturally written one node at a time, write a
`PerNodePhaseSpaceVelocity` and adapt it with `make_phase_space_velocity(...)`:

```cpp
qmomkit::PhaseSpaceVelocity velocity =
    qmomkit::make_phase_space_velocity(
        [](double time, const std::vector<double>& coordinate)
        {
            const double x = coordinate[0];
            return std::vector<double>{1.0 - x};
        });
```

The per-node callback receives one node coordinate vector, ordered by
phase-space dimension, and must return $\mathrm{d}\xi/\mathrm{d}t$ with the
same number of components. The adapted callback evaluates it once per node in
ascending node order and preserves its per-node validation, so both forms
produce identical results. It therefore also costs one call and one returned
vector per node per right-hand-side evaluation; prefer writing a
`PhaseSpaceVelocity` directly when advancing many nodes.

`PhaseSpaceVelocity` is the single velocity form accepted by
`advance_nodes(...)`, by the advection-only stages, and by the full
`advance_qmom_phase_space(...)` and `advance_cqmom_phase_space(...)`
operations. The saving from writing it directly grows with phase-space
dimension, because the node count does.

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

The full QMOM and CQMOM helpers described below keep accepted state
together with the reconstruction settings needed for the next step.

`NodeAdvectionOptions::scheme` selects how node trajectories are integrated:

- `rk23`: the adaptive Bogacki-Shampine 3(2) scheme [[Bogacki & Shampine, 1989](https://www.doi.org/10.1016/0893-9659(89)90079-7)];
- `rk45`: the adaptive Dormand-Prince RK5(4)7M scheme [[Dormand & Prince, 1980](https://www.doi.org/10.1016/0771-050X(80)90013-3)];
- `ssp_rk3`: the fixed-stage SSP-RK3 scheme [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)].

All three choices integrate the same advection quadrature-node trajectory
equation. SSP-RK3 node advection does not retry, reduce the timestep, or
perform realizability checks. Its accuracy is controlled by the timestep chosen
by the caller. Realizability remains a property of reconstructed moments and
later quadrature inversions, not of the coordinate trajectory step itself.

Configuration-driven drivers can load the complete node-advection controls
with `qmomkit::config::load_node_advection_options(...)`. For RK23 and RK45,
`NodeAdvectionOptions::adaptive_options` carries scalar relative and absolute
tolerances, an initial step that defaults to `1.0e-6 s`, maximum step size,
minimum step size, and maximum internal accepted steps. The initial step may
be overridden but must be finite and strictly positive. QMOMKit does not
estimate it from the state, derivatives, or tolerances; accepted states use
the higher-order embedded solution and the methods reuse FSAL derivatives
internally.

## Moment source terms

Moment sources represent internal or volumetric changes that act directly on
raw moments. They are integrated with an explicit third-order strong-stability-preserving
Runge-Kutta (SSP-RK3) scheme [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)].
Use moment sources when the model gives rates for moments, or when an event
changes moments without being well described as smooth coordinate
motion. Examples include aggregation, fragmentation, relaxation, chemistry, or
any other model that can be written as:
$$\frac{\mathrm{d}M}{\mathrm{d}t} = S(t, M)$$


Users provide a `QmomMomentSource` callback when moment-source evolution is part of
the model:

```cpp
qmomkit::QmomMomentSource source =
    [](double time, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        rate.values[0] = 0.0;
        rate.values[1] = -0.1 * context.value(1);
        return rate;
    };
```

The returned rate must have the same size as the current stage vector.
QMOMKit does not know what each component means; the ordering must match the
moment ordering used by your QMOM setup. CQMOM uses the
context-specific callback below so multidimensional source terms can follow the
configured CQMOM moment layout.

### QmomSourceContext

`QmomSourceContext` gives a source callback the moment values it needs. It is
valid only for the duration of the call, and the values it reports are those of
the current Runge-Kutta stage rather than the accepted state. Use:

- `context.value(order)` for the current value of a transported moment;
- `context.moment(order)` for any integer order. Transported orders return their
  stage value directly; an order you do not transport is evaluated from a
  quadrature built from that same stage;
- `context.generalized_moment(order)` for a real, non-negative order, evaluated
  from that same quadrature. It obeys the same rules as
  `evaluate_generalized_moment(...)`, and it is a query only: asking for it does
  not add that order to the moments you evolve.

Requests for an order you do not transport need a quadrature, so they work only
when you evolve a `QmomState`. The generic vector-based
`advance_qmom_moments_with_source(...)` shown below has no reconstruction settings to
work from and rejects them.

Use the advanced `advance_qmom_moments_with_source(...)` operation when you want to
advance only the moment-source part of a workflow:

```cpp
const qmomkit::QmomMomentSourceUpdateResult update =
    qmomkit::advance_qmom_moments_with_source(moments, time, time_step, source);

if (!update.ok())
{
    std::cerr << "Moment-source update failed: " << update.message << '\n';
    return 1;
}
```

After a moment-source update, reconstruct quadrature before applying another
phase-space advection step:

```cpp
const qmomkit::QmomResult post_source =
    qmomkit::compute_qmom(qmomkit::MomentSet(update.moments), node_count);
```

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
- `context.index(i)` for a `MultiIndexView` of entry `i`, which indexes and
  iterates like a `MultiIndex` without copying it;
- `context.value(i)` for the current raw moment value of entry `i`;
- `context.position(alpha)` or `context.value(alpha)` when a source term needs
  a particular retained moment by multi-index;
- `context.moment(alpha)` to obtain an ordinary integer-order raw moment by
  `MultiIndex`. Retained entries return their current SSP-RK stage value
  directly; an untracked entry is evaluated from a quadrature built from that
  same stage's moments. During intermediate SSP-RK stages, this is not the
  accepted pre-source quadrature;
- `context.generalized_moment(orders)` for a real, non-negative order in each
  phase-space dimension, evaluated from that same stage quadrature. Asking for
  it does not add that order to the moments you evolve. A negative abscissa is
  valid only at an integer order; a fractional order there has no real value and
  is rejected;
- `context.for_each_terminal_node(...)` to visit each terminal node's weight and
  coordinates directly, when a rate is easier to write over nodes than over
  moments. The coordinates it hands you are read-only and valid only inside the
  visit;
- `context.hierarchy()` for the CQMOM hierarchy associated with the state.

Write the rate over moments with `context.moment(...)`, or over nodes with
`context.for_each_terminal_node(...)`; there is no need to reconstruct a
quadrature yourself, as the context does it for you when one is needed. As in
the one-dimensional case, both require a `CqmomState` to evolve: the generic
vector-based source update rejects them.

The returned source values correspond to the context ordering. For entry `i`,
compute the derivative for `context.index(i)` from the current state and write
that derivative to return entry `i`. The multi-index uses phase-space
dimension order.

QMOMKit may request only the moments required by the configured hierarchy, so
CQMOM moment-source and boundary callbacks are most robust when they iterate
over the context entries instead of assuming a full rectangular tensor.

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
    [](double, const qmomkit::CqmomSourceContext& context)
    {
        const qmomkit::PhaseSpaceHierarchy& hierarchy = context.hierarchy();
        std::vector<double> derivatives(context.size(), 0.0);

        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndexView order = context.index(entry);
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

        return qmomkit::CqmomMomentRate{derivatives};
    };
```

This pattern is dimension-independent: it uses the hierarchy dimension count
instead of assuming a two-dimensional phase space.

## Boundary contributions

Boundary contributions are optional additive raw-moment rates caused by
phase-space boundary exchange. Use `QmomMomentBoundaryContribution` only when the
model describes nucleation, injection, removal, evaporation/condensation
exchange, or similar production/destruction $B$ through a phase-space boundary:
$$\frac{\mathrm{d}M}{\mathrm{d}t} = S(t, M) + B(t, M)$$

For one-dimensional QMOM:

```cpp
qmomkit::QmomMomentBoundaryContribution boundary =
    [](double time, const qmomkit::QmomSourceContext& context)
    {
        constexpr double rate = 0.02;
        constexpr double boundary_coordinate = 0.1;

        std::vector<double> values(context.size(), 0.0);
        double power = 1.0;
        for (double& value : values)
        {
            value = rate * power;
            power *= boundary_coordinate;
        }

        (void)time;
        return qmomkit::QmomMomentRate{values};
    };
```

Use `advance_qmom_moments_with_source(...)` when QMOMKit should perform SSP-RK3
moment-source integration, or use the composed helpers when it
should also own the configured operator-splitting sequence.

For CQMOM, use `CqmomMomentBoundaryContribution`. The two forms are shaped
alike: each boundary callback receives the same context as its own source
callback and returns one additive rate per context entry:

```cpp
qmomkit::CqmomMomentBoundaryContribution boundary =
    [](double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> values(context.size(), 0.0);

        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndexView order = context.index(entry);
            values[entry] = /* boundary exchange rate for M_order */;
        }

        (void)time;
        return qmomkit::CqmomMomentRate{values};
    };
```

Boundary contributions are not a generic PDE boundary-condition framework.
They do not define mesh faces, normals, reflection, clipping, or boundary
geometry. They are runtime physics callbacks, just like velocities and moment
sources, and are not stored in `QmomCase`, `CqmomCase`, `CqmomPreparedCase`, or
`CqmomState`.

When a boundary callback is supplied to one of these helpers, it is evaluated
with the moment source during each SSP-RK3 moment-source stage:

```text
effective_rate(t, M) = source(t, M) + boundary(t, M)
```

The boundary-aware overloads take an explicit options argument to
keep overload resolution unambiguous.

## Time integration

QMOMKit uses two user-visible integration roles:

- node-advection integration advances quadrature coordinates through
  phase-space advection velocities using default adaptive RK23, adaptive RK45,
  or selected fixed-stage SSP-RK3;
- moment-source integration advances raw moments through moment-source or source-like
  moment updates with SSP-RK3.

Most users interact with these through `NodeAdvectionOptions`,
`MomentSourceUpdateOptions`, or the higher-level step options:

```cpp
qmomkit::QmomPhaseSpaceEvolutionOptions options;
options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
options.node_options.adaptive_options.relative_tolerance = 1.0e-7;
options.source_options.ssp_rk3_options.max_retries = 4;
```

Phase-space evolution supports three splitting schemes:

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

## Full phase-space evolution step

### One-dimensional QMOM

Use `QmomCase::initial_state()` to construct an accepted one-dimensional state,
then use `advance_qmom_phase_space(...)` to advance that state by one
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
        [](double, const std::vector<double>& coordinates,
           std::vector<double>& velocities, std::size_t)
        {
            for (std::size_t node = 0; node < coordinates.size(); ++node)
            {
                velocities[node] = 1.0 - coordinates[node];
            }
        };

    const qmomkit::QmomMomentSource source =
        [](double, const qmomkit::QmomSourceContext& context)
        {
            return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
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
        moments = state.moment_values();
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

Phase-space advection stages move the accepted quadrature topology and
reconstruct the stored moment vector without another QMOM inversion.
Moment-source stages own the moment-source update and the same fixed or
adaptive post-source reconstruction policy used by the composed helper.
Application-managed workflows should pass `QmomState` objects between stages
instead of reconstructing quadrature or mutating moments manually.

If the model also has phase-space boundary exchange, define a
`QmomMomentBoundaryContribution` and call the boundary-aware overload with explicit
phase-space evolution options:

```cpp
const qmomkit::QmomMomentBoundaryContribution boundary =
    [](double, const qmomkit::QmomSourceContext& context)
    {
        constexpr double injection_rate = 0.02;
        constexpr double injection_coordinate = 0.1;

        std::vector<double> rate(context.size(), 0.0);
        double power = 1.0;
        for (double& value : rate)
        {
            value = injection_rate * power;
            power *= injection_coordinate;
        }
        return qmomkit::QmomMomentRate{rate};
    };

const qmomkit::QmomPhaseSpaceEvolutionResult result =
    qmomkit::advance_qmom_phase_space(
        state, time, step, velocity, source, boundary, options);
```

### Inspecting and writing the state

`QmomState` exposes read-only `moment_values()`, `quadrature()`,
`maximum_node_count()`, and `qmom_options()` accessors. Its
stored moments are reconstructed from its accepted quadrature, and
`QmomState::active_node_count()` is derived from the accepted quadrature size.
Users may copy and move accepted states for retry or in-process checkpointing,
but direct construction from independent moments and quadrature is intentionally
unavailable; no public restart/import factory is currently provided.
In adaptive QMOM, `QmomState::maximum_node_count()` remains the configured maximum. Node advection
preserves the active count and reconstructs the stored moment vector
from the advected quadrature without another Wheeler inversion. After an
accepted moment-source update, QMOMKit first reconstructs at the current active count
and then tries at most one additional node when the configured maximum permits
it. Fixed-count QMOM remains strict at `node_count`.

The one-dimensional `QuadratureWriter` writes one row per active QMOM node:

```text
time node_index weight x0
```

Create it with the configured maximum node count. Fixed and adaptive QMOM use
the same writer, and adaptive output writes only active nodes without padding
inactive entries.

### Multidimensional CQMOM

For CQMOM, use `advance_cqmom_phase_space(...)` with the `CqmomState` obtained
from the loaded CQMOM case when QMOMKit should own one complete configured
splitting step:

```cpp
const auto configured = qmomkit::config::load_cqmom_prepared("configuration.qmom");
if (!configured)
{
    // Report configured.error(), then stop.
}

const qmomkit::CqmomPreparedCase cqmom_prepared_case = configured.value();

const qmomkit::CqmomResult initial_quadrature = qmomkit::compute_cqmom(cqmom_prepared_case);
if (!initial_quadrature.ok())
{
    // Report initial_quadrature.message, then stop.
}

qmomkit::CqmomState state = cqmom_prepared_case.initial_state();

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
            const qmomkit::MultiIndexView order = context.index(entry);
            rate[entry] = /* phase-space boundary exchange rate for M_order */;
        }
        return qmomkit::CqmomMomentRate{rate};
    };

const qmomkit::CqmomPhaseSpaceEvolutionResult result_with_boundary =
    qmomkit::advance_cqmom_phase_space(
        state, time, time_step, velocity, cqmom_source, cqmom_boundary, options);
```

### Driving the stages yourself

Applications that own advection and moment-source stage orchestration can use
the public CQMOM stage helpers directly:

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

The boundary-aware moment-source evolution stage has the same shape as the
high-level boundary overload:

```cpp
const qmomkit::CqmomPhaseSpaceEvolutionResult source_with_boundary =
    qmomkit::advance_cqmom_moment_source(
        first.state, time + 0.5 * time_step, time_step,
        cqmom_source, cqmom_boundary, options);
```

Phase-space advection CQMOM stages preserve terminal paths and active branch
counts, move the accepted terminal quadrature nodes, and reconstruct moments
from that advected quadrature without a post-stage CQMOM inversion.
Moment-source stages advance the CQMOM phase-space evolution moment layout,
then apply the same recursive fixed or adaptive reconstruction policy used by
`advance_cqmom_phase_space(...)`. Moment-source stages ignore node-advection
options, phase-space advection stages ignore moment-source integration options,
and standalone stage helpers do not use the splitting scheme option. On
success, the accepted quadrature is `result.state.quadrature()` and
`result.diagnostic_quadrature` is empty. On failure, the returned state is the
stage input state, and `result.diagnostic_quadrature`, when present, holds the
quadrature the failing stage produced or started from, for diagnosis only.

Here `compute_cqmom(...)` computes the initial quadrature, while
`CqmomPreparedCase::initial_state()` creates the runtime state that phase-space
evolution advances.

### Adaptive node counts

In adaptive mode the node counts you configure are ceilings rather than fixed
sizes: each branch uses as many nodes as its local moments support, and the
state remembers how many that currently is. Moving nodes does not change those
counts, since advection relocates the existing nodes rather than adding any. A
branch can only gain a node at a reconstruction, described next.

Node counts shrink readily and grow cautiously. A reconstruction during a step
may reduce a branch's count if its moments no longer support the current one,
but never raises it. Growth happens only after an accepted moment-source update,
and only one node at a time: QMOMKit retries with one extra node and keeps the
result only if that larger rule actually succeeds. This asymmetry is deliberate,
since a branch that has just lost resolution is rarely ready to regain it in the
same step. A branch that appears for the first time starts from its configured
maximum and settles wherever its own moments allow.

An adaptive `CqmomState` takes its starting active counts from the initial
inversion, so if that inversion fails you get no state at all rather than an
unusable one. A failed step likewise returns the previous state, so you can
retry from the last accepted one.

### Where each option belongs

The two option groups divide cleanly, which is worth knowing when deciding
where a setting belongs. `CqmomPhaseSpaceEvolutionOptions` holds only what a
step does: node advection, moment-source integration, and the splitting scheme.
How the quadrature is reconstructed — Wheeler options, adaptive branch counts,
diagnostic level, solver choice — belongs to `CqmomPreparedCase`, and
`CqmomState` carries it forward from step to step along with the adaptive node
counts it has settled on.

Write CQMOM physics against `CqmomMomentSource` and
`CqmomMomentBoundaryContribution`. The context they receive presents the same
moments in the same order however the state stores them internally, so a
callback written once does not depend on that.

### Checking the result

`result.ok()` identifies whether returned moments or quadrature describe an
accepted phase-space evolution step. On failure, `result.message` and `result.diagnostics`
indicate whether the inversion, node advection, moment-source update, or
moment reconstruction failed.

`result.final_time` is always the time of `result.state`. On success it is
`time + time_step`. On failure the returned state is the input state, so
`result.final_time` is the input `time`, whichever stage failed.

## Configuration workflow

Configuration files can select initialization, quadrature settings, phase-space
evolution scheme, and output options for drivers that read those sections. The
phase-space evolution scheme is set in:

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

## Batch convenience wrappers

When an application owns several independent accepted states with shared
phase-space evolution controls, it can call the `*_batch(...)` accepted-state
phase-space evolution wrappers. QMOM and CQMOM provide matching phase-space
advection stage, moment-source evolution stage, and complete phase-space
evolution batch families, including the corresponding boundary-aware overloads.

A batch call is a serial convenience wrapper around the corresponding scalar
call. It processes input positions in increasing order and returns one result
per input in the same order. Each entry retains its own adaptive substeps,
retries, active node or branch topology, diagnostics, and phase-space evolution
status. An empty input returns an empty result without invoking callbacks or
validating the shared arguments.

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
    [](std::size_t, double, const std::vector<double>&,
       std::vector<double>& velocities, std::size_t)
    {
        std::fill(velocities.begin(), velocities.end(), 0.0);
    };

const qmomkit::BatchQmomMomentSource indexed_source =
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

A batch shares one time, time step, and set of options, while each state keeps
its own reconstruction settings. If entries need different controls, call them
individually instead. A failed entry keeps its place in the returned vector and
does not stop the others, but an exception escaping one entry stops the batch.

Batch wrappers do not provide threading, MPI, scheduling, or spatial
organization. Calling scalar APIs inside an application-managed parallel loop
remains fully supported.

[`examples/cqmom_2d_batch_host_data/`](../../examples/cqmom_2d_batch_host_data/README.md)
is a complete program using `advance_cqmom_phase_space_batch_with_index(...)`
with time-varying host data indexed like the state vector.

## Practical guidance

Coordinate order links each callback to the configured model. A
`PhaseSpaceVelocity` callback receives coordinates in phase-space dimension
order within each node, the same order used by CQMOM moments, node counts, and
output.

Each moment callback must return one rate per entry of the context it receives.
A `QmomMomentSource` or `QmomMomentBoundaryContribution` callback returns one
rate per `QmomSourceContext` entry; a `CqmomMomentSource` or
`CqmomMomentBoundaryContribution` callback returns one rate per
`CqmomSourceContext` entry.

For CQMOM source and boundary rates, the context order defines the return
order. For entry `i`, compute the rate for `context.index(i)` from the current
value `context.value(i)`, then store it in return entry `i`.

Moment sources update raw moments. Quadrature nodes are rebuilt from those
moments before another node-advection step.

Boundary contributions are added to the moment-source rate and integrated by
the same SSP-RK3 moment-source update path. Moment-source retries,
acceptability predicates, and post-source inversion behavior are shared.

The evolved moment set still needs to support the requested quadrature.
Phase-space evolution can fail if updated moments no longer support the
requested QMOM or CQMOM reconstruction.

Callbacks are application code. QMOMKit does not provide meshes, boundary
conditions, boundary geometry, physical source models, or units. The library
advances the moment and quadrature data provided by the surrounding workflow.
