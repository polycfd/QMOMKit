# Callbacks

Callbacks are where your application meets QMOMKit. The library owns moment
storage and ordering, quadrature reconstruction, state updates, and the
operator-split phase-space evolution. Your callbacks supply the physics: how
fast quadrature nodes move through phase space, and how raw moments change for
reasons other than that motion.

Every callback is a `std::function` alias, so a lambda, a free function, or a
functor all bind to one. Nothing in QMOMKit assigns physical meaning or units to
a coordinate or a moment; that interpretation lives entirely in your callback.

For when to use each kind of physics, see
[Phase-space evolution](phase_space_evolution.md). This page is the reference
for how they are spelled.

## The two shapes

Callbacks come in two shapes, and the difference is deliberate.

**Velocities fill a buffer the library supplies.** The buffer arrives already
sized; write every entry, and do not resize it.

```cpp
void(double time, const std::vector<double>& coordinates,
     std::vector<double>& velocities, std::size_t dimensions)
```

**Moment sources return a vector.** Build it in the natural way and return it.

```cpp
std::vector<double>(double time, const QmomSourceContext& context)
```

The reason is how often each runs. A moment source is evaluated once per SSP-RK
stage — three times per phase-space step — regardless of how large the phase
space is. A velocity is evaluated for every quadrature node, and the terminal
node count is the product of the per-dimension node counts, so the same
per-call allocation that is negligible for a source dominates for a velocity.
Write moment sources the easy way; there is nothing to gain from avoiding it.

## Reference

Physics callbacks, which most applications write:

| Callback | Signature | Purpose |
| --- | --- | --- |
| `PhaseSpaceVelocity` | `void(double, const std::vector<double>&, std::vector<double>&, std::size_t)` | Node velocity for all nodes at once |
| `PerNodePhaseSpaceVelocity` | `std::vector<double>(double, const std::vector<double>&)` | Convenience form, one node at a time |
| `QmomMomentSource` | `std::vector<double>(double, const QmomSourceContext&)` | QMOM raw-moment rates |
| `QmomMomentBoundaryContribution` | `QmomMomentRate(double, const std::vector<double>&)` | QMOM phase-space boundary exchange |
| `CqmomMomentSource` | `CqmomMomentRate(double, const CqmomSourceContext&)` | CQMOM raw-moment rates |
| `CqmomMomentBoundaryContribution` | `CqmomMomentRate(double, const CqmomSourceContext&)` | CQMOM boundary exchange |

Batch variants take a leading `std::size_t batch_index` and are otherwise
identical: `BatchPhaseSpaceVelocity`, `BatchQmomMomentSource`,
`BatchQmomMomentBoundaryContribution`, `BatchCqmomMomentSource`, and
`BatchCqmomMomentBoundaryContribution`. See
[Batch convenience wrappers](phase_space_evolution.md#batch-convenience-wrappers).

Lower-level callbacks, needed only when driving the time integrators or
realizability policies directly:

| Callback | Signature | Purpose |
| --- | --- | --- |
| `OdeRhs` | `void(double, const std::vector<double>&, std::vector<double>&)` | Generic ODE right-hand side |
| `StateAcceptancePredicate` | `bool(const std::vector<double>&)` | Accept or reject a completed SSP-RK3 step |
| `StageAcceptancePredicate` | `bool(SspRk3Stage, const std::vector<double>&)` | Accept or reject an individual stage |
| `MomentStateConverter` | `MomentSet(const std::vector<double>&)` | View a state vector as a moment set |

`QmomMomentRate` and `CqmomMomentRate` are thin structs holding a single
`values` vector.

## Phase-space velocity

The callback receives every terminal node coordinate at once, packed
node-major: node `i` occupies `[i * dimensions, (i + 1) * dimensions)`, and
within a node the values are ordered by phase-space dimension. Fill
`velocities` in exactly the same layout.

In one dimension the layout is simply one value per node:

```cpp
const qmomkit::PhaseSpaceVelocity velocity =
    [](double, const std::vector<double>& coordinates,
       std::vector<double>& velocities, std::size_t)
{
    for (std::size_t i = 0; i < coordinates.size(); ++i)
    {
        velocities[i] = -0.5 * coordinates[i];
    }
};
```

In more dimensions, step through the nodes and index within each one:

```cpp
const qmomkit::PhaseSpaceVelocity velocity =
    [](double, const std::vector<double>& coordinates,
       std::vector<double>& velocities, std::size_t dimensions)
{
    const std::size_t nodes = coordinates.size() / dimensions;
    for (std::size_t node = 0; node < nodes; ++node)
    {
        const std::size_t base = node * dimensions;
        const double diameter = coordinates[base + 0];
        const double speed = coordinates[base + 1];

        velocities[base + 0] = -1.0e-3 * diameter;
        velocities[base + 1] = -3.0 * speed / diameter;
    }
};
```

### Adapting a per-node callback

When the physics is naturally written for a single node,
`make_phase_space_velocity` adapts a `PerNodePhaseSpaceVelocity` into the
packed form:

```cpp
const qmomkit::PhaseSpaceVelocity velocity =
    qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& coordinate)
        {
            return std::vector<double>{-1.0e-3 * coordinate[0],
                                       -3.0 * coordinate[1] / coordinate[0]};
        });
```

This costs one call and one returned vector per node per evaluation, so prefer
the packed form for anything but small phase spaces or a quick experiment. See
[Adapting a per-node callback](phase_space_evolution.md#adapting-a-per-node-callback).

## Moment sources

A `QmomMomentSource` returns one rate per entry the context asks for, in context
order. `context.size()` is that count and `context.value(i)` is the current
value of entry `i` at this SSP-RK stage:

```cpp
const qmomkit::QmomMomentSource source =
    [](double, const qmomkit::QmomSourceContext& context)
{
    std::vector<double> rates(context.size(), 0.0);
    for (std::size_t order = 0; order < context.size(); ++order)
    {
        rates[order] = -0.1 * context.value(order);
    }
    return rates;
};
```

The CQMOM form is the same, wrapped in a `CqmomMomentRate`. Because CQMOM
moments are multi-indexed, use `context.index(entry)` to learn which moment an
entry is:

```cpp
const qmomkit::CqmomMomentSource source =
    [](double, const qmomkit::CqmomSourceContext& context)
{
    qmomkit::CqmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
    for (std::size_t entry = 0; entry < context.size(); ++entry)
    {
        const qmomkit::MultiIndexView order = context.index(entry);
        rate.values[entry] = -0.1 * static_cast<double>(order[0]) * context.value(entry);
    }
    return rate;
};
```

`context.index(entry)` returns a non-owning view. It indexes and iterates like
a `MultiIndex`, but it is only valid inside the callback; copy it with
`to_multi_index()` if you need to keep it. The full context interface is
documented under
[CqmomSourceContext](phase_space_evolution.md#cqmomsourcecontext).

## Boundary contributions

Boundary contributions describe exchange through a phase-space boundary. They
are added to the moment-source rates rather than forming a separate operator:

```cpp
const qmomkit::QmomMomentBoundaryContribution boundary =
    [](double, const std::vector<double>& moments)
{
    qmomkit::QmomMomentRate rate{std::vector<double>(moments.size(), 0.0)};
    rate.values[0] = 1.0e-3;
    return rate;
};
```

## Rules that apply to every callback

**Sizes are fixed by the library.** A velocity buffer arrives correctly sized
and must not be resized. A source must return exactly `context.size()` values.
A wrong size is reported as a failure, not silently padded.

**Every value must be finite.** A non-finite entry is reported with the
evaluation time and, for velocities, the offending node and dimension.

**Fill every entry of a supplied buffer.** A partially filled velocity buffer
is a contract violation. Debug builds detect it and report which node and
dimension was not written; release builds cannot detect it.

**Failures become a status, not an exception.** At the phase-space evolution
entry points, a callback that throws, returns a wrong-sized vector, or produces
a non-finite value yields a result with `ok() == false` and an explanatory
`message`, leaving the input state untouched. Check the result rather than
wrapping calls in `try`:

```cpp
const qmomkit::QmomPhaseSpaceEvolutionResult result =
    qmomkit::advance_qmom_phase_space(state, time, time_step, velocity, source, options);

if (!result.ok())
{
    std::cerr << result.message << '\n';
    return;
}
state = result.state;
```

**Callbacks are evaluated many times per step.** They run at every Runge-Kutta
stage, and rejected adaptive steps are retried, so a callback must not assume
it is called once per timestep, in increasing time order, or only at accepted
states. Keep them free of side effects that accumulate; write diagnostics from
the accepted result instead.

**Anything a callback captures must outlive the call.** The library stores the
`std::function` only for the duration of the operation, but a dangling
reference captured by `[&]` is still your responsibility.

## Where to go next

- [Phase-space evolution](phase_space_evolution.md): when to use each kind of
  physics, the splitting schemes, and the full step.
- [Configuration](configuration.md): settings that surround the callbacks.
- [Examples](examples.md): complete programs that define real callbacks.
