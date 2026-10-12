# Evolving in time

This page covers the second half of a workflow: advancing a state through time
with your physics, checking each step, and writing results. It continues the
time loop of [Getting started](2_getting_started.md).

## The time step

`advance_qmom_phase_space(...)` and `advance_cqmom_phase_space(...)` advance a
state by one step and return the result:

```cpp
const qmomkit::QmomPhaseSpaceEvolutionResult result =
    qmomkit::advance_qmom_phase_space(state, time, time_step, velocity, source, options);
```

A step combines two operations by operator splitting:

- **node advection** moves the quadrature nodes with the velocity through
  phase space and keeps their weights; the moments then follow from the moved
  nodes;
- **the source step** integrates the moments with the source rates, using the
  strong-stability-preserving SSP-RK3 scheme
  [[Shu & Osher, 1988](https://www.doi.org/10.1016/0021-9991(88)90177-5)]; the
  **post-source reconstruction** then reconstructs the quadrature from the new
  moments.

The splitting scheme sets their order, and a separate scheme integrates the
node advection:

| Setting | Values |
| --- | --- |
| Splitting, `options.scheme` | `strang` (default): half a step of node advection, the source step, half a step of node advection [[Strang, 1968](https://www.doi.org/10.1137/0705041)]. `advsrc` and `srcadv`: node advection then source, or the reverse, first-order Lie splitting. |
| Node advection, `options.node_options.scheme` | `rk23` (default), adaptive [[Bogacki & Shampine, 1989](https://www.doi.org/10.1016/0893-9659(89)90079-7)]; `rk45`, adaptive [[Dormand & Prince, 1980](https://www.doi.org/10.1016/0771-050X(80)90013-3)]; `ssp_rk3`, three stages per step, no step control. |

Only node advection with `rk23` or `rk45` adapts its steps to a tolerance. The
source step is explicit and has no error control: its accuracy and stability
depend on the time step you choose, so a stiff source needs a small step.

The options can be loaded from a `.qmom` file and changed in code:

```cpp
const auto loaded = qmomkit::config::load_qmom_phase_space_evolution_options("qmom.qmom");
if (!loaded)
{
    std::cerr << loaded.error() << '\n';
    return 1;
}
qmomkit::QmomPhaseSpaceEvolutionOptions options = loaded.value();
options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
options.source_options.ssp_rk3_options.max_retries = 4;
```

The evolution options describe only the step. How quadratures are
reconstructed, adaptivity and bounds included, belongs to the case, and every
state carries it forward. The adaptive schemes take tolerances and step limits;
the [Reference](7_reference.md#evolution-options) lists them with their
defaults.

A source step that a realizability policy rejects is retried with half the
previous step size, at most `max_retries` times per attempt (none by default)
and not below `min_time_step`. After an accepted part, a new attempt covers the
rest of the interval, so a whole step can retry more often than `max_retries`.

## Callbacks

Your physics enters through three kinds of callback. Supply the ones your model
has; the [Reference](7_reference.md#callbacks) lists their signatures.

### Velocity

The phase-space velocity defines how fast the nodes move through phase space,
for example by growth or drag. One call receives every node, packed node by
node, and fills a buffer of the same layout:

```cpp
const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinates, std::vector<double>& velocities,
                                                std::size_t dimensions)
{
    for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
    {
        const double size = coordinates[node * dimensions + 0];
        const double speed = coordinates[node * dimensions + 1];
        velocities[node * dimensions + 0] = -0.1 * size;
        velocities[node * dimensions + 1] = -speed;
    }
};
```

Write every entry, and do not resize the buffer.
`make_phase_space_velocity(...)` adapts a function of one node instead, as in
[Getting started](2_getting_started.md); it is simpler, but costs one call per
node.

### Source

A source returns the rate of change of each moment, for example from breakage,
aggregation, or removal. Its context holds the moments of the current
Runge–Kutta stage:

| QMOM context | CQMOM context | Gives |
| --- | --- | --- |
| `size()` | `size()` | the number of rates to return |
| `value(k)` | `value(i)`, `index(i)` | the value of entry `k` or `i`, and for CQMOM its multi-index |
| `moment(k)` | `moment(index)` | any integer-order moment, tracked or not |
| `generalized_moment(order)` | `generalized_moment(orders)` | a moment of non-negative real order |
| `for_each_node(...)` | `for_each_terminal_node(...)` | each node's weight and coordinates |

Moments the state does not track, generalized moments, and nodes come from a
quadrature reconstructed from the stage, at most once per stage. Within a
source step, a branch of the state never gains a node.

A rate is often easiest to write over the nodes. Here particles are removed at a
rate proportional to their size, $\mathrm{d}m_k/\mathrm{d}t = -\sum_i 0.1\,\xi_i\,w_i\,\xi_i^k$:

```cpp
const qmomkit::QmomMomentSource removal = [](double, const qmomkit::QmomSourceContext& context)
{
    qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
    context.for_each_node(
        [&rate](double weight, double xi)
        {
            double power = 1.0; // xi^k
            for (double& value : rate.values)
            {
                value -= 0.1 * xi * weight * power;
                power *= xi;
            }
        });
    return rate;
};
```

The same removal for CQMOM, written over moments, is
$\mathrm{d}M_{k,l}/\mathrm{d}t = -0.1\,M_{k+1,l}$. Return one rate per context
entry, in context order:

```cpp
const qmomkit::CqmomMomentSource cqmom_removal = [](double, const qmomkit::CqmomSourceContext& context)
{
    qmomkit::CqmomMomentRate rate{std::vector<double>(context.size())};
    for (std::size_t i = 0; i < context.size(); ++i)
    {
        qmomkit::MultiIndex order = context.index(i).to_multi_index();
        order[0] += 1; // dM_{k,l}/dt = -0.1 M_{k+1,l}
        rate.values[i] = -0.1 * context.moment(order);
    }
    return rate;
};
```

A model without a source passes one that returns zero rates. A model without a
velocity passes one that writes zeros, or calls
`advance_qmom_moment_source(...)` alone, as in
[Coupling QMOMKit](6_coupling_qmomkit.md#driving-the-stages-yourself).

### Boundary contribution

Exchange across a phase-space boundary, such as nucleation or injection, is a
boundary contribution: a callback shaped like a source, whose rates are added
to it. Here particles are injected at size 0.1 at a rate of 0.02:

```cpp
const qmomkit::QmomMomentBoundaryContribution injection = [](double, const qmomkit::QmomSourceContext& context)
{
    qmomkit::QmomMomentRate rate{std::vector<double>(context.size())};
    double power = 1.0; // 0.1^k
    for (double& value : rate.values)
    {
        value = 0.02 * power;
        power *= 0.1;
    }
    return rate;
};
const qmomkit::QmomPhaseSpaceEvolutionResult result =
    qmomkit::advance_qmom_phase_space(state, time, time_step, velocity, removal, injection, options);
```

### Rules for every callback

- Return exactly the number of values requested, all finite.
- A callback runs many times per step, at every Runge–Kutta stage and for
  retried steps, so it must not count calls or assume increasing time.
- A callback that throws a `std::exception`, returns the wrong size, or
  returns a non-finite value makes the step fail with a message, without an
  exception reaching you. Other exceptions propagate.
- A source context, and the coordinates a node visitor receives, are valid
  only during the call; do not keep references to them.
- A Debug build of QMOMKit detects an unwritten velocity entry; a Release build
  cannot, so develop callbacks against a Debug build of QMOMKit.

## Adaptive node counts

In adaptive mode a state remembers how many nodes each branch currently has.
Node advection never changes these counts; only the reconstruction after a
source step does. There a branch may lose nodes, or gain one if the moments
support it. In CQMOM, when a branch's node count changes, the branches below it
condition on other nodes and start from their configured maximum. Node counts
therefore shrink readily and grow by at most one per step.

## Reading a state

A QMOM state gives its moments $m_0, m_1, \ldots$ with `moment_values()` and its
quadrature with `quadrature()`; `active_node_count()` and `maximum_node_count()`
give the current and the configured number of nodes. A CQMOM state lists its
moments in the order of `moment_indices()`, and `moment(index)` returns one of
them if the state holds it:

```cpp
const std::vector<qmomkit::MultiIndex> indices = state.moment_indices();
const std::vector<double>& values = state.moment_values();
for (std::size_t i = 0; i < values.size(); ++i)
{
    std::cout << "M[" << indices[i][0] << "," << indices[i][1] << "] = " << values[i] << '\n';
}
const std::optional<double> m10 = state.moment(qmomkit::MultiIndex{1, 0});
```

`quadrature()` gives its nodes, and `active_branch_counts()` the current node
count of each branch.

## Checking the result

```cpp
if (!result.ok())
{
    std::cerr << result.message << " (stage " << static_cast<int>(result.diagnostics.failed_stage) << ")\n";
    return 1;
}
state = result.state;
time = result.final_time;
```

On failure, `result.state` is the state the step started from and
`result.final_time` its time. `status` says what failed and
`diagnostics.failed_stage` in which part of the step; the
[Reference](7_reference.md#statuses) lists both. An invalid option, such as a
negative tolerance, is reported as `invalid_input` before anything runs, so a
step failure always concerns the numerics or your callbacks. Retrying with a
smaller step, or stopping, is your decision. A failed CQMOM step also returns,
for diagnosis only, a `diagnostic_quadrature`: the quadrature a failed
reconstruction produced, or the input of a failed advection or source stage.

## Writing output

`MomentWriter` writes the moments of a state, and `QuadratureWriter` its nodes,
as text tables with one row per time, or per node and time:

```cpp
const auto output = qmomkit::config::load_output_options("qmom.qmom");
qmomkit::MomentWriter moment_writer = qmomkit::MomentWriter::for_qmom("moments.txt", state, output.value());
qmomkit::QuadratureWriter node_writer = qmomkit::QuadratureWriter::for_qmom("nodes.txt", state.maximum_node_count());

// After each step:
moment_writer.write(time, state, state.quadrature());
node_writer.write(time, state.quadrature());
```

For CQMOM, use `MomentWriter::for_cqmom("moments.txt", state, output_options)`
and `QuadratureWriter::for_cqmom("nodes.txt", state.hierarchy().number_of_dimensions())`.
`[output] generalized_moments` adds columns of generalized moments to the
moment file. When to write is up to your loop; `[output] interval` and
`directory` only carry your choice from the file to it. A writer throws if a
file cannot be written. The [Reference](7_reference.md#output-files) shows the
file formats.
