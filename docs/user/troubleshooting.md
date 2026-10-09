# Troubleshooting

QMOMKit reports problems rather than throwing or aborting. Almost every
operation returns a result carrying three things: `ok()`, a `status` naming the
category of failure, and a `message` describing the specific one. This page maps
the failures you are most likely to meet onto what to do about them.

## Reading a failure

Two checks are separate and both matter:

```cpp
const auto configured = qmomkit::config::load_qmom("configuration.qmom");
if (!configured)
{
    std::cerr << configured.error() << '\n';   // the file was not understood
    return 1;
}

const qmomkit::QmomResult result = qmomkit::compute_qmom(configured.value());
if (!result.ok())
{
    std::cerr << result.message << '\n';       // the numerics could not proceed
    return 1;
}
```

A file can be perfectly valid and still describe moments from which no
quadrature can be built. Configuration errors name a file and line; numerical
failures name the quantity that went wrong.

A failed evolution step does not corrupt anything. The returned result carries
the last accepted state, so you can inspect it, adjust, and retry from there
rather than restarting the run.

The status tells you *which stage* failed, which is the fastest way to narrow
down where to look:

| Status | Stage that failed | Usual cause |
| --- | --- | --- |
| `invalid_input` | Nothing ran | Sizes, counts, or non-finite values in what you supplied |
| `not_realizable` | QMOM inversion | The moments cannot come from a real distribution |
| `qmom_failed` | A CQMOM marginal or conditional inversion | As above, inside one branch |
| `conditional_moment_failed` | CQMOM conditioning | The conditional system could not be solved |
| `inversion_failed` | Reconstruction during a step | Evolved moments drifted out of the realizable set |
| `advection_failed` | Node advection | The velocity callback, or the node integrator |
| `moment_source_failed` | The moment-source update | The source callback, or its integrator |
| `moment_reconstruction_failed` | Rebuilding moments from the moved nodes | The moved quadrature could not be summed back into moments |
| `velocity_failed` | The velocity callback itself | A non-finite, wrong-sized, or partly filled result |
| `source_failed` | The source callback itself | As above, for rates |
| `rejected_by_predicate` | Source integration | A realizability or acceptance check rejected every retry |
| `minimum_timestep_reached` | Source integration | A rejected step could not be halved: retries reached `min_time_step`, or the step was already shorter |
| `step_size_underflow` | Adaptive node integration | The step at the floor produced a non-finite state or error estimate |
| `max_steps_exceeded` | Adaptive node integration | The interval needs more steps than `max_steps` allows |

## Reconstruction problems

### `Moment sequence has fewer than 2*n entries.`

An `n`-node quadrature needs at least `2n` raw moments, with or without
abscissa bounds and in fixed or adaptive mode. Supply more moments, or request
fewer nodes. When moments are generated from a distribution, raise
`maximum_order` to at least `2 * node_count - 1`.

### `Moment sequence has negative recurrence b = ...`

The moments are not realizable: no non-negative distribution has them. This is
the most common genuine failure, and it usually means the moments are wrong
rather than that QMOMKit needs different settings.

Check the moments themselves first. A frequent cause is mixing normalised and
unnormalised moments — QMOMKit uses raw, unnormalised moments throughout, with
`m0` the total number density. Another is moments assembled from separate
sources that are not mutually consistent.

To find out *which* moment is the problem, run the dedicated check:

```cpp
const auto check = qmomkit::check_realizability(moments, qmomkit::MomentSupport::real_line);
std::cerr << check.message << '\n';
```

which reports the offending order directly, for example:

```text
Moment sequence is not realizable: H_1 has a negative eigenvalue beyond
tolerance. Support: real_line.
```

See [Realizability](realizability.md) for supports and tolerances.

### `Moment sequence supports at most k nodes: ...`

The moment sequence carries less information than the requested node count
needs. The clearest case is a population concentrated at a single point: its
moments are consistent, but they support only one node, and asking for three
fails. The message names the number of nodes the moments do support.

**Use adaptive mode,** the default. It reconstructs with as many nodes as the
moments actually support instead of failing, so this message only appears when
the case asks for a fixed count with `[qmom] adaptive = false`. Removing that
setting lets the same case that fails at three fixed nodes succeed with one, and
lets individual CQMOM branches reduce independently.

### `Wheeler recurrence produced non-finite values.`

The recurrence overflowed: the moments, or the products the recurrence forms
from them, exceed the range of a double. This happens only for moments within
about a decade of the largest double, about $10^{308}$. Express the coordinate
in units that keep its moments smaller.

### `phase_space.adaptive_branch_counts was removed; use [qmom] adaptive.`

`[qmom] adaptive` decides whether a CQMOM reconstruction is adaptive, branch
node counts included. Replace the key with that setting; a fixed CQMOM case
uses:

```ini
[qmom]
adaptive = false
```

### `Moment sequence has m0 < cutoff.`

The zeroth moment is zero or negligible, so there is no population to represent.
Usually the initial condition was not set, or a source term removed all of it.
If a vanishing population is physically expected, stop evolving that case rather
than reconstructing from it.

### `Bounded-support Wheeler inversion returned nodes outside the configured interval.`

A configured support was declared but the moments imply abscissas outside it.
Either the bounds are wrong, or the moments are. Note that QMOMKit never clips
nodes into the interval — a candidate outside it is rejected, not projected. In
adaptive mode the reconstruction retries at lower node counts before giving up.

## Evolution problems

### The velocity callback is reported as non-finite

```text
Phase-space velocity callback returned a non-finite value:
  time: 0
  dimension: 0
  value: nan
  ...
```

The diagnostic names the time, the phase-space dimension, the offending value,
and the coordinates of the node that produced it. That node's coordinates are
usually the clue: a velocity dividing by a coordinate will produce this the
moment an abscissa reaches zero, and a bounded support or adaptive mode is
often the real fix rather than guarding the division.

A partially filled velocity buffer produces a similar report in debug builds,
naming the node and dimension that was never written. Release builds cannot
detect it, so develop callbacks against a debug build. See
[Callbacks](callbacks.md).

### `... exceeded max_steps.`

The adaptive node integrator needed more steps than allowed for the interval.
Either the dynamics are stiffer than the tolerances assume, or `max_steps` is
set too low:

```ini
[time_integration.node_advection.adaptive]
max_steps = 200000
```

Before raising it, check whether the velocity is far larger than intended —
exhausting a large step budget over a short interval usually means the model,
not the budget.

### `... produced a non-finite state or error estimate at the minimum step.`

The node integrator reached its smallest permitted step and the step taken there
produced a non-finite state or error estimate. This is divergence rather than a
tolerance problem, so relaxing `relative_tolerance` and `absolute_tolerance`
will not help; check the velocity field first, then whether the quadrature
feeding it still holds sensible coordinates.

Note that the *accuracy* case does not fail here. When the error controller
would ask for a step below the floor but the step at the floor is finite, the
integrator accepts it and continues, so a successful integration does not by
itself certify that the requested local-error tolerance was met. Lower `min_step`
if you need the controller to keep refining, and compare against a tighter run
to confirm the result is converged.

### Evolution succeeds for a while, then fails with `inversion_failed`

This is the characteristic pattern of a case whose moments slowly leave the
realizable set. The step itself is fine; the reconstruction after it is not.

In order of usefulness:

1. Turn on adaptive mode, so branches that lose information reduce their node
   count instead of failing.
2. Turn on [runtime realizability](realizability.md) checks so the source update
   rejects and retries a step that would produce unusable moments, rather than
   committing it:

   ```ini
   [time_integration.realizability]
   check_mode = final_only
   support = nonnegative
   ```

3. Reduce the timestep. Operator splitting has a finite order, and a source
   strong enough to break realizability in one step is usually being integrated
   too coarsely.
4. Declare a coordinate support if the population is physically bounded, so
   reconstruction cannot place nodes outside it.

### `rejected_by_predicate` or `minimum_timestep_reached`

A realizability or acceptance predicate rejected the source step, and halving
the timestep, as far as `min_time_step` allowed, did not help. The predicate is
doing its job: the source is driving the moments somewhere unusable. Treat this
as a physics or timestep problem rather than raising `max_retries`.

### A source or boundary callback is rejected

Messages of the form `... rate values must match the source context size.` or
`... rate values must be finite.` mean the callback's contract was broken. Return
exactly `context.size()` values, in context order, all finite. See
[Callbacks](callbacks.md).

## Getting more information

**Check realizability directly** rather than inferring it from a reconstruction
failure — `check_realizability(...)` names the offending moment order.

**Turn up CQMOM diagnostics** when a multidimensional case fails and you need to
know which branch:

```ini
[cqmom]
diagnostic_level = detailed
```

`basic` reports branch structure; `detailed` adds residual and condition
estimates. Both are observability only and change no numerical decision.

**Compare against a reference.** When results are suspicious rather than
failing, the [Monte Carlo companion](monte_carlo.md) computes moment histories
from a sampled ensemble for comparison. Remember it carries sampling
uncertainty; it is a reference, not exact truth.

**Check your node counts against your moments.** Many intermittent failures are
a case asking for more resolution than its moments carry. Fewer nodes are
cheaper and more robust, and adaptive mode makes that trade automatically.
