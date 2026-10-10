# Coupling QMOMKit

QMOMKit advances one population at a time. In a simulation code, your code owns
the rest: cells and meshes, transport in physical space, parallelism, the time
loop, and what to do when a step fails. This page shows the QMOMKit side of
such a coupling. Examples of complete couplings to host codes are maintained in
[QMOMKit-Couplings](https://github.com/polycfd/QMOMKit-Couplings).

## Moments changed by your code

When your code changes the moments itself, for example by transporting them
from cell to cell, it needs a new state for the changed moments: a state cannot
be edited. Keep the moments of each cell between steps, not a state, and in
every step:

1. transport the moments;
2. create a state from them with `state_candidate(...)`;
3. advance the state;
4. store the moments of the advanced state.

Start each cell from your own moments, in the case's moment order, or from
`initial_state().moment_values()`.

```cpp
const qmomkit::QmomCase qmom_case = qmomkit::config::load_qmom("qmom.qmom").value();
const qmomkit::QmomPhaseSpaceEvolutionOptions options =
    qmomkit::config::load_qmom_phase_space_evolution_options("qmom.qmom").value();

for (std::vector<double>& moments : cell_moments)
{
    qmomkit::QmomStateCandidate candidate = qmom_case.state_candidate(qmomkit::MomentSet(moments), options.realizability_policy);
    if (candidate.status == qmomkit::StateCandidateStatus::empty)
    {
        continue; // nothing to evolve in an empty cell
    }
    if (!candidate.ok())
    {
        std::cerr << candidate.message << '\n'; // your decision: retry the transport step, or stop
        return 1;
    }
    const qmomkit::QmomPhaseSpaceEvolutionResult evolved =
        qmomkit::advance_qmom_phase_space(*candidate.state, time, time_step, velocity, source, options);
    if (evolved.ok())
    {
        moments = evolved.state.moment_values();
    }
}
```

The candidate uses the node count, options, and bounds of the case, with the
moments it is given. It reconstructs the quadrature once and does not throw:
`status` says whether a state was created, the cell was `empty`, or the moments
were rejected, and the [Reference](7_reference.md#statuses) lists all cases.
Passing a realizability policy is optional; it checks the moments before the
reconstruction.

For CQMOM, create the moment data once and replace only its values, in the
order of its layout's `required_indices`. Every state then shares the prepared
case's structure:

```cpp
const qmomkit::CqmomPreparedCase prepared_case = qmomkit::config::load_cqmom_prepared("configuration.qmom").value();
const qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_options =
    qmomkit::config::load_cqmom_phase_space_evolution_options("configuration.qmom").value();
qmomkit::CqmomMomentData moment_data{qmomkit::create_cqmom_moment_layout(prepared_case.hierarchy()), {}};

for (std::vector<double>& moments : cqmom_cells)
{
    moment_data.values = moments;
    qmomkit::CqmomStateCandidate candidate =
        prepared_case.state_candidate(moment_data, cqmom_options.source_options.coordinate_realizability_policies);
    if (candidate.ok())
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult evolved =
            qmomkit::advance_cqmom_phase_space(*candidate.state, time, time_step, velocity, cqmom_source, cqmom_options);
        if (evolved.ok())
        {
            moments = evolved.state.moment_values();
        }
    }
}
```

A new state stores the moments of its quadrature. They equal the given moments
up to order $2N - 1$, or for CQMOM the required moments, unless an adaptive
reconstruction dropped nodes the moments do not support; every other moment is
that of the quadrature. Storing `moment_values()` back keeps each cell
consistent with what was evolved.
The same pattern restarts a run: store the moments of each state, and create
the states again from them with `state_candidate(...)`.

Keeping the moments realizable during transport is your code's part. A scheme
whose new cell values are convex combinations of realizable neighbour values,
such as first-order upwinding with a Courant number of at most one, keeps them
realizable; higher-order schemes need a realizability-preserving
reconstruction. Without a policy, an adaptive reconstruction quietly reduces
non-realizable moments: `{1, 0, -1, 0}`, with a negative variance, becomes a
single node at 0. A policy rejects them instead. It checks moments, not node
positions, and for CQMOM only each coordinate's own moments; see
[Moments and quadrature](4_moments_and_quadrature.md#check-realizability).

The [`qmom_state_candidate`](../../examples/qmom_state_candidate/README.md)
example shows the outcomes.

## Many states at once

The batch functions advance a vector of states with one time, time step, and set of
options, and return one result per state, in order:

```cpp
std::vector<qmomkit::QmomState> states(3, qmom_case.initial_state());
const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> results =
    qmomkit::advance_qmom_phase_space_batch(states, time, time_step, velocity, source, options);
```

When the physics differs between states, the `_with_index` variants pass each
callback the position of the state in the vector:

```cpp
const qmomkit::BatchQmomMomentSource indexed_source =
    [&local_rate](std::size_t batch_index, double, const qmomkit::QmomSourceContext& context)
{
    qmomkit::QmomMomentRate rate{std::vector<double>(context.size())};
    for (std::size_t k = 0; k < context.size(); ++k)
    {
        rate.values[k] = -local_rate[batch_index] * context.value(k);
    }
    return rate;
};
const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> indexed =
    qmomkit::advance_qmom_phase_space_batch_with_index(states, time, time_step, indexed_velocity, indexed_source, options);
```

Batches run serially; a failed state does not stop the others. To run states in
parallel, call the single-state functions from your own parallel loop. QMOMKit
keeps no global mutable state, so calls on different states can run at the same
time; your callbacks must then be safe to call concurrently. The
[`cqmom_2d_batch_host_data`](../../examples/cqmom_2d_batch_host_data/README.md)
example uses host data per state.

## Driving the stages yourself

When your code must interleave its own work within a step, call the stages
directly. This is the Strang step:

```cpp
const qmomkit::QmomPhaseSpaceEvolutionResult first =
    qmomkit::advance_qmom_advection(state, time, 0.5 * time_step, velocity, options);
if (!first.ok())
{
    return 1;
}
const qmomkit::QmomPhaseSpaceEvolutionResult middle =
    qmomkit::advance_qmom_moment_source(first.state, time, time_step, source, options);
if (!middle.ok())
{
    return 1;
}
const qmomkit::QmomPhaseSpaceEvolutionResult last =
    qmomkit::advance_qmom_advection(middle.state, time + 0.5 * time_step, 0.5 * time_step, velocity, options);
```

`advance_cqmom_advection(...)` and `advance_cqmom_moment_source(...)` do the
same for CQMOM. The
[`qmom_1d_explicit_strang_splitting`](../../examples/qmom_1d_explicit_strang_splitting/README.md)
example runs this loop.

## Without a state

The operators beneath the stages work on a quadrature or a moment vector
directly, without a state:

```cpp
const qmomkit::QmomResult quadrature = qmomkit::compute_qmom(qmomkit::MomentSet({1.0, 0.8, 2.8, 4.4}), 2);
const qmomkit::QmomNodeAdvectionResult moved = qmomkit::advance_nodes(quadrature, time, time_step, velocity);
const std::vector<double> moved_moments = qmomkit::reconstruct_raw_moments(moved.quadrature, 4);

const qmomkit::QmomMomentSourceUpdateResult updated =
    qmomkit::advance_qmom_moments_with_source(moved_moments, time, time_step, source);
```

Without a state there are no reconstruction settings, so the source context
cannot reconstruct a quadrature: `moment(...)` works only for tracked orders,
and `generalized_moment(...)` and the node visitors are unavailable.

## Building QMOMKit with your code

QMOMKit can be built as part of your project instead of installed; see
[Installation](1_installation.md#build-qmomkit-inside-your-project) for how its
defaults then apply to your project.
