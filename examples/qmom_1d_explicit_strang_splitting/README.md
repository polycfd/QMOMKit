# Explicit 1D QMOM Strang Splitting

This example demonstrates application-managed Strang splitting for adaptive
one-dimensional QMOM transport. The application owns the stage sequence:

```text
C(dt/2) -> D(dt) -> C(dt/2)
```

and calls the public stage API directly.

## API Levels

`advance_transport_qmom(...)` is the library-managed complete transport
workflow. It owns the configured Strang or Lie stage ordering for one macro
step.

`advance_qmom_continuous(...)` and `advance_qmom_discontinuous(...)` let the
surrounding application own the splitting sequence while still passing accepted
`QmomState` objects between stages.

## Numerical Distinction

The continuous stage directly advects the current quadrature nodes and
reconstructs the stored raw moments from those same active nodes. It preserves
the active node count.

The discontinuous stage integrates raw moments through the source callback and
then performs fixed or adaptive QMOM reconstruction. Adaptive reduction and
controlled one-node growth occur only after this discontinuous reconstruction.

## Model

The continuous phase-space model is linear relaxation:

```text
du/dt = -a u
```

with `a = 500`. For zero mean, the variance decays as
`Var[u](t) = Var[u](0) exp(-2*a*t)`.

The discontinuous source is a manufactured symmetric variance-injection source:

```text
dm0/dt = 0
dm1/dt = 0
dm2/dt = S(t)
dm3/dt = 0
```

`S(t)` is positive only over the finite interval `[0.16, 0.25)`. The source
preserves total mass and zero mean while allowing an adaptively reduced
one-node state to regrow to two nodes.

## Explicit Strang Sequence

The time loop keeps the three public stage calls visible:

```cpp
const qmomkit::QmomTransportResult first_continuous =
    qmomkit::advance_qmom_continuous(
        trial_state, time, 0.5 * time_step, velocity, transport_options);

const qmomkit::QmomTransportResult discontinuous =
    qmomkit::advance_qmom_discontinuous(
        first_continuous.state, maximum_node_count,
        time + 0.5 * time_step, time_step, source, transport_options);

const qmomkit::QmomTransportResult second_continuous =
    qmomkit::advance_qmom_continuous(
        discontinuous.state, time + 0.5 * time_step,
        0.5 * time_step, velocity, transport_options);
```

For a macro step beginning at `time`, the first continuous half-step starts at
`time`, the discontinuous source step starts at `time + 0.5 * time_step`, and
the second continuous half-step also starts at `time + 0.5 * time_step`.

## Transactional Commit

The host application advances a trial state through all three stages. The
accepted macro-step state is updated only after the second continuous
half-step succeeds. If a standalone stage fails, it rolls back only to its own
stage input, so the application keeps the previous macro-step state unchanged.

## Lifecycle

The shipped constants produce a deterministic adaptive lifecycle:

```text
2 -> 1 -> 2 -> 1
```

The exact transition times depend on selected tolerances and source
parameters. With this configuration, continuous relaxation collapses the
variance while retaining the current active topology, a source-free
discontinuous reconstruction reduces two nodes to one, the finite source pulse
regrows the state to two nodes, and later source-free reconstruction reduces it
again.

## Build And Run

From the repository root, install QMOMKit, then configure this standalone
example against the installed package:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cd examples/qmom_1d_explicit_strang_splitting
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_1d_explicit_strang_splitting_example
```

## Output

The executable writes:

- `output/quadrature_history.txt`

The file is written through `qmomkit::QuadratureWriter::for_qmom(...)` after
the initial accepted state and after each completed Strang macro-step.

`quadrature_history.txt` records one row per active node:

```text
time node_index weight x0
```

Inactive adaptive nodes are not written as physical quadrature nodes.
