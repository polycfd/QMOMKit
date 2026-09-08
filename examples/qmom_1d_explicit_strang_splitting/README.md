# Explicit 1D QMOM Strang Splitting

This example demonstrates application-managed Strang splitting for adaptive
one-dimensional QMOM phase-space evolution. The application owns the stage sequence:

```text
C(dt/2) -> D(dt) -> C(dt/2)
```

and calls the public stage API directly.

The example-local velocity and source callbacks are defined in [physics.cpp](./physics.cpp);
the application-owned stage sequence remains visible in [main.cpp](./main.cpp).

## API Levels

`advance_qmom_phase_space(...)` is the library-managed complete phase-space evolution
workflow. It owns the configured Strang or Lie stage ordering for one macro
step.

`advance_qmom_advection(...)` and `advance_qmom_moment_source(...)` let the
surrounding application own the splitting sequence while still passing accepted
`QmomState` objects between stages.

## Numerical Distinction

The advection stage directly advects the current quadrature nodes and
reconstructs the stored raw moments from those same active nodes. It preserves
the active node count.

The moment-source-evolution stage integrates raw moments through the source callback and
then performs fixed or adaptive QMOM reconstruction. Adaptive reduction and
controlled one-node growth occur only after this moment-source reconstruction.

## Model

The phase-space advection model is linear relaxation:

$$
\frac{\mathrm{d}u}{\mathrm{d}t} = -a u.
$$

with $a = 500$. For zero mean, the variance decays as

$$
\operatorname{Var}[u](t) = \operatorname{Var}[u](0)\exp(-2at).
$$

The moment-source is a manufactured symmetric variance-injection source:

$$
\begin{aligned}
\frac{\mathrm{d}m_0}{\mathrm{d}t} &= 0, &
\frac{\mathrm{d}m_1}{\mathrm{d}t} &= 0, \\
\frac{\mathrm{d}m_2}{\mathrm{d}t} &= S(t), &
\frac{\mathrm{d}m_3}{\mathrm{d}t} &= 0.
\end{aligned}
$$

$S(t)$ is positive only over the finite interval $[0.16, 0.25)$. The source
preserves total mass and zero mean while allowing an adaptively reduced
one-node state to regrow to two nodes.

## Explicit Strang Sequence

The time loop keeps the three public stage calls visible:

```cpp
const qmomkit::QmomPhaseSpaceEvolutionResult first_advection =
    qmomkit::advance_qmom_advection(
        trial_state, time, 0.5 * time_step, velocity, phase_space_evolution_options);

const qmomkit::QmomPhaseSpaceEvolutionResult moment_source =
    qmomkit::advance_qmom_moment_source(
        first_advection.state, maximum_node_count,
        time + 0.5 * time_step, time_step, source, phase_space_evolution_options);

const qmomkit::QmomPhaseSpaceEvolutionResult second_advection =
    qmomkit::advance_qmom_advection(
        moment_source.state, time + 0.5 * time_step,
        0.5 * time_step, velocity, phase_space_evolution_options);
```

For a macro step beginning at `time`, the first advection half-step starts at
`time`, the moment-source stage starts at `time + 0.5 * time_step`, and
the second advection half-step also starts at `time + 0.5 * time_step`.

## Transactional Commit

The host application advances a trial state through all three stages. The
accepted macro-step state is updated only after the second advection
half-step succeeds. If a standalone stage fails, it rolls back only to its own
stage input, so the application keeps the previous macro-step state unchanged.

## Lifecycle

The shipped constants produce a deterministic adaptive lifecycle:

```text
2 -> 1 -> 2 -> 1
```

The exact transition times depend on selected tolerances and source
parameters. With this configuration, advection relaxation collapses the
variance while retaining the current active topology, a source-free
moment-source reconstruction reduces two nodes to one, the finite source pulse
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
