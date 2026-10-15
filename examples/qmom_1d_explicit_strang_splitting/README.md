# 1D QMOM with explicit Strang splitting

This example drives the stages of a Strang step itself instead of calling
`advance_qmom_phase_space(...)`, as described in [Driving the stages
yourself](../../docs/user/6_coupling_qmomkit.md#driving-the-stages-yourself).
Each step of length $\Delta t$ runs

```text
advection(Δt/2) -> source(Δt) -> advection(Δt/2)
```

with `advance_qmom_advection(state, t, Δt/2, ...)`,
`advance_qmom_moment_source(..., t, Δt, ...)` and
`advance_qmom_advection(..., t + Δt/2, Δt/2, ...)`, the sequence of the
`strang` scheme, so the results equal those of `advance_qmom_phase_space(...)`.
Between the stages it reports when the adaptive reconstruction after the
source changes the number of active nodes.

## Model

The phase-space coordinate is a particle velocity $u$; $u$ and the time $t$ are
dimensionless. Node advection relaxes the velocities toward
$u_\mathrm{eq} = 0.5$,

```math
\frac{\mathrm{d}u}{\mathrm{d}t} = -a\left(u - u_\mathrm{eq}\right),
\qquad a = 500,
```

so the variance decays as $\exp(-2at)$. The source adds variance at constant
number density and mean,

```math
\frac{\mathrm{d}m_0}{\mathrm{d}t} = \frac{\mathrm{d}m_1}{\mathrm{d}t}
= \frac{\mathrm{d}m_3}{\mathrm{d}t} = 0,
\qquad
\frac{\mathrm{d}m_2}{\mathrm{d}t} = S(t),
```

with $S = 12$ for $0.16 \le t < 0.25$ and $S = 0$ otherwise. $a$,
$u_\mathrm{eq}$ and $S$ are set in [physics.cpp](./physics.cpp). The nonzero
$u_\mathrm{eq}$ matters: relaxation toward zero would only shrink the
distribution without changing its shape, and adaptive reduction does not
depend on the unit of the coordinate, so the state would keep both nodes.

[configuration.qmom](./configuration.qmom) gives the moments
$m_0, \ldots, m_3 = 1, 0, 1, 0$, a zero-mean unit-variance population, on two
nodes with `adaptive = true`, so the reconstruction after the source may drop
or regain a node. It also sets the `strang` scheme, which this example
reproduces by hand, and the tolerances of the `rk23` node advection.
[main.cpp](./main.cpp) advances from $t = 0$ to $0.5$ with steps of $0.1$ and
replaces its state only after all three stages of a step succeed.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_1d_explicit_strang_splitting_example
```

## Output

The example prints the initial number of active nodes, then, for each step in
which the reconstruction after the source changes it, the time $t + \Delta t/2$
of the source stage and the counts before and after, and finally the final
number of active nodes and the output directory. It writes the quadrature
after the initial state and after each step to
`build/output/quadrature_history.txt`; setting `write_moments = true` in the
`[output]` section of `configuration.qmom` adds
`build/output/moment_history.txt`.

## What to expect

The active node count goes $2 \to 1 \to 2 \to 1$, with the changes reported at
$t = 0.05$, $0.15$ and $0.35$. Relaxation collapses the two nodes onto
$u = 0.5$ during the first step, and the reconstruction after the source
reduces them to one. The source pulse adds variance and the state regains a
second node, which the reconstruction drops again once relaxation has removed
that variance. The run ends with one node of weight 1 at $u = 0.5$.
