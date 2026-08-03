# Continuous Node-Advection Analytical Verification

These verification cases exercise the public continuous phase-space operator:
configuration loading, QMOM/CQMOM inversion, `qmomkit::advance_nodes(...)`,
and raw-moment reconstruction. Initial moment states are loaded from `.qmom`
files. Configuration files provide reproducible case setup and may define
supported numerical scheme choices such as operator-splitting scheme and
the node-advection scheme. The ODE coefficients, velocity callbacks, timestep
plan, and detailed RK23 tolerance settings live in the verifier because those
are application runtime choices rather than case-definition fields.

The continuous operator advances abscissas only. Quadrature weights represent
the amount of mass carried by each node, so they are invariant under these
deterministic phase-space trajectories.

`node_advection_analytical.cpp` is the canonical verification driver for this
area. It uses named `.qmom` configurations, closed-form node trajectories, and
reported error metrics.

Direct overload and diagnostic regression checks live under
`tests/component_regression/node_advection.cpp`.

## Constant Velocity

Configuration: `constant_velocity.qmom`

Initial 1D QMOM quadrature:

```text
nodes   = {-0.8, 0.4, 1.7}
weights = {0.2, 0.35, 0.45}
```

The ODE is:

```text
dxi/dt = 0.6
xi(t0 + tau) = xi(t0) + 0.6 tau
m_k(t0 + tau) = sum_a w_a xi_a(t0 + tau)^k
```

The verifier runs one step from `t=0.2` with `dt=0.15`, and seven outer
steps from `t=0.2` with `dt=0.12`.

## Linear Growth

Configuration: `linear_growth.qmom`

Initial 1D QMOM quadrature:

```text
nodes   = {-1.2, 0.3, 1.4}
weights = {0.15, 0.55, 0.30}
```

The ODE is:

```text
dxi/dt = 0.55 xi
xi(t0 + tau) = xi(t0) exp(0.55 tau)
m_k(t0 + tau) = exp(0.55 k tau) m_k(t0)
```

The verifier runs one step from `t=-0.1` with `dt=0.2`, and eight outer
steps from `t=-0.1` with `dt=0.1`. It also repeats a single `dt=0.9` run
with loose and tight RK23 tolerances and requires the tighter setting to reduce
the final node error.

## Affine Evolution

Configuration: `affine_evolution.qmom`

Initial 1D QMOM quadrature:

```text
nodes   = {0.2, 1.6, 2.8}
weights = {0.25, 0.45, 0.30}
```

The ODE is:

```text
dxi/dt = -0.4 xi + 0.7
xi_eq = 1.75
xi(t0 + tau) = (xi(t0) + 0.7 / -0.4) exp(-0.4 tau) - 0.7 / -0.4
```

The equilibrium lies inside the initial node range. The verifier runs one step
from `t=0.4` with `dt=0.2`, and six outer steps from `t=0.4` with `dt=0.15`.
Moments are reconstructed from the exact evolved nodes.

## Two-Dimensional Linear System

Configuration: `linear_system_2d.qmom`

The initial state is a 2D CQMOM moment tensor with
`conditioning_order = {1,0}`, `qmom.adaptive = true`, and
`adaptive_branch_counts = true`. The exact initial terminal quadrature is
ragged:

```text
(x, y) = (-0.7, -0.4), W = 0.45
(x, y) = ( 0.2,  1.1), W = 0.33
(x, y) = ( 1.4,  1.1), W = 0.22
```

The ODE is:

```text
dx/dt = -0.3 x
dy/dt = 0.8 x + 0.2 y
```

The exact solution is:

```text
x(t) = x0 exp(-0.3 t)
y(t) = y0 exp(0.2 t)
     + 0.8 x0 (exp(-0.3 t) - exp(0.2 t)) / (-0.3 - 0.2)
```

The verifier runs one step from `t=0.3` with `dt=0.1`, and eight outer steps
from `t=0.3` with `dt=0.08`. It checks terminal coordinates in phase-space
dimension order, unchanged terminal weights, the adaptive/ragged structure,
and all configured raw moments including cross moments.

## Numerical Settings

The standard RK23 tolerance and step-size settings are supplied by the
verifier. `.qmom` files may select the node-advection scheme, but these
verification cases keep the detailed adaptive embedded-RK controls in code together
with the velocity callback and timestep plan. The settings are:

```text
relative_tolerance = 1e-10
absolute_tolerance = 1e-12
max_step = 0.05
```

The linear-growth convergence sanity check uses a loose run with
`relative_tolerance = 1e-3`, `absolute_tolerance = 1e-6`, and a tight run with
the standard tolerances.

The pass tolerance is `5e-8` for node positions, weight changes, and
unit-floored relative moment errors. Absolute moment errors are also reported;
near-zero moments are governed by absolute error through the unit denominator
floor in the relative metric.

## Running

Configure and build the standalone test project first:

```bash
cd tests
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../install"
cmake --build build
```

Run all node-advection analytical cases:

```bash
ctest --test-dir build -R verification_node_advection_ --output-on-failure
```

Run individual cases:

```bash
ctest --test-dir build -R verification_node_advection_constant_velocity --output-on-failure
ctest --test-dir build -R verification_node_advection_linear_growth --output-on-failure
ctest --test-dir build -R verification_node_advection_affine --output-on-failure
ctest --test-dir build -R verification_node_advection_linear_system_2d --output-on-failure
```
