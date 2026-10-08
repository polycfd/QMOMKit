# Node-Advection Analytical Verification

These verification cases exercise the public phase-space advection operator:
configuration loading, QMOM/CQMOM inversion, `qmomkit::advance_nodes(...)`,
and raw-moment reconstruction. Initial moment states and reconstruction
settings are loaded from `.qmom` files. The node-advection scheme (RK23), its
adaptive controls, the ODE coefficients, the velocity callbacks, and the
timestep plan are set in the verifier, so each case controls them explicitly.
The configuration format can carry the node-advection scheme and adaptive
controls in `[time_integration.node_advection]`, but these cases do not use
it.

The advection operator advances abscissas only. Quadrature weights represent
the amount of mass carried by each node, so they are invariant under these
deterministic phase-space trajectories.

`node_advection_analytical.cpp` is the canonical verification driver for this
area. It uses named `.qmom` configurations, closed-form node trajectories, and
error metrics reported on failure.

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
Moments reconstructed from the advanced quadrature are compared with moments
computed from the exact evolved nodes.

## Two-Dimensional Linear System

Configuration: `linear_system_2d.qmom`

The initial state is a 2D CQMOM moment tensor with
`conditioning_order = {1,0}` and `qmom.adaptive = true`. The exact initial
terminal quadrature is
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

The standard RK23 tolerance and step-size settings are set in the verifier,
together with the velocity callback and timestep plan. The settings are:

```text
relative_tolerance = 1e-10
absolute_tolerance = 1e-12
max_step = 0.05
```

The linear-growth convergence sanity check uses a loose run with
`relative_tolerance = 1e-3`, `absolute_tolerance = 1e-6`, and a tight run with
the standard tolerances.

The pass tolerance is `5e-8` for node positions, weight changes, and
unit-floored relative moment errors, and the absolute moment error must not
exceed `1e-5`. Near-zero moments are governed by absolute error through the unit
denominator floor in the relative metric. The node error also includes the gap
between the integrator's final time and `t + dt`.

## Running

Configure and build the test suite from the repository root, in Debug as
described in [Testing](../../../docs/developer/testing.md#debug-build-requirement):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build --parallel
```

Run all node-advection analytical cases:

```bash
ctest --test-dir build -R '^verification_node_advection_' --output-on-failure
```

Run individual cases:

```bash
ctest --test-dir build -R '^verification_node_advection_constant_velocity$' --output-on-failure
ctest --test-dir build -R '^verification_node_advection_linear_growth$' --output-on-failure
ctest --test-dir build -R '^verification_node_advection_affine$' --output-on-failure
ctest --test-dir build -R '^verification_node_advection_linear_system_2d$' --output-on-failure
```
