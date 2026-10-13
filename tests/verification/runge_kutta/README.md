# Runge–Kutta verification

## Purpose

These cases verify the generic ODE integrators of
[`runge_kutta.hpp`](../../../include/qmomkit/runge_kutta.hpp): SSP-RK3 and the
embedded Bogacki–Shampine 3(2) (`rk23`) and Dormand–Prince 5(4) (`rk45`)
pairs that node advection and the moment source build on
([Phase-space evolution](../../../docs/developer/phase_space_evolution.md#node-advection)).
Each method must reach its formal order under step refinement, and the
controller must deliver the requested accuracy. Every reference is a
closed-form ODE solution.

## Cases

`verification_runge_kutta_analytical_<case>` (`runge_kutta_analytical.cpp`)
compares the final state with the exact solution:

- `ssp_rk3_constant_rhs`, `rk23_constant_rhs`: a constant right-hand side,
  integrated exactly by one SSP-RK3 step, whose error estimate must be zero,
  and by RK23;
- `rk23_exponential_decay`: $y' = -y$ to $t = 2$, also checking the step and
  evaluation counts;
- `rk23_vector_coupled_linear_system`, `rk45_vector_coupled_linear_system`:
  $y_0' = y_1$, $y_1' = -y_0$ to $t = 1$, against $(\cos t, -\sin t)$.

`verification_<case>_convergence` (`runge_kutta_convergence.cpp`) measures
the absolute final error of $y' = -y$, $y(0) = 1$ against $e^{-T}$. Slopes come
from the finest pair of runs, since order is asymptotic.

- `ssp_rk3`: 5 to 80 fixed steps over $T = 1$; the observed order must exceed
  2.85.
- `rk23_fixed_step`, `rk45_fixed_step`: the controller is pinned to one step
  size (initial, minimum and maximum step equal, tolerances `1e-15`), 8 to 64
  steps over $T = 1$. The orders of the final error and of the single-step
  embedded error estimate must lie within 0.15 of 3 or 5; the estimate is set
  by the lower-order solution, whose local error is $O(h^p)$.
- `rk23_adaptive`, `rk45_adaptive`: relative tolerances `1e-4` to `1e-7`
  (absolute `1e-3` times that) over $T = 2$. The controller chooses the steps,
  so no order is asserted: the error must fall at every tightening, not exceed
  the relative tolerance from the second run on, and scale with it at a slope
  in $(0.85, 1.25)$; the accepted step count must not decrease.

## Tolerances

The constant right-hand sides are compared at `1e-15` and `1e-13`, roundoff.
The other closed-form cases use `5e-7`, `2e-7` and `2e-10` for relative
tolerances `1e-7`, `1e-8` and `1e-10`, since a controller bounds local error
estimates, not the global error; observed errors are `4e-8`, `8e-9` and
`1e-11`, so these are safety bounds.
Observed orders lie within 0.04 of nominal and slopes near 1.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_(runge_kutta_analytical_|(ssp_rk3|rk23|rk45)_.*convergence$)' --output-on-failure
```

See [Building and running](../../../docs/developer/testing.md#building-and-running).
