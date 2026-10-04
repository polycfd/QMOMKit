# Runge-Kutta Analytical Verification

These checks verify the standalone ODE integration utilities used by
phase-space evolution-related components:

- SSP-RK3 on a constant right-hand side, where one fixed step is exact;
- adaptive RK23 on scalar exponential decay;
- adaptive RK23 on a constant vector right-hand side;
- adaptive RK23 and RK45 on a two-component harmonic oscillator over a short
  interval.

Expected values come from closed-form ODE solutions. The comparison tolerances
sit somewhat above the requested controller tolerances, after accounting for
the tested final time and solution scale, because an adaptive integrator
controls its local error estimate rather than the global error.

`runge_kutta_convergence.cpp` contains five independently registered studies
using the scalar closed-form problem

```text
y' = -y
y(0) = 1
y(T) = exp(-T)
```

- `verification_ssp_rk3_convergence` verifies fixed-step SSP-RK3 formal
  convergence;
- `verification_rk23_adaptive_convergence` studies tolerance-controlled
  adaptive RK23 behavior;
- `verification_rk23_fixed_step_convergence` verifies constrained fixed-step
  RK23 formal convergence;
- `verification_rk45_fixed_step_convergence` verifies constrained fixed-step
  RK45 formal convergence;
- `verification_rk45_adaptive_convergence` studies tolerance-controlled
  adaptive RK45 behavior.

For SSP-RK3, the verifier advances from `T=0` to `T=1` with repeated fixed
steps and halves the timestep over several refinement levels. The error norm is
the absolute final-state error `|y_num(T) - exp(-T)|`. The observed order from
the finest pair of refinements must exceed `2.85`.

For adaptive RK23 and RK45, the verifier integrates the same ODE from `T=0` to
`T=2` with tolerances tightened tenfold per run. An explicit initial step is
supplied for each run, while the solver chooses later accepted timesteps. These
studies deliberately assert no convergence order, because the controller, not
a refined step size, chooses the steps. They require the achieved error to stay
within the requested tolerance, to fall at every tightening, the accepted step
count not to decrease, and the error-versus-tolerance slope from the finest
pair of runs to lie between `0.85` and `1.25`.

The constrained fixed-step studies pin the controller to a single step size and
halve it exactly. From the finest pair of refinements, both the final-error
order and the order of the single-step embedded error estimate must lie within
`0.15` of the nominal order: three for RK23 and five for RK45.

Nothing is printed on success. On failure the reason and the full refinement
table go to standard error.

Run with:

```bash
ctest --test-dir build -R '^verification_runge_kutta_analytical_' --output-on-failure
ctest --test-dir build -R '^verification_(ssp_rk3|rk23_adaptive|rk23_fixed_step|rk45_adaptive|rk45_fixed_step)_convergence$' --output-on-failure
```

Build the test suite first as described in
[Testing](../../../docs/developer/testing.md#debug-build-requirement).
