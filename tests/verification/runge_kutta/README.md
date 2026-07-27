# Runge-Kutta Analytical Verification

These checks verify the standalone ODE integration utilities used by
transport-related components:

- SSP-RK3 on a constant right-hand side, where one fixed step is exact;
- adaptive RK23 on scalar exponential decay;
- adaptive RK23 on a constant vector right-hand side;
- adaptive RK23 on a two-component harmonic oscillator over a short interval.

Expected values come from closed-form ODE solutions. The tolerances are chosen
to be tighter than the configured adaptive tolerances after accounting for the
tested final time and solution scale.

`runge_kutta_convergence.cpp` is the formal order-of-accuracy verification
case. It uses the scalar closed-form problem

```text
y' = -y
y(0) = 1
y(T) = exp(-T)
```

For SSP-RK3, the verifier advances from `T=0` to `T=1` with repeated fixed
steps and halves the timestep over several refinement levels. The error norm is
the absolute final-state error `|y_num(T) - exp(-T)|`. The printed table reports
timestep size, error, and observed order, and the finest refinements are
required to approach third-order global convergence.

For adaptive RK23, the verifier integrates the same ODE from `T=0` to `T=2`
with successively tighter relative and absolute tolerances. The solver chooses
its own accepted timesteps. The printed table reports tolerance, accepted and
rejected steps, function evaluations, average accepted step, final error, and
observed order against the average accepted step. The check requires tighter
tolerances to reduce the error, increase or preserve accepted-step counts, and
show a high-order trend for the finest refinements.

Run with:

```bash
ctest --test-dir build -R "verification_runge_kutta_(analytical|convergence)" --output-on-failure
```
