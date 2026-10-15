# Moment-source verification

## Purpose

The moment source advances raw moments through a caller-supplied rate with
fixed-step SSP-RK3
([Phase-space evolution](../../../docs/developer/phase_space_evolution.md#moment-source)).
These cases verify that `advance_qmom_moments_with_source` reproduces the
closed-form solutions of linear moment ODEs, for rates acting on each order
alone and for rates coupling orders, that it converges at third order, and
that the advanced moments remain invertible. The `common_rate_relaxation` case
advances the CQMOM moment vector through the vector overload of
`advance_qmom_moments_with_source`. Separate cases
call `advance_cqmom_moments_with_source` itself: with a boundary contribution
and with realizability policies for single coordinates.

## Cases

`verification_moment_source_<case>` (`moment_source_analytical.cpp`). The
initial moments, loaded from the fixture, are those of the weighted points
below; rates and step plans are set in the code. Each case runs one step and
a multi-step plan.

- `order_dependent_relaxation` (`order_dependent_relaxation_qmom.qmom`):
  $m_0$ to $m_5$ of $\xi = \{-0.6, 0.4, 1.3\}$, $w = \{0.25, 0.35, 0.40\}$
  relax toward the moments $m_k^*$ of $\xi = \{-0.4, 0.7, 1.5\}$ with the
  same weights, $\dot m_k = -0.15\,k\,(m_k - m_k^*)$, so
  $m_k(t_0+\tau) = m_k^* + (m_k(t_0) - m_k^*)\,e^{-0.15k\tau}$. One step of
  0.08 from $t = 0.2$; six steps of 0.04.
- `common_rate_relaxation` (`common_rate_relaxation_cqmom.qmom`, conditioning
  order $\{1, 0\}$, fixed node counts): the moments of the weighted points
  $(\xi_0, \xi_1) = (-0.5, -1.0), (-0.5, 0.8), (1.2, -1.0), (1.2, 0.8)$ with
  weights $0.20, 0.30, 0.15, 0.35$ relax toward those with weights
  $0.25, 0.20, 0.20, 0.35$, $\dot M_\alpha = -0.45\,(M_\alpha - M_\alpha^*)$
  for $\alpha \ne (0, 0)$. One step of 0.06 from $t = 0.1$; five steps of
  0.06.
- `linear_coupled` (`linear_coupled_source_qmom.qmom`): $m_0$ to $m_3$ of
  $\xi = \{0.5, 1.5\}$, $w = \{0.45, 0.55\}$ under a lower-triangular rate whose
  coupling of $m_2$ to $m_1$ guards against component-wise rate handling:

  ```math
  \dot m_0 = 0,\quad \dot m_1 = -0.2\,m_1,\quad \dot m_2 = 0.15\,m_1 + 0.1\,m_2,\quad \dot m_3 = 0.6\,m_3,
  ```

  with
  $m_2(t_0+\tau) = m_2(t_0)\,e^{0.1\tau} + 0.15\,m_1(t_0)\,(e^{-0.2\tau} - e^{0.1\tau})/(-0.3)$
  and exponentials for the others. One step of 0.08 from $t = -0.2$; eight
  steps of 0.06. A refinement study then advances $\tau = 0.48$ from $t = 0$
  with steps 0.12, 0.06, 0.03 and 0.015: the error must fall at every
  refinement, and the two finest observed orders must be at least 2.8.

After each plan the test checks:

- every advanced moment against the exact solution;
- that $m_0$, the number density, is unchanged, since its rate is zero;
- that every step is accepted at the requested size, without rejection or
  retry;
- that the advanced and the exact moments both invert, and that the inverted
  quadrature reproduces them; the CQMOM case also keeps its conditioning order
  and fixed storage.

`verification_moment_source_cqmom_update_<case>` (`cqmom_source_update.cpp`)
uses node counts $\{2, 2\}$ and references that SSP-RK3 reproduces exactly:

- `relaxation_matches_the_ssp_rk3_step`: the weighted points of
  `common_rate_relaxation` relax toward their target,
  $\dot M_\alpha = -r\,(M_\alpha - M^*_\alpha)$ with $r = 0.45$. One SSP-RK3
  step of $\Delta t$ is exactly $M^* + P(-r\Delta t)\,(M - M^*)$ with
  $P(z) = 1 + z + z^2/2 + z^3/6$; one step and five steps of 0.06 from
  $t = 0.1$. The rate reads each entry's multi-index from the source context.
- `source_and_boundary_are_integrated_exactly`: a constant source
  $0.5\,c^\alpha$ plus a boundary birth $(0.8 - 0.6t + 1.5t^2)\,b^\alpha$ at
  $c = (0.7, -0.2)$, $b = (0.3, -0.4)$. SSP-RK3 weighs its stage times
  $t$, $t+\Delta t$, $t+\Delta t/2$ like Simpson's rule, so one step of 0.3 from
  $t = 0.2$ adds the exact integral.
- Coordinate policies, on independent coordinates
  $M_{(k,j)} = A_k B_j$ with $B = (1, 0, 1, 0)$ and the rate $-r M_\alpha$ on
  every entry with second exponent 2, which multiplies $B_2$ by $P(-r\Delta t)$
  per step and by $1 - r\Delta t$ in the first stage; a real-line policy on the
  second coordinate rejects a negative $B_2$, over a step of 1.5:
  - `final_only_coordinate_policy_halves_a_rejected_step`: FinalOnly,
    $r = 2$, one retry: $B_2 = P(-1.5)^2$;
  - `every_stage_coordinate_policy_halves_a_rejected_step`: EveryStage,
    $r = 1$, one retry: $B_2 = P(-0.75)^2$;
  - `rejected_step_without_retries_returns_the_input`: EveryStage, $r = 1$, no
    retries: `rejected_by_predicate` and the input moments;
  - `failed_update_discards_an_accepted_part_of_the_interval`: EveryStage with
    $r = 1$ before $t = 0.75$ and $r = 2$ from then on, one retry, minimum
    time step 0.5: the half step to 0.75 is accepted, the rest fails with
    `minimum_timestep_reached`, and the update returns its input, not the
    moments at 0.75;
  - `policy_checks_only_its_own_coordinate`: the FinalOnly policy on the first
    coordinate accepts the full step, $B_2 = P(-3) = -2$.

## Tolerances

- Advanced moments: `2e-6` absolute; the unit-floored relative error checked
  alongside never exceeds it. The bound covers the SSP-RK3 truncation error at
  these steps, of order `1e-6` for the QMOM cases and `1e-8` for the CQMOM case.
  The two finest observed orders, the ones asserted, are within 0.03 of 3.
- Reconstructed moments and $m_0$: `2e-8`, a safety bound; the observed
  errors are at roundoff or zero.
- `cqmom_update` cases: `1e-13` relative to $\max(|\text{reference}|, 1)$,
  since the references are exact up to rounding; the observed errors are at
  most `2.2e-16`. Statuses and retry counts are compared exactly, and a failed
  update must return its input moments exactly.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_moment_source_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing.md#building-and-running).
