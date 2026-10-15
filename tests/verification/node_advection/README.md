# Node-advection verification

## Purpose

Node advection moves every quadrature node along the phase-space velocity,
$\mathrm{d}\boldsymbol{\xi}_j/\mathrm{d}t = \dot{\boldsymbol{\xi}}(t, \boldsymbol{\xi}_j)$,
and keeps the weights, which carry the number density of each node, and the
quadrature structure
([Phase-space evolution](../../../docs/developer/phase_space_evolution.md#node-advection)).
These cases verify that `advance_nodes` reproduces closed-form node
trajectories in QMOM and CQMOM, leaves the weights and node count unchanged,
and that the moments of the advanced quadrature are those of the exact
trajectories.

## Cases

`verification_node_advection_<case>` (`node_advection_analytical.cpp`), with
fixture `<case>.qmom` (`affine` uses `affine_evolution.qmom`). Each fixture
holds the moments of the weighted points below, and the test first checks that
the fixture moments and the initial inversion match them. Velocities, step
plans and RK23 settings are set in the code. Each case runs one step and a
multi-step plan.

| Case | Initial weighted points | Velocity | $t_0$ | Steps |
| --- | --- | --- | --- | --- |
| `constant_velocity` | $\xi = \{-0.8, 0.4, 1.7\}$, $w = \{0.2, 0.35, 0.45\}$ | $\dot\xi = 0.6$ | 0.2 | 1 × 0.15; 7 × 0.12 |
| `linear_growth` | $\xi = \{-1.2, 0.3, 1.4\}$, $w = \{0.15, 0.55, 0.30\}$ | $\dot\xi = 0.55\,\xi$ | −0.1 | 1 × 0.2; 8 × 0.1 |
| `affine` | $\xi = \{0.2, 1.6, 2.8\}$, $w = \{0.25, 0.45, 0.30\}$ | $\dot\xi = -0.4\,\xi + 0.7$ | 0.4 | 1 × 0.2; 6 × 0.15 |
| `linear_system_2d` | $(\xi_0, \xi_1) = (-0.7, -0.4), (0.2, 1.1), (1.4, 1.1)$, $w = \{0.45, 0.33, 0.22\}$ | $\dot\xi_0 = -0.3\,\xi_0$, $\dot\xi_1 = 0.8\,\xi_0 + 0.2\,\xi_1$ | 0.3 | 1 × 0.1; 8 × 0.08 |

After a time $\tau$, the exact positions are $\xi(0) + 0.6\tau$,
$\xi(0)\,e^{0.55\tau}$ and $(\xi(0) - 1.75)\,e^{-0.4\tau} + 1.75$, with the
equilibrium inside the initial nodes, and, for the 2-D case,

```math
\xi_0(\tau) = \xi_0(0)\,e^{-0.3\tau},\qquad
\xi_1(\tau) = \xi_1(0)\,e^{0.2\tau} + 0.8\,\xi_0(0)\,\frac{e^{-0.3\tau} - e^{0.2\tau}}{-0.5}.
```

`linear_system_2d` uses conditioning order $\{1, 0\}$ and adaptive node
counts, so its initial inversion must give a ragged three-node quadrature.

After each plan the test checks node positions, in phase-space dimension order
for CQMOM, against the exact trajectories; that weights are unchanged and
positive and, in CQMOM, the terminal node count kept; that the integrator
reached $t + \Delta t$; and that all moments of the advanced quadrature,
including $m_0$ and the cross moments, match those of the exact nodes.
`linear_growth` also repeats a single step of 0.9 from $t = 0$ with loose
(`1e-3`, `1e-6`) tolerances and with those below, and requires the second run
to be more accurate.

## Tolerances

RK23 runs with relative tolerance `1e-10`, absolute `1e-12`, initial step
`0.01` and maximum step `0.05`. Node positions, the final-time gap, weight
changes and unit-floored relative moment errors must stay within `5e-8`, and
absolute moment errors within `1e-5`. The moments stay below 60, where the
relative bound is the tighter one, so the absolute bound is only a backstop.
Observed node errors are of order `1e-10` or below, so these are safety bounds
above the accuracy requested from the controller.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_node_advection_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing.md#building-and-running).
