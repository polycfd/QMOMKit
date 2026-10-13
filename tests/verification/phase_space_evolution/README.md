# Phase-space evolution verification

## Purpose

These cases verify the
[phase-space evolution](../../../docs/developer/phase_space_evolution.md)
steps of QMOM and CQMOM states. Commuting operators must split without
splitting error; `strang` must reach second order and `advsrc` and `srcadv`
first order; with time-dependent callbacks, each scheme must reproduce its
exact split solution; source and boundary callbacks must see and return
moments in the documented order and stage, and source callbacks the moments
they request beyond the tracked ones from the quadrature of the current stage;
adaptive QMOM node counts must change as documented; realizability policies
and retries must govern the source step; states created from moments the host
changed must hold the quadrature of those moments, after any realizability
check; and failures must report their status and failed stage.

The error is
$\max_\alpha |M_\alpha - M_\alpha^{\mathrm{ref}}| / \max(1, |M_\alpha^{\mathrm{ref}}|)$
over the state's moments, the moments reconstructed from its quadrature, or
both, as each case states. The observed order is
$\log_2(e_{\Delta t} / e_{\Delta t/2})$.

QMOM references are closed forms of the configured or chosen moments. A CQMOM
case keeps only the moments its conditioning order needs (the loader reduces the
fixture's `4 x 4` tensor), so most CQMOM references apply the exact operators to
the quadrature `compute_cqmom` recovers from it; `commuting_cqmom_2d` and the
unit-invariance case start from the configured moments or weighted points.

Unless stated otherwise, the fixture cases advect nodes with RK23 at
tolerances `1e-12` (relative) and `1e-14` (absolute) and maximum step `0.01`,
the cases without a fixture with the default node-advection options, and the
moment source takes one SSP-RK3 step per macro step.

## Cases

Test names omit `verification_phase_space_evolution_` except the three
`verification_qmom_*` tests; "same" repeats the test name, and `<case>`
stands for the case name.

| Test | Case | Fixture |
| --- | --- | --- |
| `commuting_qmom` | same | `<case>.qmom` |
| `commuting_cqmom_2d` | same | `<case>.qmom` |
| `commuting_cqmom_state_adaptive_2d` | same | `<case>.qmom` |
| `noncommuting_convergence` | `noncommuting_qmom` | `<case>.qmom` |
| `noncommuting_cqmom_convergence` | `noncommuting_cqmom` | `noncommuting_cqmom_2d.qmom` |
| `verification_qmom_phase_space_source_evolution_matches_manufactured_solution` | `phase_space_source` | none |
| `verification_qmom_advection_matches_exponential_phase_evolution` | `direct_advection` | none |
| `verification_qmom_node_source_matches_exact_translation` | `node_source` | none |
| `cqmom_state_is_unit_invariant` | same | `unit_invariance_cqmom_state_adaptive_{reference,rescaled}.qmom` |
| `cqmom_native_source` | same | `<case>.qmom` |
| `cqmom_combined_native_source` | same | `<case>.qmom` |
| `qmom_boundary_contribution` | same | `commuting_qmom.qmom` |
| `cqmom_boundary_contribution` | same | `commuting_cqmom_2d.qmom` |
| `qmom_phase_space_evolution_failures` | same | `<case>.qmom` |
| `qmom_phase_space_evolution_inversion_failures` | same | `<case>.qmom` |
| `cqmom_phase_space_evolution_failures` | same | `<case>.qmom` |
| `cqmom_case_phase_space_evolution_failures` | same | `cqmom_phase_space_evolution_failures.qmom` |
| `cqmom_phase_space_evolution_inversion_failures` | same | `<case>.qmom` |
| `cqmom_native_source_failures` | same | `<case>.qmom` |

`qmom_phase_space_evolution_failures.qmom` repeats `commuting_qmom.qmom`. The
non-adaptive fixtures with conditioning order `{0,1}` hold the moments of the
product of $\xi_0 \in \{1, 3\}$ (weights $0.4, 0.6$) and
$\xi_1 \in \{-1, 2\}$ (weights $0.25, 0.75$). The inversion-failure fixtures
differ from their failure counterparts only in $m_0 = -1$ or $M_{0,0} = -1$.

The cases of `time_dependent_splitting.cpp`, `adaptive_qmom_evolution.cpp`,
`source_step_retries.cpp`, `state_candidate.cpp` and
`source_context_moments.cpp` build their states in the test, without a
fixture, and are listed in their sections below.

## Commuting operators

`commuting_qmom` uses $\dot\xi = 0.2\,\xi$ and $dm_k/dt = -0.05\,k\,m_k$.
Both act diagonally, so $m_k(t) = e^{0.15\,k\,t} m_k(0)$ and Strang splitting
adds no error. It runs one step of `0.08`, six of `0.04`, and the refinement
$\Delta t = 0.12, 0.06, 0.03$ to $T = 0.24$. Beforehand, `advance_nodes`
alone must keep the weights and scale the abscissas by $e^{0.2t}$, and
`advance_qmom_moments_with_source` alone must match the decay.

`commuting_cqmom_2d` uses $\dot\xi_0 = 0.16\,\xi_0$,
$\dot\xi_1 = -0.11\,\xi_1$ and $dM_{ij}/dt = -0.025\,(i+j)\,M_{ij}$, so
$M_{ij}(t) = e^{(0.135\,i - 0.135\,j)\,t} M_{ij}(0)$, with conditioning order
`{1,0}`. Each of five steps of `0.04` rebuilds the state from the moments
through `CqmomCase::initial_state()`.

`commuting_cqmom_state_adaptive_2d` applies the same advection with a zero
source to a `CqmomPreparedCase` with adaptive node counts, capacity
`2 x 2` and conditioning order `{0,1}`. The initial quadrature must be ragged
with fewer than four nodes; the fixture needs three. After five steps of
`0.04` it must still be ragged, with at most four nodes and one terminal path
per node, and its layout moments must match the advection reference.

## Noncommuting operators: splitting order

`noncommuting_qmom` evolves $m_0, \dots, m_3$ under $\dot\xi = 0.45$, which
gives $dm_k/dt = 0.45\,k\,m_{k-1}$, and $dm_k/dt = -0.17\,k\,m_k$. The
reference solves this bidiagonal system in closed form as
$m_k(t) = \sum_{r \le k} c_{k,r}\, e^{-0.17\,r\,t}$. Each scheme runs to
$T = 1.2$ with $\Delta t = 0.3, 0.15, 0.075, 0.0375, 0.01875$. The finest run
must conserve $m_0$ to `1e-12`, and `advsrc` and `srcadv` must differ by more
than `1e-4` at $\Delta t = 0.3$, showing that the scheme option takes effect.

`noncommuting_cqmom` uses a `CqmomPreparedCase`,
$\dot\xi_0 = \xi_0$, $\dot\xi_1 = 0$ and $dM_{ij}/dt = j\,M_{i+1,j-1}$, the
moment form of $\dot\xi_1 = \xi_0$. The weighted points follow
$\xi_0(t) = e^t \xi_0(0)$ and $\xi_1(t) = \xi_1(0) + (e^t - 1)\,\xi_0(0)$. The
rate needs moments outside the stored set, so the source callback
reconstructs a quadrature from its context once per step and evaluates the
three SSP-RK3 stage rates on the resulting `5 x 4` tensor. The plan matches
the QMOM case; the error is the larger of the layout and `4 x 4` tensor
errors, and `advsrc` and `srcadv` must differ by more than `1e-5`.

`phase_space_source` evolves $m_0, \dots, m_5$ of the weighted points
$\xi = 0.5, 1.5, 3$ (weights $0.2, 0.3, 0.5$) under $\dot\xi = \xi$ and a
source adding `0.1` to every moment, so $m_0(t) = m_0(0) + 0.1\,t$ and
$m_k(t) = m_k(0)\,e^{kt} + 0.1\,(e^{kt} - 1)/k$. Each scheme runs to
$T = 0.1$ with $\Delta t = 0.02$ down to `0.0025`, using RK45 node advection
at step `0.0025`; the error must fall at every refinement.

`time_dependent_splitting_<case>` (`time_dependent_splitting.cpp`) make both
operators depend on time. The weighted points $\xi = 0.5, 2$ (weights
$0.4, 0.6$) move under $\dot\xi = \xi/(1+t) + (1+t)\cos t$, so
$\xi(t) = (1+t)\,(\xi(0) + \sin t)$. The first term is the velocity; the
second is the drift $dm_k/dt = k\,(1+t)\cos t\; m_{k-1}$. Each scheme runs to
$T = 1$ with 10, 20, 40 and 80 steps and RK45 node advection at `1e-13`
relative.

- `schemes_match_their_exact_split_solutions`: each operator alone has a
  closed-form flow over $[s, s+h]$, $\xi \mapsto \xi\,(1+s+h)/(1+s)$ and
  $\xi \mapsto \xi + F(s+h) - F(s)$ with $F(t) = (1+t)\sin t + \cos t$.
  Composed in the order and on the intervals of each scheme, they give its
  exact split solution, which `strang`, `advsrc` and `srcadv` must reproduce
  up to the third-order error of the SSP-RK3 source step. A stage on the wrong
  interval adds a deviation of order $\Delta t$; for Lie splitting it leaves
  the order at one, so only this comparison detects it.
- `strang_is_second_order`: against the exact solution, the error must fall
  at every refinement and reach second order, which also checks the Strang
  split model; the Lie split solutions are first order by construction.

## Adaptive QMOM node counts

`adaptive_qmom_{grows_from_one_to_two_nodes,drops_a_vanishing_node}`
(`adaptive_qmom_evolution.cpp`) evolve an adaptive QMOM state with at most two
nodes, zero velocity and a constant moment rate over four steps of `0.25`.
The weights move linearly between two rules on the nodes $0.5$ and $1.5$, and
SSP-RK3 integrates the constant rate exactly, so after each step the reference
is the rule with weights $(1-t)\,w_a + t\,w_b$:

- growth, from $1$ at $1.5$ to $0.5$ at each node: one node initially, two
  after every step;
- reduction, the reverse: two nodes until the weight at $0.5$ vanishes at
  $t = 1$, then one node at $1.5$ with weight $1$.

The active count after each step must match, and the nodes and weights must
equal the reference.

## Source-step policies and retries

`source_retries_<case>` (`source_step_retries.cpp`) take one Strang step of
`1.5` for the moments $(1, 0, m_2, 0)$ of two points at $\pm\sqrt{m_2}$,
with $m_2 = 1$ and $dm_2/dt = -r\,m_2$ under a real-line policy. An SSP-RK3
step of size $\Delta t$ multiplies $m_2$ by $P(-r\,\Delta t)$, and its first
stage gives $m_2\,(1 - r\,\Delta t)$. For $r = 1$ the first stage of the full
step is $-0.5$ while $P(-1.5) = 0.0625$; for $r = 2$ the full step gives
$P(-3) = -2$, and each half step has the first stage $-0.5$ and the result
$P(-1.5)$.

| Case | Policy and options | Required result |
| --- | --- | --- |
| `final_only_policy_checks_only_the_step_result` | `FinalOnly`, $r = 2$, `max_retries = 1` | result rejected, two steps of `0.75` whose negative first stages are ignored, $m_2 = P(-1.5)^2$, one retry |
| `every_stage_policy_halves_a_rejected_step` | `EveryStage`, $r = 1$, `max_retries = 1` | first stage rejected, two steps of `0.75`, $m_2 = P(-0.75)^2$, one retry |
| `rejected_step_without_retries_fails` | `EveryStage`, $r = 1$ | `moment_source_failed`, `rejected_by_predicate` |
| `minimum_time_step_stops_halving` | `EveryStage`, $r = 1$, `max_retries = 4`, `min_time_step = 1` | `moment_source_failed`, `minimum_timestep_reached` |

The accepted cases use zero velocity and start at $t = 0$; the failing ones
use $\dot\xi = 0.1\,\xi$ and start at $t = 1$, so that their first advection
half-step changes the moments. An accepted step must also give the nodes
$\pm\sqrt{m_2}$ with weights $0.5$; a failed one must return its input
moments, not those after the first half-step, and its start time. The cases
check the first stage and the result of a step; the second stage stays
positive.

## States from changed moments

`state_candidate_<case>` (`state_candidate.cpp`) create states with
`QmomCase::state_candidate(...)` and `CqmomPreparedCase::state_candidate(...)`
from the moments of chosen weighted points, as a conservative host update
produces when it mixes populations. The cases are set up with other moments,
which a candidate does not use; the QMOM cases are adaptive with at most two
nodes, the CQMOM cases have two nodes per coordinate.

| Case | Moments | Required result |
| --- | --- | --- |
| `qmom_state_holds_the_quadrature_of_the_moments` | $0.4$ at $0.5$, $0.6$ at $1.5$, with $m_4$ and $m_5$ changed by $+10$ and $-10$ percent | that rule; the stored $m_4$ and $m_5$ are its moments |
| `qmom_policy_support_decides_the_verdict` | the same, $m_0$ to $m_3$ | accepted on the real line and on $[0, 2]$; rejected on $[0, 1]$, where $m_1 \le m_0$ is necessary |
| `qmom_state_keeps_the_maximum_node_count` | one point at $0.5$ | one node, then two after a step of `0.25` with zero velocity and the constant rate towards $0.5$ at each of $0.5$ and $1.5$: $0.875$ and $0.125$ |
| `cqmom_state_holds_the_quadrature_of_the_moments` | $0.15$ at $(0.5, 1)$, $0.2$ at $(0.5, 3)$, $0.3$ at $(1.5, 0.5)$, $0.35$ at $(1.5, 2)$ | those four points |
| `cqmom_policy_names_the_rejected_coordinate` | the same | checks on $[0, 2]$ for both coordinates: rejected at coordinate 1, where $2m_1 - m_2 = -0.225$; a check on $\xi_0$ only: accepted |

A rejection must leave `reconstruction_status` unset, since the check runs
before the reconstruction.

## Exact single operators

`direct_advection` advects the same weighted points under $\dot\xi = \xi$
with 400 SSP-RK3 steps of `2.5e-4`. The abscissas must equal
$e^{0.1}\xi(0)$, the weights must be bitwise unchanged, and $m_0, \dots,
m_{10}$ rebuilt from the quadrature must equal $e^{0.1k} m_k(0)$.

`node_source` moves the weighted points $0.5, 2$ (weights $0.4, 0.6$) at
$\dot\xi = 0.75$ through a source written with
`QmomSourceContext::for_each_node`, so $dm_k/dt = 0.75\,k\,m_{k-1}$. Each stage
quadrature reproduces $m_0, \dots, m_3$ and SSP-RK3 integrates this nilpotent
system exactly, so ten source-only steps of `0.01` must translate the points
by $0.75\,t$ up to roundoff.

`cqmom_state_is_unit_invariant` runs the weighted points $(-1, -1; 0.1)$,
$(1, -1; 0.3)$ and $(0.5, 2; 0.6)$ in units $(s_0, s_1)$ = `(1, 1)` and
`(1e6, 1e-6)`, with adaptive capacity `2 x 2` and conditioning order `{1,0}`.
One Strang step of `0.2` applies the velocity $(0.4, -0.25)$ times the units
with SSP-RK3 node advection, exact for a constant velocity, and
$dM_\alpha/dt = 0.1\,M_\alpha$, which multiplies the weights by $P(0.02)$, where
$P(z) = 1 + z + z^2/2 + z^3/6$ is the SSP-RK3 polynomial. Terminal paths,
weights, coordinates divided by the units and the `4 x 4` moments $M_{ij}$
divided by $s_0^i s_1^j$ are checked, and the source callback confirms through
`for_each_terminal_node` that it sees coordinates in the fixture's units after
the first half step.

## CQMOM source callbacks

`cqmom_native_source` takes one Strang step of `0.08` with zero velocity and
the rate $\gamma M_\alpha + B_\alpha$, $\gamma = -0.12$, where $B_\alpha$ are
the moments of the weight rates $0.02, -0.01, 0.03, -0.015$ at the four
fixture points. The reference is one SSP-RK3 step,
$P(z)\,M_\alpha + \Delta t\,(1 + z/2 + z^2/6)\,B_\alpha$ with
$z = \gamma\,\Delta t$. The callback checks the context size, hierarchy and
index order `{0,0}, {1,0}, {2,0}, {3,0}, {0,1}, {0,2}, {0,3}, {1,1}, {1,2},
{1,3}`, and on its first call that the values are the initial moments.

`cqmom_combined_native_source` adds the advection of `commuting_cqmom_2d` and
$dM_\alpha/dt = -0.09\,M_\alpha$ for one Strang step of `0.08`, with reference
$e^{(0.16\,i - 0.11\,j)\,\Delta t} P(-0.09\,\Delta t)\, M_{ij}(0)$. The
callback's first context must hold the moments after the first advection half
step.

## Untracked moments in source callbacks

`source_context_moments_<case>` (`source_context_moments.cpp`) write a node
motion as a moment source that needs moments the state does not track. The
velocity is zero, so the weights stay constant and the nodes follow closed
forms. The context must evaluate those moments from the quadrature of the
current SSP-RK3 stage; the quadrature at the start of the step makes the step
first order.

- `qmom_untracked_moment`: $0.4$ at $0.5$ and $0.6$ at $1$ to $T = 0.5$ under
  $\dot\xi = \xi^2$, written as $dm_k/dt = k\,m_{k+1}$ with `moment(k + 1)`,
  so $\xi(t) = \xi(0)/(1 - \xi(0)\,t)$; the rate of $m_3$ needs $m_4$.
- `qmom_generalized_moment`: $0.4$ at $0.5$ and $0.6$ at $2$ to $T = 1$ under
  $\dot\xi = \sqrt{\xi}$, written as $dm_k/dt = k\,m_{k-1/2}$ with
  `generalized_moment`, so $\xi(t) = (\sqrt{\xi(0)} + t/2)^2$.
- `cqmom_untracked_moment`, `cqmom_generalized_moment`: the four points of the
  CQMOM state-candidate cases to $T = 0.5$, with $\xi_0$ constant and
  $\dot\xi_1 = \xi_0^2$ or $\sqrt{\xi_0}$, written as
  $dM_{ij}/dt = j\,M_{i+2,j-1}$ or $j\,M_{i+1/2,j-1}$. The layout holds
  $M_{i0}$ for $i \le 3$ and $M_{ij}$ for $i \le 1$, $j \le 3$, so $M_{21}$,
  $M_{22}$, $M_{31}$ and $M_{32}$ are untracked.

The QMOM cases run 20, 40, 80 and 160 steps with a fixed count of two nodes,
and the error must fall at the third order of the SSP-RK3 source step. The
CQMOM cases are exact: the $\xi_0$ nodes do not move, each stage quadrature
reproduces the conditional moments of orders 0 to 3, and the rates need orders
up to 2, so the moment equations are linear and nilpotent and SSP-RK3 equals
the exact flow. Four steps of `0.125` must match the reference up to roundoff.

## Boundary contributions

`qmom_boundary_contribution` uses zero velocity, `srcadv`, injection at
$\xi = 0.25$ at rate `0.08` and removal at $\xi = 1.1$ at rate `-0.03`.
SSP-RK3 integrates a constant rate exactly, so one step of `0.04` must add
the step times the rate; a second run of `0.03` adds the source
$0.01\,(k+1)$. A wrong-sized, NaN or infinite boundary rate must fail in the
moment-source stage with a "boundary contribution" message, and a rate of
`-100` on $m_0$ must fail the post-source inversion.

`cqmom_boundary_contribution` uses a `CqmomPreparedCase` with
conditioning order `{1,0}`, zero velocity and one Strang step of `0.02`. The
rate is the moments of the weight rates $0.01\,(n+1)$ at the initial nodes, so
the reference is the initial `4 x 4` tensor plus the step times the rate. The
callback checks the context size, hierarchy and index order.

## Failure propagation

Except for the inversion failures and the `CqmomCase` policy probe, each case
takes one Strang step of `0.08`.

| Trigger | Required result |
| --- | --- |
| Velocity of the wrong length, or non-finite | `advection_failed` at `first_advection_half_step`, `velocity_failed`, message kept, final time `0`, no later stage completed; QMOM and `CqmomCase` states return unchanged moments |
| RK23 step `0.04`, minimum `0.02`, $\dot\xi = \xi$ (CQMOM: $(\xi_0, -\xi_1)$) | success with rejected steps in the first half step, endpoints `0.04` and `0.08`, moments matching the advection reference |
| Source rate `-20` on $m_0$ or $M_{0,0}$ | source step succeeds, `inversion_failed` at `post_source_inversion`, message names $m_0$ or realizability, no final reconstruction |
| Wrong-sized, NaN or infinite CQMOM source rate | `moment_source_failed` at `moment_source_step`, `source_failed`, message kept, no post-source reconstruction |
| $m_0 = -1$ or $M_{0,0} = -1$ | `initial_state()` throws; for CQMOM, `state_candidate(...)` reports `reconstruction_failed` and a realizability message |

The QMOM and prepared-case failure tests cover the first three rows.
`cqmom_case_phase_space_evolution_failures` covers the first and third rows
for a `CqmomCase`, adds a wrong-sized source rate, and, in a zero-rate step
of `0.02`, checks that basic diagnostics, the dense Vandermonde solver and
adaptive node counts set on the case reach the post-source
reconstructions. `cqmom_native_source_failures` covers the fourth
row. The CQMOM inversion-failure test reads its moments through an internal
accessor.

## Tolerances

| Case | Bound | Observed | Why |
| --- | --- | --- | --- |
| `commuting_qmom` | `2e-6`, `8e-6` in the refinement; order $\ge 2.4$; component checks: weights `1e-13`, abscissas `1e-10`, source alone `2e-9` | `1e-9` to `1e-7`; order 3 | safety bounds; the order shows the third-order SSP-RK3 source step dominates |
| `commuting_cqmom_2d`, adaptive | `3e-6` | `1e-11`, `1e-13` | safety bounds |
| `noncommuting_qmom` | `8e-4` Strang, `3e-2` Lie; order $\ge 1.75, 1.85$ and $\ge 0.75, 0.85$ on the two finest pairs | `1e-4`, `2e-2` at $\Delta t = 0.3$; orders 1.97, 1.00 | the bounds hold the splitting error at the coarsest step; the order thresholds sit 7 to 25 percent below 2 and 1 |
| `noncommuting_cqmom` | order $\ge 1.8$, $\ge 0.8$ on the two finest pairs | 2.00; 0.97 to 1.03 | only the order is checked |
| `phase_space_source` | order $> 1.7$; Lie in $(0.7, 1.3)$ on the last two pairs | 2.0, 1.0 | an order from two finite steps is approximate |
| exact split solutions | order $\ge 2.8$ on the two finest pairs; deviation `1e-6` at 80 steps | 2.97 to 3.01; `1e-7` | the SSP-RK3 source error is third order; the bound leaves a factor of 10 |
| time-dependent Strang order | order $\ge 1.75, 1.85$ on the two finest pairs | 2.07, 2.03 | the thresholds of `noncommuting_qmom` |
| adaptive QMOM, source retries | nodes, weights and moments `1e-12` | below `1e-14`, `1e-16` | exact references; rounding of the step and the reconstruction |
| state candidates | weights, positions and moments `1e-12` | below `1e-14` | exact references; rounding of the reconstruction |
| untracked source moments, QMOM | order $\ge 2.8$ on the two finest pairs; error `1e-5` at 160 steps | 2.94 to 3.00 on those pairs; `2e-6`, `2e-8` | the SSP-RK3 source error is third order; the bound leaves a factor of 6 and lies far below the `7e-3` and `6e-2` of the step-start quadrature, which gives order 1 |
| untracked source moments, CQMOM | moments `1e-12` | below `1e-14` | exact reference; roundoff |
| `direct_advection` | positions `2e-12`, moments `2e-11` | `1e-13`, `1e-12` | SSP-RK3 node error |
| RK23 minimum step | `2e-6` QMOM, `3e-6` CQMOM | `1e-7` | RK23 error at the minimum step |
| `node_source` | moments `1e-13`, points and weights `1e-12` | below `1e-15` | the step is exact; roundoff |
| unit invariance | `2e-9` | `1e-15` | exact reference; safety bound. Positions and moments are compared after dividing by the units, so both runs meet the same relative bound |
| native source, combined | `5e-10`, `5e-9` | below `1e-15`, `1e-13` | safety bounds |
| QMOM, CQMOM boundary | `5e-11`, `5e-9` | `1e-15`, below `1e-15` | exact reference; safety bounds |

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_(phase_space_evolution_|qmom_advection_matches|qmom_phase_space_source|qmom_node_source)' --output-on-failure
```

See [Building and running](../../../docs/developer/testing.md#building-and-running).

The verifiers are silent on success; on failure they print the failed checks
and, for the analytical cases, a per-run error table to standard error.
