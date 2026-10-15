# Reference

This page lists the configuration keys, the C++ entry points and callbacks,
the statuses and exceptions, troubleshooting messages, and the output formats.

## Configuration file

A `.qmom` file holds `[sections]` of `key = value` lines.

- A comment runs from `#` or `;` to the end of the line, wherever it appears,
  so no value can contain either character.
- Numbers are plain decimals, such as `2`, `-0.5`, or `1.0e-6`, whatever the
  locale.
- Indices in section names have no leading zeros: `[initial.dimension.0]`.
- Unknown keys are errors, and so are keys that contradict the other settings,
  such as a parameter of another distribution family; both are reported with
  the file and line.
- Explicit `[moments]` take precedence over `[initial]`, which is then neither
  read nor checked.
- An optional top-level `format_version = 1` states the format version.

The case loaders check every key they read; the other loaders read only their
own sections:

| Function | Returns |
| --- | --- |
| `config::load_qmom(path)`, `config::load_cqmom_prepared(path)`, `config::load_cqmom(path)` | the case |
| `config::load_qmom_phase_space_evolution_options(path)`, `config::load_cqmom_phase_space_evolution_options(path)` | the evolution options |
| `config::load_output_options(path)` | the output options |
| `config::parse_qmom_initial_condition(text)`, `config::parse_cqmom_initial_condition(text)` | the `[initial]` distribution, for example to sample it |

Each `load_*(path)` has a `parse_*(text)` counterpart that reads the same from a
string.

### `[phase_space]`

| Key | Default | Meaning | C++ |
| --- | --- | --- | --- |
| `node_counts` | required | Maximum nodes per coordinate; one coordinate selects QMOM, several CQMOM. | `QmomCase::node_count`, `PhaseSpaceHierarchy` |
| `names` | none | One label per coordinate, for the reader; QMOMKit does not use them. | — |
| `conditioning_order` | `0, 1, …` | CQMOM conditioning order. | `PhaseSpaceHierarchy` |
| `support_lower`, `support_upper` | none | Interval the nodes must lie in, per coordinate; `none` for no bound. Both or neither, finite, `lower < upper`. | `abscissa_bounds` |

### `[moments]`

| Key | Meaning |
| --- | --- |
| `values` | QMOM: raw moments $m_0, m_1, \ldots$, at least $2N$. |
| `extents`, `M[i,j,…]` | CQMOM: tensor extents, at least twice each coordinate's node count, then one `M[i,j,…] = value` line for every entry within them. |

### `[initial]`

| Key | Default | Meaning |
| --- | --- | --- |
| `type` | required | `distribution` or `mixture` for QMOM, `independent_product` for CQMOM. |
| `family`, parameters | — | `type = distribution`: the family and its parameters, below. |
| `number_density` | `1.0` | `type = distribution`: the total, $m_0$. Required for `independent_product`. |
| `maximum_order` | $2N - 1$ | QMOM: highest moment order generated. |
| `maximum_orders` | — | CQMOM: highest order per coordinate, only for code that builds a full tensor with `initialize_cqmom_moment_tensor(...)`. |

- **Mixture:** sections `[initial.mixture.component.0]`, `[...1]`, …, each with
  `family`, its parameters, and a required `number_density`, which adds to the
  total.
- **Independent product:** sections `[initial.dimension.0]`, `[...1]`, …, one
  per coordinate in declared order, each with `family` and parameters, or
  `type = mixture` with components `[initial.dimension.N.mixture.component.M]`.
  Each coordinate's distribution is normalized; `number_density` scales the
  product.

### Distribution families

| `family` | Parameters | Constraints | C++ type |
| --- | --- | --- | --- |
| `normal` | `mean`, `standard_deviation` | `standard_deviation > 0` | `NormalDistribution` |
| `lognormal` | `mu_log`, `sigma_log` | `sigma_log > 0` | `LogNormalDistribution` |
| `beta` | `alpha`, `beta` | both `> 0`; support $[0, 1]$ | `BetaDistribution` |
| `gamma` | `shape`, `scale` | both `> 0` | `GammaDistribution` |
| `weibull` | `shape`, `scale` | both `> 0` | `WeibullDistribution` |
| `chi` | `degrees_of_freedom` | `> 0` | `ChiDistribution` |
| `chi-squared` | `degrees_of_freedom` | `> 0` | `ChiSquaredDistribution` |

`log-normal` and `log_normal` are accepted for `lognormal`, and `chi_squared`
and `chisquared` for `chi-squared`. All parameters must be finite.

### Reconstruction options

| Key | Default | Meaning | C++ |
| --- | --- | --- | --- |
| `[qmom] adaptive` | `true` | Use fewer nodes where the moments support no more. | `QmomOptions::adaptive` |
| — | `1e-10` | Adaptive: smallest accepted ratio of the smallest to the largest weight. | `QmomOptions::rmin` |
| — | `1e-10` | Adaptive: smallest accepted relative node separation. | `QmomOptions::eabs` |
| — | `1e-30` | Smallest $m_0$ that is reconstructed. | `QmomOptions::cutoff` |

### `[cqmom]` and `[vandermonde]`

| Key | Default | Meaning | C++ |
| --- | --- | --- | --- |
| `[cqmom] diagnostic_level` | `none` | `basic` records branch structure, `detailed` adds the relative residual and condition estimate of each conditional system, in `CqmomResult::diagnostics`. Results are unchanged. | `CqmomCase::diagnostic_level` |
| `[vandermonde] solver` | `automatic` | Solver of the conditional systems: `automatic` (Björck–Pereyra), `bjorck_pereyra`, `dense`, `pivoted_qr`, `svd`, `moore_penrose`. | `CqmomCase::vandermonde_solver` |
| `[vandermonde] rcond` | `1e-14` | Singular values at or below `rcond` times the largest are discarded by `moore_penrose`; the other solvers ignore it. `0 <= rcond < 1`. | `ConditionalMomentOptions::moore_penrose_rcond` |

### Evolution options

| Key | Default | Meaning | C++ (`QmomPhaseSpaceEvolutionOptions`, `CqmomPhaseSpaceEvolutionOptions`) |
| --- | --- | --- | --- |
| `[time_integration] scheme` | `strang` | Splitting: `strang`, `advsrc`, or `srcadv`. | `scheme` |
| `[time_integration.node_advection] scheme` | `rk23` | Node advection: `rk23`, `rk45`, or `ssp_rk3`. | `node_options.scheme` |
| `[time_integration.node_advection.adaptive]` `relative_tolerance`, `absolute_tolerance` | `1e-6`, `1e-9` | Error tolerances of `rk23` and `rk45`. | `node_options.adaptive_options` |
| `initial_step` | `1e-6` | First step tried. | 〃 |
| `max_step` | whole interval | Largest step. | 〃 |
| `min_step` | automatic | Smallest step, never below ten floating-point spacings of the elapsed time. A step at this floor is accepted even above the tolerance. | 〃 |
| `max_steps` | `1000000000` | Most accepted steps per interval. | 〃 |
| — | `0` | Source step: retries with half the step after a rejection. | `source_options.ssp_rk3_options.max_retries` |
| — | `0` | Source step: shortest retried step. | `source_options.ssp_rk3_options.min_time_step` |

### `[time_integration.realizability]`

For QMOM; a CQMOM file gives one section per checked coordinate,
`[time_integration.realizability.N]`.

| Key | Default | Meaning | C++ (`RealizabilityPolicy`) |
| --- | --- | --- | --- |
| `check_mode` | `none` | `none`, `final_only` (the result of the source step), or `every_stage` (also its stages). | `check_mode` |
| `support` | — | `real_line`, `nonnegative`, or `bounded_interval`. Required unless `check_mode = none`. | `support` |
| `lower`, `upper` | — | Required for `bounded_interval`, an error otherwise. | `interval` |
| `accept_boundary` | `true` | Accept moments on the boundary of the realizable set. | `accept_boundary` |
| `accept_indeterminate` | `false` | Accept moments that carry too few digits to decide (`numerically_indeterminate`). | `accept_indeterminate` |
| `tolerance` | `1e-12` | Relative tolerance of the check. | `options.tolerance` |

### `[output]`

| Key | Default | Meaning |
| --- | --- | --- |
| `directory` | `.` | For your code: where to write. |
| `write_moments`, `write_quadrature` | `false` | For your code: what to write. |
| `interval` | none | For your code: when to write; positive. |
| `generalized_moments` | none | Extra moment columns: orders separated by `\|`, for example `0.5 \| 1.5`, or for CQMOM `0.5, 0 \| 1, 1.5`. Non-negative. |

## C++ entry points

`#include <qmomkit/qmomkit.hpp>` includes everything. The headers document each
function.

| Task | Functions | Header |
| --- | --- | --- |
| Load a configuration | `config::load_qmom`, `load_cqmom_prepared`, `load_cqmom`, `load_*_options`, their `parse_*` forms, `with_initial_moments` | `configuration.hpp` |
| Generate moments | `compute_raw_moments` | `distributions.hpp` |
| Set up CQMOM moments | `create_cqmom_moment_layout`, `extract_cqmom_moments`, `create_cqmom_moment_data` | `cqmom.hpp` |
| Reconstruct | `compute_qmom`, `compute_cqmom` | `qmom.hpp`, `cqmom.hpp` |
| Create a state | `initial_state()`, `state_candidate(...)` | `qmom.hpp`, `cqmom.hpp` |
| Read a state | `moment_values()`, `moment_indices()`, `moment(index)`, `quadrature()`, `active_node_count()`, `active_branch_counts()` | `qmom.hpp`, `cqmom.hpp` |
| Advance a state | `advance_{qmom,cqmom}_phase_space`, `_advection`, `_moment_source`, and their `_batch` and `_batch_with_index` forms | `phase_space_evolution.hpp` |
| Check realizability | `check_realizability` | `realizability.hpp` |
| Derived moments | `evaluate_generalized_moment`, `reconstruct_raw_moments` | `moment_reconstruction.hpp` |
| Write output | `MomentWriter`, `QuadratureWriter` | `output.hpp` |
| Without a state | `advance_nodes`, `advance_qmom_moments_with_source`, `advance_cqmom_moments_with_source` | `node_advection.hpp`, `moment_source.hpp` |
| Integrate an ODE | `ssp_rk3_step`, `integrate_embedded_rk` | `runge_kutta.hpp` |

## Callbacks

| Type | Signature |
| --- | --- |
| `PhaseSpaceVelocity` | `void(double time, const std::vector<double>& coordinates, std::vector<double>& velocities, std::size_t dimensions)` |
| `PerNodePhaseSpaceVelocity` | `std::vector<double>(double time, const std::vector<double>& coordinate)`, adapted by `make_phase_space_velocity` |
| `QmomMomentSource`, `QmomMomentBoundaryContribution` | `QmomMomentRate(double time, const QmomSourceContext& context)` |
| `CqmomMomentSource`, `CqmomMomentBoundaryContribution` | `CqmomMomentRate(double time, const CqmomSourceContext& context)` |

The `Batch` forms, such as `BatchPhaseSpaceVelocity`, take a leading
`std::size_t batch_index`. For driving the integrators directly:

| Type | Signature |
| --- | --- |
| `OdeRhs` | `void(double time, const std::vector<double>& state, std::vector<double>& derivative)` |
| `StateAcceptancePredicate` | `bool(const std::vector<double>& state)`, set in `ssp_rk3_options.accept` |
| `StageAcceptancePredicate` | `bool(SspRk3Stage stage, const std::vector<double>& state)`, set in `ssp_rk3_options.stage_accept` |

## Statuses

| Result | Status | Meaning |
| --- | --- | --- |
| `QmomResult` | `invalid_input` | Too few or non-finite moments, no nodes, or invalid options or bounds. |
| | `not_realizable` | The moments give no quadrature with the requested nodes. |
| `CqmomResult` | `invalid_input` | The moments do not fit the hierarchy, or invalid options. |
| | `qmom_failed` | A one-dimensional reconstruction of a branch failed. |
| | `conditional_moment_failed` | A conditional system could not be solved. |
| Evolution results | `invalid_input` | Invalid state, time, callbacks, or options; nothing ran. |
| | `inversion_failed` | The post-source reconstruction failed. |
| | `advection_failed` | The velocity or the node integration failed. |
| | `moment_reconstruction_failed` | The moments of the moved or reconstructed nodes do not fit in a double. |
| | `moment_source_failed` | The source, or its integration or check, failed. |
| `QmomStateCandidate`, `CqmomStateCandidate` | `empty` | $0 \le m_0 <$ cutoff, whatever a policy says: nothing to evolve. |
| | `realizability_rejected` | The policy rejected the moments. |
| | `reconstruction_failed` | No quadrature from these moments. |
| | `invalid_input` | Non-finite moments, too few, a wrong layout, or an invalid policy. |
| `QmomNodeAdvectionResult`, `CqmomNodeAdvectionResult` | `invalid_input` | An invalid quadrature, time, callback, or options. |
| | `velocity_failed` | The velocity threw, or returned the wrong size or a non-finite value. |
| | `integration_failed` | The node integration failed. |
| `QmomMomentSourceUpdateResult`, `CqmomMomentSourceUpdateResult` | `invalid_input` | Invalid moments, time, callbacks, or options. |
| | `source_failed` | The source or boundary contribution failed. |
| | `rejected_by_predicate`, `minimum_timestep_reached` | A rejected step could not be retried. |
| `SspRk3Result` | `invalid_input`, `rejected_by_predicate`, `minimum_timestep_reached` | As for a moment-source update. |
| `AdaptiveRkResult` | `invalid_input` | Invalid options or input. |
| | `step_size_underflow`, `max_steps_exceeded` | A non-finite state at the minimum step, or more than `max_steps` steps. |
| `RawMomentResult` and the other moment-generation results | `invalid_input` | Invalid distribution parameters or orders. |
| | `overflow`, `underflow`, `non_finite` | A moment does not fit in a double. |
| `RealizabilityResult` | `strictly_realizable`, `boundary` | Realizable; `realizable()` is true. |
| | `not_realizable` | Not realizable on the support. |
| | `numerically_indeterminate` | The moments carry too few digits to decide. |
| | `invalid_input` | Empty or non-finite moments, an invalid tolerance, or an invalid support or interval. |

Every result type in this table also has a `success` status and `ok()`, except
`RealizabilityResult`; `AdaptiveRkResult` has `success` but no `ok()`. An
evolution result's `diagnostics.failed_stage` names the part of the step that
failed: `first_advection_half_step`, `first_moment_reconstruction`,
`moment_source_step`, `post_source_inversion`, `second_advection_half_step`, or
`final_moment_reconstruction`.

These throw instead of returning a status:

| Function | Throws |
| --- | --- |
| `initial_state()` | `std::runtime_error` when the moments of the case cannot be reconstructed or any of them is not finite; for a `CqmomCase` also `std::invalid_argument` for a layout that is not the required layout of its hierarchy |
| `PhaseSpaceHierarchy` constructor | `std::invalid_argument` for no node counts, a zero node count, an invalid conditioning order, or too many nodes in total |
| `CqmomPreparedCase` constructor | `std::invalid_argument` for a layout that is not the required layout of its hierarchy; wrong values show up when the case is used |
| `extract_cqmom_moments(...)` | `std::invalid_argument` when the tensor's rank or value count is wrong; `std::out_of_range` when it lacks a required moment |
| `create_cqmom_moment_data(...)` | `std::invalid_argument` when the result failed, or its indices, count or values do not match the layout |
| `reconstruct_raw_moments(...)`, `reconstruct_cqmom_moments(...)`, `evaluate_generalized_moment(...)` | `std::invalid_argument` for malformed input, `std::domain_error` for a fractional order with a negative abscissa, `std::overflow_error` for a value beyond the range of a double |
| `weight(index)`, `abscissa(dimension, index)` of an adaptive `CqmomQuadrature` | `std::logic_error`; iterate over `weights()` and `abscissas()` instead |
| `MomentWriter`, `QuadratureWriter` | `std::runtime_error` when a file cannot be written; `MomentWriter` also `std::domain_error` for a fractional generalized order with a negative abscissa |
| `value()` of a failed configuration result | `std::bad_optional_access` |
| `ssp_rk3_step(...)`, `integrate_embedded_rk(...)` | exceptions from your callbacks, unchanged; `std::invalid_argument` when `rhs` changes the size of its buffer, returns a non-finite value or, in a Debug build, leaves an entry unwritten, and from `ssp_rk3_step` also when a candidate state is not finite |

## Troubleshooting

| Message | Cause and remedy |
| --- | --- |
| `Moment sequence has fewer than 2*n entries.` | $N$ nodes need $2N$ moments. Supply more, request fewer nodes, or raise `maximum_order`. |
| `Moment sequence has negative recurrence b = …` | The moments are not realizable. Check them first: mixing normalized and raw moments is a common cause. `check_realizability(...)` names the failing matrix. |
| `Moment sequence supports at most k nodes …` | A fixed reconstruction asked for more nodes than the moments carry. Use adaptive mode or fewer nodes. |
| `Moment sequence has m0 < cutoff.` | There is no population left; stop evolving it. |
| `… returned nodes outside the configured interval.` | The moments imply nodes outside the bounds: the bounds or the moments are wrong. |
| `Wheeler recurrence produced non-finite values.` | The moments are near the range of a double. Use units that keep coordinates near unity. |
| `Phase-space velocity callback returned a non-finite value …` | The message names the node; a velocity dividing by a coordinate that reached zero is typical. |
| `… exceeded max_steps.` | The velocity is stiffer than expected, or `max_steps` is too low. Check the velocity first. |
| `… non-finite state or error estimate at the minimum step.` | The node advection diverged; check the velocity. Relaxing tolerances does not help. |
| `inversion_failed` after many good steps | The moments drift out of the realizable set. Use adaptive mode, a realizability policy with retries, or a smaller step. |
| `Moment-source update was rejected by the acceptability predicate.`, `… reached the minimum timestep before an acceptable state was found.` | The source drives the moments out of the realizable set. Allow retries (`max_retries`), use a smaller step, or treat it as a physics problem. |
| `… must match the source context size.`, `… must match the moment state size.`, `RHS returned a vector with the wrong size.`, `… must be finite.` | Return exactly `context.size()` finite values, in context order. |

For more detail on a CQMOM failure, set `[cqmom] diagnostic_level = detailed`
and read `CqmomResult::diagnostics`, or
`diagnostics.post_source_cqmom_diagnostics` of an evolution result.

## Output files

Writers produce text tables: one header line, `# time` and a label per column,
then columns separated by whitespace, in full double precision.

| Writer | Columns |
| --- | --- |
| `MomentWriter` | `time`, then one column per moment: `M[0] M[1] …` for QMOM, `M[0,0] M[1,0] …` for CQMOM, in the order of `moment_indices()`, then generalized moments, such as `G[0.5]`. |
| `QuadratureWriter`, QMOM | `time node_index weight x0`, one row per node. |
| `QuadratureWriter`, CQMOM | `time node_index path_0 … weight x0 …`, one row per terminal node; `path_0 …` are the node's index at each conditioning level. |

`x0 …` are the node's abscissas, one per coordinate. An adaptive quadrature can
have a different number of rows at different times; group rows by `time`.
