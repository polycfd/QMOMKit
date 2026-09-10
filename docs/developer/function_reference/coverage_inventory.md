# Definition coverage and navigation inventory

This inventory routes core implementation definitions to functional explanations.
It complements the algorithm sections and responsibility tables; a symbol listing
does not by itself establish intent or test coverage.

## Scope and method

The inventory covers the 23 core `.cpp` files and installed/private headers.
It lists explicit functions, methods, constructors, and operators by defining
source range. Defaulted or implicit special members are omitted; consequential
lambdas are documented with their owning function. Instrumentation definitions
are included for its enabled preprocessor branch. Line links are navigation aids,
not permanent API identifiers.

The inventory includes **940 definition locations**, grouped by named function
within each defining file. Overloads and separately compiled instrumentation
branches share a row with distinct line links. The count includes trivial helpers
and is not a count of numerical algorithms. Detailed pages cover the important
algorithm families; entry, accessor, diagnostic, parsing, formatting, and
forwarding functions are grouped by responsibility. No coverage percentage or
exhaustive branch-coverage claim is made.

The core numerical backbone receives Level A treatment: Wheeler and scaling;
support-specific realizability; weighted conditional systems and numerical
kernels; required/conditional layouts and recursive branch policy; index/topology
transformations; raw/generalized moments and rates; analytical initialization;
stage-local source closure; SSP/embedded RK; and accepted evolution transitions.
Configuration policy resolution and public workflows combine Levels A/B; ordinary
parsing, ownership, getters, output, and counters receive lighter B/C coverage.
Specific tests and remaining gaps are attached to the functional explanations,
not inferred from proximity of filenames.

Tests, examples, configuration fixtures, package metadata, and build files are
linked where they clarify a contract, but this inventory is limited to the core
library. The optional `tools/monte_carlo/` companion, plotting scripts, benchmark
harnesses, and application physics in examples remain separate consumer/package
responsibilities. The inventory does not claim to document every function in
those consumers.

## Definition ledger

Each file heading links to its definitions and owning functional reference. Local
scope means an anonymous namespace; public/internal ownership and class contracts
are explained on the destination page. For out-of-class methods, the owner is
retained when present in the definition signature. Trivial defaulted ownership
operations are covered by the state/result helper tables, even when excluded from
the executable-body count.

### `include/qmomkit/conditional_moments.hpp`

[Source](../../../include/qmomkit/conditional_moments.hpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `ConditionalMomentResult::ok` | [80](../../../include/qmomkit/conditional_moments.hpp#L80) |

### `include/qmomkit/configuration.hpp`

[Source](../../../include/qmomkit/configuration.hpp) · [Responsibilities, contracts, and evidence](initialization_and_configuration.md)

| Function or overload family | Definition lines |
|---|---|
| `ConfigurationResult` | [46](../../../include/qmomkit/configuration.hpp#L46) |
| `failure` | [48](../../../include/qmomkit/configuration.hpp#L48) |
| `ok` | [55](../../../include/qmomkit/configuration.hpp#L55) |
| `operator bool` | [60](../../../include/qmomkit/configuration.hpp#L60) |
| `value` | [65](../../../include/qmomkit/configuration.hpp#L65), [75](../../../include/qmomkit/configuration.hpp#L75) |
| `ConfigurationResult::value` | [70](../../../include/qmomkit/configuration.hpp#L70) |
| `error` | [80](../../../include/qmomkit/configuration.hpp#L80) |

### `include/qmomkit/cqmom.hpp`

[Source](../../../include/qmomkit/cqmom.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `CqmomResult::ok` | [303](../../../include/qmomkit/cqmom.hpp#L303) |

### `include/qmomkit/distributions.hpp`

[Source](../../../include/qmomkit/distributions.hpp) · [Responsibilities, contracts, and evidence](initialization_and_configuration.md)

| Function or overload family | Definition lines |
|---|---|
| `RawMomentResult::ok` | [118](../../../include/qmomkit/distributions.hpp#L118) |
| `SelectedRawMomentResult::ok` | [154](../../../include/qmomkit/distributions.hpp#L154) |
| `TensorMomentResult::ok` | [168](../../../include/qmomkit/distributions.hpp#L168) |

### `include/qmomkit/moment_source.hpp`

[Source](../../../include/qmomkit/moment_source.hpp) · [Responsibilities, contracts, and evidence](moment_sources.md)

| Function or overload family | Definition lines |
|---|---|
| `MomentSourceUpdateResult::ok` | [235](../../../include/qmomkit/moment_source.hpp#L235) |
| `CqmomMomentSourceUpdateResult::ok` | [248](../../../include/qmomkit/moment_source.hpp#L248) |

### `include/qmomkit/node_advection.hpp`

[Source](../../../include/qmomkit/node_advection.hpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomNodeAdvectionResult::ok` | [66](../../../include/qmomkit/node_advection.hpp#L66) |
| `CqmomNodeAdvectionResult::ok` | [79](../../../include/qmomkit/node_advection.hpp#L79) |

### `include/qmomkit/phase_space_evolution.hpp`

[Source](../../../include/qmomkit/phase_space_evolution.hpp) · [Responsibilities, contracts, and evidence](phase_space_evolution.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomPhaseSpaceEvolutionResult::ok` | [139](../../../include/qmomkit/phase_space_evolution.hpp#L139) |
| `CqmomPhaseSpaceEvolutionResult::ok` | [154](../../../include/qmomkit/phase_space_evolution.hpp#L154) |

### `include/qmomkit/qmom.hpp`

[Source](../../../include/qmomkit/qmom.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomResult::ok` | [40](../../../include/qmomkit/qmom.hpp#L40) |
| `QmomState::moments` | [69](../../../include/qmomkit/qmom.hpp#L69) |
| `QmomState::quadrature` | [74](../../../include/qmomkit/qmom.hpp#L74) |
| `QmomState::active_node_count` | [79](../../../include/qmomkit/qmom.hpp#L79) |
| `QmomState::abscissa_bounds` | [84](../../../include/qmomkit/qmom.hpp#L84) |

### `include/qmomkit/quadrature.hpp`

[Source](../../../include/qmomkit/quadrature.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `Quadrature1D::ok` | [30](../../../include/qmomkit/quadrature.hpp#L30) |

### `include/qmomkit/realizability_policy.hpp`

[Source](../../../include/qmomkit/realizability_policy.hpp) · [Responsibilities, contracts, and evidence](realizability.md)

| Function or overload family | Definition lines |
|---|---|
| `accepts_realizability_policy` | [56](../../../include/qmomkit/realizability_policy.hpp#L56) |
| `make_realizability_acceptance_predicate` | [95](../../../include/qmomkit/realizability_policy.hpp#L95) |
| `make_realizability_stage_acceptance_predicate` | [113](../../../include/qmomkit/realizability_policy.hpp#L113) |

### `include/qmomkit/runge_kutta.hpp`

[Source](../../../include/qmomkit/runge_kutta.hpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `SspRk3Result::ok` | [64](../../../include/qmomkit/runge_kutta.hpp#L64) |

### `src/conditional_moments.cpp`

[Source](../../../src/conditional_moments.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `invalid_result` | [25](../../../src/conditional_moments.cpp#L25) |
| `valid_options` | [40](../../../src/conditional_moments.cpp#L40) |
| `extract_rhs` | [46](../../../src/conditional_moments.cpp#L46) |
| `ConditionalMomentResult::moment` | [71](../../../src/conditional_moments.cpp#L71) |
| `reconstruct_conditional_moments` | [80](../../../src/conditional_moments.cpp#L80) |

### `src/conditioning_path.cpp`

[Source](../../../src/conditioning_path.cpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `ConditioningPath::ConditioningPath` | [14](../../../src/conditioning_path.cpp#L14) |
| `ConditioningPath::empty` | [16](../../../src/conditioning_path.cpp#L16) |
| `ConditioningPath::size` | [21](../../../src/conditioning_path.cpp#L21) |
| `ConditioningPath::node_indices_by_level` | [26](../../../src/conditioning_path.cpp#L26) |
| `ConditioningPath::node_index_at_level` | [31](../../../src/conditioning_path.cpp#L31) |
| `ConditioningPath::with_child` | [40](../../../src/conditioning_path.cpp#L40) |
| `ConditioningPath::phase_space_node_index` | [47](../../../src/conditioning_path.cpp#L47) |
| `ConditioningPath::to_string` | [68](../../../src/conditioning_path.cpp#L68) |
| `ConditioningPath::describe` | [84](../../../src/conditioning_path.cpp#L84) |
| `operator==` | [102](../../../src/conditioning_path.cpp#L102) |
| `operator!=` | [107](../../../src/conditioning_path.cpp#L107) |
| `operator<` | [112](../../../src/conditioning_path.cpp#L112) |
| `recursive_path_to_node_index` | [118](../../../src/conditioning_path.cpp#L118) |

### `src/configuration.cpp`

[Source](../../../src/configuration.cpp) · [Responsibilities, contracts, and evidence](initialization_and_configuration.md)

| Function or overload family | Definition lines |
|---|---|
| `trim` | [54](../../../src/configuration.cpp#L54) |
| `add_error` | [65](../../../src/configuration.cpp#L65) |
| `parse_size` | [71](../../../src/configuration.cpp#L71) |
| `parse_coordinate_realizability_section` | [83](../../../src/configuration.cpp#L83) |
| `is_coordinate_realizability_key` | [92](../../../src/configuration.cpp#L92) |
| `is_configured_coordinate_realizability_key` | [107](../../../src/configuration.cpp#L107) |
| `parse_double` | [117](../../../src/configuration.cpp#L117) |
| `parse_size_list` | [131](../../../src/configuration.cpp#L131) |
| `parse_double_list` | [147](../../../src/configuration.cpp#L147) |
| `parse_optional_double_list` | [163](../../../src/configuration.cpp#L163) |
| `validate_names` | [185](../../../src/configuration.cpp#L185) |
| `parse_moment_key` | [206](../../../src/configuration.cpp#L206) |
| `mixture_component_prefix` | [227](../../../src/configuration.cpp#L227) |
| `is_mixture_component_section` | [232](../../../src/configuration.cpp#L232) |
| `is_independent_product_dimension_section` | [243](../../../src/configuration.cpp#L243) |
| `parse_file` | [253](../../../src/configuration.cpp#L253) |
| `require_entry` | [319](../../../src/configuration.cpp#L319) |
| `find_entry` | [330](../../../src/configuration.cpp#L330) |
| `reject_unknown_keys` | [340](../../../src/configuration.cpp#L340) |
| `reject_unknown_sections` | [358](../../../src/configuration.cpp#L358) |
| `validate_format_version` | [377](../../../src/configuration.cpp#L377) |
| `parse_phase_space_node_counts` | [387](../../../src/configuration.cpp#L387) |
| `validate_phase_space_names` | [399](../../../src/configuration.cpp#L399) |
| `validate_qmom_phase_space` | [409](../../../src/configuration.cpp#L409) |
| `validate_cqmom_phase_space` | [422](../../../src/configuration.cpp#L422) |
| `validate_optional_qmom_phase_space` | [442](../../../src/configuration.cpp#L442) |
| `validate_optional_cqmom_phase_space` | [458](../../../src/configuration.cpp#L458) |
| `parse_required_double` | [474](../../../src/configuration.cpp#L474) |
| `parse_optional_double` | [489](../../../src/configuration.cpp#L489) |
| `parse_optional_bool` | [504](../../../src/configuration.cpp#L504) |
| `parse_optional_size` | [525](../../../src/configuration.cpp#L525) |
| `parse_optional_positive_double` | [540](../../../src/configuration.cpp#L540) |
| `infer_qmom_initial_maximum_order` | [557](../../../src/configuration.cpp#L557) |
| `parse_optional_positive_size` | [584](../../../src/configuration.cpp#L584) |
| `parse_distribution` | [601](../../../src/configuration.cpp#L601), [688](../../../src/configuration.cpp#L688) |
| `component_key_allowed` | [693](../../../src/configuration.cpp#L693) |
| `marginal_key_allowed` | [705](../../../src/configuration.cpp#L705) |
| `has_explicit_qmom_moments` | [721](../../../src/configuration.cpp#L721) |
| `has_explicit_cqmom_moments` | [726](../../../src/configuration.cpp#L726) |
| `qmom_allowed_keys` | [731](../../../src/configuration.cpp#L731) |
| `append_time_integration_allowed_keys` | [756](../../../src/configuration.cpp#L756) |
| `append_output_allowed_keys` | [775](../../../src/configuration.cpp#L775) |
| `with_runtime_allowed_keys` | [783](../../../src/configuration.cpp#L783) |
| `time_integration_sections` | [790](../../../src/configuration.cpp#L790) |
| `reject_unknown_qmom_sections` | [795](../../../src/configuration.cpp#L795) |
| `reject_unknown_qmom_keys` | [808](../../../src/configuration.cpp#L808) |
| `parse_mixture` | [822](../../../src/configuration.cpp#L822) |
| `has_mixture_component_sections` | [890](../../../src/configuration.cpp#L890) |
| `cqmom_explicit_allowed_keys` | [902](../../../src/configuration.cpp#L902) |
| `cqmom_initial_allowed_keys` | [910](../../../src/configuration.cpp#L910) |
| `reject_unknown_cqmom_sections` | [919](../../../src/configuration.cpp#L919) |
| `reject_unknown_cqmom_keys` | [932](../../../src/configuration.cpp#L932) |
| `parse_independent_product` | [942](../../../src/configuration.cpp#L942) |
| `checked_product` | [1039](../../../src/configuration.cpp#L1039) |
| `parse_explicit_tensor` | [1050](../../../src/configuration.cpp#L1050) |
| `format_configuration_error` | [1112](../../../src/configuration.cpp#L1112) |
| `format_configuration_errors` | [1128](../../../src/configuration.cpp#L1128) |
| `failure_result` | [1137](../../../src/configuration.cpp#L1137), [1142](../../../src/configuration.cpp#L1142) |
| `file_failure` | [1147](../../../src/configuration.cpp#L1147) |
| `validate_replacement_tensor_structure` | [1152](../../../src/configuration.cpp#L1152) |
| `format_moment_index` | [1175](../../../src/configuration.cpp#L1175) |
| `read_configuration_file` | [1191](../../../src/configuration.cpp#L1191) |
| `has_realizability_section` | [1216](../../../src/configuration.cpp#L1216) |
| `has_unnumbered_realizability_section` | [1226](../../../src/configuration.cpp#L1226) |
| `coordinate_realizability_sections` | [1236](../../../src/configuration.cpp#L1236) |
| `is_realizability_key` | [1250](../../../src/configuration.cpp#L1250) |
| `realizability_policy_fields` | [1255](../../../src/configuration.cpp#L1255) |
| `is_node_advection_key` | [1261](../../../src/configuration.cpp#L1261) |
| `node_advection_allowed_keys` | [1266](../../../src/configuration.cpp#L1266) |
| `parse_qmom_options_from_parsed` | [1275](../../../src/configuration.cpp#L1275) |
| `parse_automatic_scaling_from_parsed` | [1303](../../../src/configuration.cpp#L1303) |
| `reject_unknown_qmom_option_keys` | [1322](../../../src/configuration.cpp#L1322) |
| `parse_realizability_check_mode_value` | [1340](../../../src/configuration.cpp#L1340) |
| `parse_realizability_check_mode` | [1359](../../../src/configuration.cpp#L1359) |
| `parse_moment_support` | [1370](../../../src/configuration.cpp#L1370) |
| `parse_phase_space_evolution_scheme_value` | [1393](../../../src/configuration.cpp#L1393) |
| `parse_node_advection_scheme_value` | [1412](../../../src/configuration.cpp#L1412) |
| `parse_node_advection_options_from_parsed` | [1431](../../../src/configuration.cpp#L1431) |
| `parse_realizability_policy_from_parsed` | [1479](../../../src/configuration.cpp#L1479) |
| `parse_coordinate_realizability_policies_from_parsed` | [1560](../../../src/configuration.cpp#L1560) |
| `parse_cqmom_diagnostic_level_value` | [1580](../../../src/configuration.cpp#L1580) |
| `parse_cqmom_diagnostic_level` | [1598](../../../src/configuration.cpp#L1598) |
| `parse_cqmom_moment_storage_value` | [1608](../../../src/configuration.cpp#L1608) |
| `parse_cqmom_moment_storage` | [1622](../../../src/configuration.cpp#L1622) |
| `parse_vandermonde_solver_value` | [1632](../../../src/configuration.cpp#L1632) |
| `parse_vandermonde_solver` | [1659](../../../src/configuration.cpp#L1659) |
| `has_output_section` | [1669](../../../src/configuration.cpp#L1669) |
| `is_output_key` | [1679](../../../src/configuration.cpp#L1679) |
| `validate_optional_phase_space_evolution_scheme` | [1684](../../../src/configuration.cpp#L1684) |
| `validate_optional_node_advection_options` | [1693](../../../src/configuration.cpp#L1693) |
| `validate_optional_realizability_policy` | [1698](../../../src/configuration.cpp#L1698) |
| `validate_cqmom_coordinate_realizability_ranks` | [1710](../../../src/configuration.cpp#L1710) |
| `reject_unknown_phase_space_evolution_sections` | [1724](../../../src/configuration.cpp#L1724) |
| `reject_unknown_phase_space_evolution_keys` | [1737](../../../src/configuration.cpp#L1737) |
| `parse_phase_space_evolution_scheme_from_parsed` | [1756](../../../src/configuration.cpp#L1756) |
| `validate_optional_cqmom_diagnostic_level` | [1776](../../../src/configuration.cpp#L1776) |
| `validate_optional_cqmom_moment_storage` | [1781](../../../src/configuration.cpp#L1781) |
| `validate_optional_vandermonde_solver` | [1786](../../../src/configuration.cpp#L1786) |
| `parse_coordinate_support` | [1791](../../../src/configuration.cpp#L1791) |
| `parse_cqmom_case_settings` | [1853](../../../src/configuration.cpp#L1853) |
| `initialize_moments` | [1974](../../../src/configuration.cpp#L1974) |
| `initialize_cqmom_dense_moments` | [2009](../../../src/configuration.cpp#L2009) |
| `initialize_cqmom_moment_data` | [2039](../../../src/configuration.cpp#L2039) |
| `parse_qmom_initial_condition` | [2086](../../../src/configuration.cpp#L2086) |
| `parse_cqmom_initial_condition` | [2182](../../../src/configuration.cpp#L2182) |
| `create_qmom_case` | [2295](../../../src/configuration.cpp#L2295) |
| `create_cqmom_dense_case` | [2316](../../../src/configuration.cpp#L2316) |
| `create_cqmom_lean_case` | [2373](../../../src/configuration.cpp#L2373) |
| `with_initial_moments` | [2392](../../../src/configuration.cpp#L2392) |
| `parse_qmom` | [2448](../../../src/configuration.cpp#L2448) |
| `parse_cqmom_dense` | [2503](../../../src/configuration.cpp#L2503) |
| `parse_cqmom_lean` | [2532](../../../src/configuration.cpp#L2532) |
| `parse_cqmom` | [2561](../../../src/configuration.cpp#L2561) |
| `parse_realizability_policy` | [2601](../../../src/configuration.cpp#L2601) |
| `parse_coordinate_realizability_policies` | [2630](../../../src/configuration.cpp#L2630) |
| `parse_phase_space_evolution_scheme` | [2650](../../../src/configuration.cpp#L2650) |
| `parse_node_advection_options` | [2669](../../../src/configuration.cpp#L2669) |
| `parse_qmom_phase_space_evolution_options` | [2688](../../../src/configuration.cpp#L2688) |
| `parse_cqmom_phase_space_evolution_options` | [2718](../../../src/configuration.cpp#L2718) |
| `parse_output_options` | [2744](../../../src/configuration.cpp#L2744) |
| `load_qmom` | [2820](../../../src/configuration.cpp#L2820) |
| `load_cqmom_dense` | [2830](../../../src/configuration.cpp#L2830) |
| `load_cqmom_lean` | [2840](../../../src/configuration.cpp#L2840) |
| `load_cqmom` | [2850](../../../src/configuration.cpp#L2850) |
| `load_output_options` | [2860](../../../src/configuration.cpp#L2860) |
| `load_node_advection_options` | [2870](../../../src/configuration.cpp#L2870) |
| `load_realizability_policy` | [2880](../../../src/configuration.cpp#L2880) |
| `load_qmom_phase_space_evolution_options` | [2890](../../../src/configuration.cpp#L2890) |
| `load_cqmom_phase_space_evolution_options` | [2900](../../../src/configuration.cpp#L2900) |

### `src/distributions.cpp`

[Source](../../../src/distributions.cpp) · [Responsibilities, contracts, and evidence](initialization_and_configuration.md)

| Function or overload family | Definition lines |
|---|---|
| `is_positive_finite` | [17](../../../src/distributions.cpp#L17) |
| `make_failure` | [22](../../../src/distributions.cpp#L22) |
| `invalid_input` | [33](../../../src/distributions.cpp#L33) |
| `make_tensor_failure` | [38](../../../src/distributions.cpp#L38) |
| `make_selected_failure` | [49](../../../src/distributions.cpp#L49) |
| `make_selected_component_failure` | [61](../../../src/distributions.cpp#L61) |
| `invalid_tensor_input` | [71](../../../src/distributions.cpp#L71) |
| `failure_at_order` | [76](../../../src/distributions.cpp#L76) |
| `selected_failure_at_order` | [81](../../../src/distributions.cpp#L81) |
| `finite_failure_status` | [86](../../../src/distributions.cpp#L86) |
| `check_common_inputs` | [95](../../../src/distributions.cpp#L95) |
| `check_maximum_order` | [108](../../../src/distributions.cpp#L108) |
| `finish_moments` | [117](../../../src/distributions.cpp#L117) |
| `MomentValueResult::ok` | [129](../../../src/distributions.cpp#L129) |
| `moment_value_failure` | [135](../../../src/distributions.cpp#L135) |
| `normalize_marginal_moment` | [140](../../../src/distributions.cpp#L140) |
| `evaluate_log_expectation` | [160](../../../src/distributions.cpp#L160) |
| `evaluate_log_normal_moment` | [192](../../../src/distributions.cpp#L192) |
| `evaluate_beta_moment` | [199](../../../src/distributions.cpp#L199) |
| `evaluate_gamma_moment` | [207](../../../src/distributions.cpp#L207) |
| `evaluate_chi_moment` | [214](../../../src/distributions.cpp#L214) |
| `evaluate_chi_squared_moment` | [221](../../../src/distributions.cpp#L221) |
| `evaluate_weibull_moment` | [228](../../../src/distributions.cpp#L228) |
| `compute_from_log_moments` | [236](../../../src/distributions.cpp#L236) |
| `check_selected_common_inputs` | [257](../../../src/distributions.cpp#L257) |
| `finish_selected_orders` | [289](../../../src/distributions.cpp#L289) |
| `zero_selected_moments` | [294](../../../src/distributions.cpp#L294) |
| `selected_failure_from_moment` | [299](../../../src/distributions.cpp#L299) |
| `NormalExpectationsResult::ok` | [311](../../../src/distributions.cpp#L311) |
| `normal_expectations` | [317](../../../src/distributions.cpp#L317) |
| `compute_normal_selected` | [343](../../../src/distributions.cpp#L343) |
| `compute_log_normal_selected` | [377](../../../src/distributions.cpp#L377) |
| `compute_beta_selected` | [410](../../../src/distributions.cpp#L410) |
| `compute_gamma_selected` | [437](../../../src/distributions.cpp#L437) |
| `compute_chi_selected` | [465](../../../src/distributions.cpp#L465) |
| `compute_chi_squared_selected` | [493](../../../src/distributions.cpp#L493) |
| `compute_weibull_selected` | [521](../../../src/distributions.cpp#L521) |
| `compute_normal` | [548](../../../src/distributions.cpp#L548) |
| `compute_log_normal` | [577](../../../src/distributions.cpp#L577) |
| `compute_beta` | [595](../../../src/distributions.cpp#L595) |
| `compute_gamma` | [607](../../../src/distributions.cpp#L607) |
| `compute_chi` | [620](../../../src/distributions.cpp#L620) |
| `compute_chi_squared` | [633](../../../src/distributions.cpp#L633) |
| `compute_weibull` | [647](../../../src/distributions.cpp#L647) |
| `RawMomentVisitor::operator()` | [664](../../../src/distributions.cpp#L664), [669](../../../src/distributions.cpp#L669), [674](../../../src/distributions.cpp#L674), [679](../../../src/distributions.cpp#L679), [684](../../../src/distributions.cpp#L684), [689](../../../src/distributions.cpp#L689), [694](../../../src/distributions.cpp#L694) |
| `SelectedRawMomentVisitor::operator()` | [707](../../../src/distributions.cpp#L707), [712](../../../src/distributions.cpp#L712), [717](../../../src/distributions.cpp#L717), [722](../../../src/distributions.cpp#L722), [727](../../../src/distributions.cpp#L727), [732](../../../src/distributions.cpp#L732), [737](../../../src/distributions.cpp#L737) |
| `finite_failure_status_from_value` | [743](../../../src/distributions.cpp#L743) |
| `add_compensated` | [752](../../../src/distributions.cpp#L752) |
| `add_selected_compensated` | [783](../../../src/distributions.cpp#L783) |
| `contains_order` | [817](../../../src/distributions.cpp#L817) |
| `order_position` | [829](../../../src/distributions.cpp#L829) |
| `compute_independent_product_marginal_moments` | [841](../../../src/distributions.cpp#L841), [851](../../../src/distributions.cpp#L851) |
| `first_requested_index_with_order` | [861](../../../src/distributions.cpp#L861) |
| `make_selected_raw_moment_result` | [875](../../../src/distributions.cpp#L875) |
| `compute_raw_moments` | [919](../../../src/distributions.cpp#L919), [935](../../../src/distributions.cpp#L935), [945](../../../src/distributions.cpp#L945), [1015](../../../src/distributions.cpp#L1015), [1096](../../../src/distributions.cpp#L1096), [1203](../../../src/distributions.cpp#L1203) |

### `src/internal/bjorck_pereyra.cpp`

[Source](../../../src/internal/bjorck_pereyra.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_all_finite` | [16](../../../src/internal/bjorck_pereyra.cpp#L16) |
| `failure_result` | [31](../../../src/internal/bjorck_pereyra.cpp#L31) |
| `node_scale` | [36](../../../src/internal/bjorck_pereyra.cpp#L36) |
| `denominator_is_singular` | [41](../../../src/internal/bjorck_pereyra.cpp#L41) |
| `solve_bjorck_pereyra` | [48](../../../src/internal/bjorck_pereyra.cpp#L48) |

### `src/internal/cqmom_active_counts.hpp`

[Source](../../../src/internal/cqmom_active_counts.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `parent_path_at_level` | [14](../../../src/internal/cqmom_active_counts.hpp#L14) |
| `find_active_branch_accumulator` | [33](../../../src/internal/cqmom_active_counts.hpp#L33) |
| `contains_node_index` | [49](../../../src/internal/cqmom_active_counts.hpp#L49) |
| `validate_contiguous_child_indices` | [54](../../../src/internal/cqmom_active_counts.hpp#L54) |
| `active_branch_counts_from_quadrature` | [67](../../../src/internal/cqmom_active_counts.hpp#L67) |

### `src/internal/cqmom_validation.hpp`

[Source](../../../src/internal/cqmom_validation.hpp) · [Responsibilities, contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `raw_moment_tensor_value_count` | [17](../../../src/internal/cqmom_validation.hpp#L17) |
| `same_cqmom_layout` | [53](../../../src/internal/cqmom_validation.hpp#L53) |
| `multi_index_within_extents` | [59](../../../src/internal/cqmom_validation.hpp#L59) |
| `validate_cqmom_moment_layout_matches_expected` | [75](../../../src/internal/cqmom_validation.hpp#L75) |
| `validate_fixed_cqmom_moment_layout` | [86](../../../src/internal/cqmom_validation.hpp#L86) |
| `validate_cqmom_moment_values_match_layout` | [92](../../../src/internal/cqmom_validation.hpp#L92) |
| `validate_fixed_cqmom_moment_data` | [104](../../../src/internal/cqmom_validation.hpp#L104) |
| `validate_cqmom_moment_data_matches_hierarchy` | [110](../../../src/internal/cqmom_validation.hpp#L110) |
| `validate_dense_cqmom_moment_tensor` | [118](../../../src/internal/cqmom_validation.hpp#L118) |
| `validate_dense_cqmom_reconstruction_tensor` | [135](../../../src/internal/cqmom_validation.hpp#L135) |
| `validate_cqmom_quadrature_storage` | [149](../../../src/internal/cqmom_validation.hpp#L149) |
| `validate_cqmom_moment_rate` | [179](../../../src/internal/cqmom_validation.hpp#L179) |

### `src/internal/embedded_runge_kutta.cpp`

[Source](../../../src/internal/embedded_runge_kutta.cpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `require_same_size` | [46](../../../src/internal/embedded_runge_kutta.cpp#L46) |
| `rhs_checked` | [54](../../../src/internal/embedded_runge_kutta.cpp#L54) |
| `embedded_rk_tableau` | [67](../../../src/internal/embedded_runge_kutta.cpp#L67) |
| `take_embedded_rk_step` | [79](../../../src/internal/embedded_runge_kutta.cpp#L79) |
| `embedded_rk_controller_exponent` | [130](../../../src/internal/embedded_runge_kutta.cpp#L130) |

### `src/internal/householder_qr.cpp`

[Source](../../../src/internal/householder_qr.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `finite_matrix` | [13](../../../src/internal/householder_qr.cpp#L13) |
| `stable_column_norm` | [28](../../../src/internal/householder_qr.cpp#L28) |
| `swap_columns` | [38](../../../src/internal/householder_qr.cpp#L38) |
| `qr_failure_result` | [50](../../../src/internal/householder_qr.cpp#L50) |
| `solve_failure_result` | [62](../../../src/internal/householder_qr.cpp#L62) |
| `numerical_rank_from_factor` | [68](../../../src/internal/householder_qr.cpp#L68) |
| `valid_permutation` | [84](../../../src/internal/householder_qr.cpp#L84) |
| `valid_qr_status` | [102](../../../src/internal/householder_qr.cpp#L102) |
| `finite_values` | [114](../../../src/internal/householder_qr.cpp#L114) |
| `rank_metadata_is_consistent` | [119](../../../src/internal/householder_qr.cpp#L119) |
| `diagnostic_rank_from_metadata` | [139](../../../src/internal/householder_qr.cpp#L139) |
| `factorization_metadata_is_consistent` | [150](../../../src/internal/householder_qr.cpp#L150) |
| `apply_householder_to_rhs` | [160](../../../src/internal/householder_qr.cpp#L160) |
| `factorize_qr` | [187](../../../src/internal/householder_qr.cpp#L187) |
| `solve_qr` | [281](../../../src/internal/householder_qr.cpp#L281), [342](../../../src/internal/householder_qr.cpp#L342) |

### `src/internal/jacobi_svd.cpp`

[Source](../../../src/internal/jacobi_svd.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `finite_matrix` | [16](../../../src/internal/jacobi_svd.cpp#L16) |
| `finite_values` | [31](../../../src/internal/jacobi_svd.cpp#L31) |
| `max_abs_entry` | [36](../../../src/internal/jacobi_svd.cpp#L36) |
| `scale_matrix` | [49](../../../src/internal/jacobi_svd.cpp#L49) |
| `stable_column_norm` | [60](../../../src/internal/jacobi_svd.cpp#L60) |
| `stable_column_dot` | [70](../../../src/internal/jacobi_svd.cpp#L70) |
| `rotate_columns` | [80](../../../src/internal/jacobi_svd.cpp#L80) |
| `copy_column` | [91](../../../src/internal/jacobi_svd.cpp#L91) |
| `resolved_jacobi_tolerance` | [99](../../../src/internal/jacobi_svd.cpp#L99) |
| `resolved_max_sweeps` | [109](../../../src/internal/jacobi_svd.cpp#L109) |
| `resolved_scale_exponent` | [118](../../../src/internal/jacobi_svd.cpp#L118) |
| `valid_options` | [128](../../../src/internal/jacobi_svd.cpp#L128) |
| `factorization_failure_result` | [134](../../../src/internal/jacobi_svd.cpp#L134) |
| `solve_failure_result` | [149](../../../src/internal/jacobi_svd.cpp#L149) |
| `numerical_rank_from_singular_values` | [155](../../../src/internal/jacobi_svd.cpp#L155) |
| `singular_values_are_nonnegative_descending` | [171](../../../src/internal/jacobi_svd.cpp#L171) |
| `valid_svd_status` | [187](../../../src/internal/jacobi_svd.cpp#L187) |
| `rank_metadata_is_consistent` | [200](../../../src/internal/jacobi_svd.cpp#L200) |
| `diagnostic_rank_from_metadata` | [221](../../../src/internal/jacobi_svd.cpp#L221) |
| `full_rank_singular_values_are_safe` | [232](../../../src/internal/jacobi_svd.cpp#L232) |
| `factorization_metadata_is_consistent` | [242](../../../src/internal/jacobi_svd.cpp#L242) |
| `pair_is_converged` | [256](../../../src/internal/jacobi_svd.cpp#L256) |
| `jacobi_tangent` | [265](../../../src/internal/jacobi_svd.cpp#L265) |
| `rotate_pair` | [282](../../../src/internal/jacobi_svd.cpp#L282) |
| `sort_singular_triplets` | [301](../../../src/internal/jacobi_svd.cpp#L301) |
| `completed_factorization` | [324](../../../src/internal/jacobi_svd.cpp#L324) |
| `factorize_svd` | [380](../../../src/internal/jacobi_svd.cpp#L380) |
| `solve_svd` | [444](../../../src/internal/jacobi_svd.cpp#L444), [508](../../../src/internal/jacobi_svd.cpp#L508) |
| `solve_svd_pseudoinverse` | [524](../../../src/internal/jacobi_svd.cpp#L524) |

### `src/internal/linear_algebra.cpp`

[Source](../../../src/internal/linear_algebra.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_infinity_norm` | [14](../../../src/internal/linear_algebra.cpp#L14) |
| `off_diagonal_norm` | [29](../../../src/internal/linear_algebra.cpp#L29) |
| `solve_dense` | [45](../../../src/internal/linear_algebra.cpp#L45) |
| `make_identity` | [116](../../../src/internal/linear_algebra.cpp#L116) |
| `residual_norm` | [126](../../../src/internal/linear_algebra.cpp#L126) |
| `condition_estimate` | [146](../../../src/internal/linear_algebra.cpp#L146) |
| `rotate` | [160](../../../src/internal/linear_algebra.cpp#L160) |
| `symmetric_jacobi_eigen` | [205](../../../src/internal/linear_algebra.cpp#L205) |

### `src/internal/linear_algebra.hpp`

[Source](../../../src/internal/linear_algebra.hpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `DenseMatrix::DenseMatrix` | [14](../../../src/internal/linear_algebra.hpp#L14) |
| `DenseMatrix::rows` | [16](../../../src/internal/linear_algebra.hpp#L16) |
| `DenseMatrix::cols` | [21](../../../src/internal/linear_algebra.hpp#L21) |
| `DenseMatrix::operator()` | [26](../../../src/internal/linear_algebra.hpp#L26), [35](../../../src/internal/linear_algebra.hpp#L35) |

### `src/internal/numeric_utils.hpp`

[Source](../../../src/internal/numeric_utils.hpp) · [Responsibilities, contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `all_finite` | [16](../../../src/internal/numeric_utils.hpp#L16) |
| `format_diagnostic_vector` | [21](../../../src/internal/numeric_utils.hpp#L21) |
| `phase_space_velocity_wrong_dimension_message` | [37](../../../src/internal/numeric_utils.hpp#L37) |
| `phase_space_velocity_non_finite_message` | [46](../../../src/internal/numeric_utils.hpp#L46) |
| `validate_phase_space_velocity_result` | [59](../../../src/internal/numeric_utils.hpp#L59) |

### `src/internal/performance.hpp`

[Source](../../../src/internal/performance.hpp) · [Responsibilities, contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `CounterSession::snapshot` | [81](../../../src/internal/performance.hpp#L81), [370](../../../src/internal/performance.hpp#L370) |
| `CounterSession::reset` | [86](../../../src/internal/performance.hpp#L86), [375](../../../src/internal/performance.hpp#L375) |
| `active_counters` | [100](../../../src/internal/performance.hpp#L100) |
| `CounterScope::CounterScope` | [108](../../../src/internal/performance.hpp#L108), [381](../../../src/internal/performance.hpp#L381) |
| `CounterScope::~CounterScope` | [116](../../../src/internal/performance.hpp#L116) |
| `enabled` | [125](../../../src/internal/performance.hpp#L125), [384](../../../src/internal/performance.hpp#L384) |
| `record_qmom_reconstruction` | [130](../../../src/internal/performance.hpp#L130), [389](../../../src/internal/performance.hpp#L389) |
| `record_cqmom_reconstruction` | [138](../../../src/internal/performance.hpp#L138), [390](../../../src/internal/performance.hpp#L390) |
| `record_cqmom_branch_solve` | [146](../../../src/internal/performance.hpp#L146), [391](../../../src/internal/performance.hpp#L391) |
| `record_wheeler_call` | [154](../../../src/internal/performance.hpp#L154), [392](../../../src/internal/performance.hpp#L392) |
| `record_wheeler_adaptive_reduction` | [162](../../../src/internal/performance.hpp#L162), [393](../../../src/internal/performance.hpp#L393) |
| `record_terminal_node_assembled` | [170](../../../src/internal/performance.hpp#L170), [394](../../../src/internal/performance.hpp#L394) |
| `record_moment_layout_build` | [178](../../../src/internal/performance.hpp#L178), [395](../../../src/internal/performance.hpp#L395) |
| `record_cqmom_prepared_metadata_build` | [186](../../../src/internal/performance.hpp#L186), [396](../../../src/internal/performance.hpp#L396) |
| `record_cqmom_prepared_metadata_use` | [194](../../../src/internal/performance.hpp#L194), [397](../../../src/internal/performance.hpp#L397) |
| `record_cqmom_dynamic_metadata_build` | [202](../../../src/internal/performance.hpp#L202), [398](../../../src/internal/performance.hpp#L398) |
| `record_vandermonde_solve` | [210](../../../src/internal/performance.hpp#L210), [399](../../../src/internal/performance.hpp#L399) |
| `record_bjorck_pereyra_transform_update` | [230](../../../src/internal/performance.hpp#L230), [400](../../../src/internal/performance.hpp#L400) |
| `record_dense_linear_solve` | [238](../../../src/internal/performance.hpp#L238), [401](../../../src/internal/performance.hpp#L401) |
| `record_residual_check` | [246](../../../src/internal/performance.hpp#L246), [402](../../../src/internal/performance.hpp#L402) |
| `record_condition_estimate` | [254](../../../src/internal/performance.hpp#L254), [403](../../../src/internal/performance.hpp#L403) |
| `record_phase_space_evolution_macro_step` | [262](../../../src/internal/performance.hpp#L262), [404](../../../src/internal/performance.hpp#L404) |
| `record_advection_operator_call` | [270](../../../src/internal/performance.hpp#L270), [405](../../../src/internal/performance.hpp#L405) |
| `record_moment_source_operator_call` | [278](../../../src/internal/performance.hpp#L278), [406](../../../src/internal/performance.hpp#L406) |
| `record_moment_reconstruction_call` | [286](../../../src/internal/performance.hpp#L286), [407](../../../src/internal/performance.hpp#L407) |
| `record_post_source_inversion` | [294](../../../src/internal/performance.hpp#L294), [408](../../../src/internal/performance.hpp#L408) |
| `record_active_count_derivation` | [302](../../../src/internal/performance.hpp#L302), [409](../../../src/internal/performance.hpp#L409) |
| `record_active_count_validation` | [310](../../../src/internal/performance.hpp#L310), [410](../../../src/internal/performance.hpp#L410) |
| `record_active_count_transition_build` | [318](../../../src/internal/performance.hpp#L318), [411](../../../src/internal/performance.hpp#L411) |
| `record_source_rk_stage` | [326](../../../src/internal/performance.hpp#L326), [412](../../../src/internal/performance.hpp#L412) |
| `record_source_retry` | [334](../../../src/internal/performance.hpp#L334), [413](../../../src/internal/performance.hpp#L413) |
| `record_node_rhs_evaluation` | [342](../../../src/internal/performance.hpp#L342), [414](../../../src/internal/performance.hpp#L414) |
| `record_embedded_rk_step` | [350](../../../src/internal/performance.hpp#L350), [415](../../../src/internal/performance.hpp#L415) |

### `src/internal/qmom_state_access.hpp`

[Source](../../../src/internal/qmom_state_access.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomStateAccess::from_quadrature` | [20](../../../src/internal/qmom_state_access.hpp#L20) |

### `src/internal/recursive_quadrature.hpp`

[Source](../../../src/internal/recursive_quadrature.hpp) · [Responsibilities, contracts, and evidence](cqmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `RecursiveQuadratureResult::ok` | [300](../../../src/internal/recursive_quadrature.hpp#L300) |

### `src/internal/vandermonde.cpp`

[Source](../../../src/internal/vandermonde.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_all_finite` | [19](../../../src/internal/vandermonde.cpp#L19) |
| `failure_result` | [34](../../../src/internal/vandermonde.cpp#L34) |
| `select_automatic_solver` | [39](../../../src/internal/vandermonde.cpp#L39) |
| `resolve_solver` | [44](../../../src/internal/vandermonde.cpp#L44) |
| `vandermonde_result` | [61](../../../src/internal/vandermonde.cpp#L61) |
| `qr_linear_solve_result` | [68](../../../src/internal/vandermonde.cpp#L68) |
| `svd_linear_solve_result` | [85](../../../src/internal/vandermonde.cpp#L85) |
| `residual_is_accepted` | [102](../../../src/internal/vandermonde.cpp#L102) |
| `make_conditional_vandermonde_matrix` | [115](../../../src/internal/vandermonde.cpp#L115) |
| `solve_conditional_vandermonde_system` | [129](../../../src/internal/vandermonde.cpp#L129) |

### `src/moment_reconstruction.cpp`

[Source](../../../src/moment_reconstruction.cpp) · [Responsibilities, contracts, and evidence](moment_evaluation.md)

| Function or overload family | Definition lines |
|---|---|
| `integer_power` | [20](../../../src/moment_reconstruction.cpp#L20) |
| `validate_reconstructed_value` | [40](../../../src/moment_reconstruction.cpp#L40) |
| `is_integer_exponent` | [56](../../../src/moment_reconstruction.cpp#L56) |
| `validate_generalized_value` | [61](../../../src/moment_reconstruction.cpp#L61) |
| `validate_generalized_exponent` | [69](../../../src/moment_reconstruction.cpp#L69) |
| `generalized_power` | [77](../../../src/moment_reconstruction.cpp#L77) |
| `cqmom_reconstruction_quadrature_messages` | [89](../../../src/moment_reconstruction.cpp#L89) |
| `reconstruct_cqmom_moment_value` | [149](../../../src/moment_reconstruction.cpp#L149) |
| `reconstruct_cqmom_moments_impl` | [197](../../../src/moment_reconstruction.cpp#L197) |
| `reconstruct_raw_moments` | [236](../../../src/moment_reconstruction.cpp#L236), [296](../../../src/moment_reconstruction.cpp#L296), [313](../../../src/moment_reconstruction.cpp#L313) |
| `evaluate_generalized_moment` | [269](../../../src/moment_reconstruction.cpp#L269), [332](../../../src/moment_reconstruction.cpp#L332) |
| `reconstruct_cqmom_moments` | [359](../../../src/moment_reconstruction.cpp#L359) |
| `reconstruct_trusted_cqmom_moments` | [367](../../../src/moment_reconstruction.cpp#L367) |

### `src/moment_source.cpp`

[Source](../../../src/moment_source.cpp) · [Responsibilities, contracts, and evidence](moment_sources.md)

| Function or overload family | Definition lines |
|---|---|
| `status_from_ssp_rk3` | [23](../../../src/moment_source.cpp#L23) |
| `message_from_status` | [37](../../../src/moment_source.cpp#L37) |
| `diagnostics_from_ssp_rk3` | [55](../../../src/moment_source.cpp#L55) |
| `invalid_result` | [66](../../../src/moment_source.cpp#L66) |
| `invalid_moment_data_result` | [73](../../../src/moment_source.cpp#L73) |
| `failed_moment_data_result` | [80](../../../src/moment_source.cpp#L80) |
| `validate_source_context` | [87](../../../src/moment_source.cpp#L87) |
| `validate_fixed_cqmom_moment_data` | [102](../../../src/moment_source.cpp#L102) |
| `validate_rate` | [109](../../../src/moment_source.cpp#L109) |
| `validate_source_rate` | [116](../../../src/moment_source.cpp#L116) |
| `validate_boundary_rate` | [128](../../../src/moment_source.cpp#L128) |
| `validate_cqmom_boundary_rate` | [140](../../../src/moment_source.cpp#L140) |
| `validate_cqmom_stage_state` | [147](../../../src/moment_source.cpp#L147) |
| `validate_coordinate_realizability_policies` | [159](../../../src/moment_source.cpp#L159) |
| `has_coordinate_realizability_policy` | [171](../../../src/moment_source.cpp#L171) |
| `extract_coordinate_marginal` | [178](../../../src/moment_source.cpp#L178) |
| `accepts_coordinate_realizability_policies` | [210](../../../src/moment_source.cpp#L210) |
| `boundary_callback_error` | [231](../../../src/moment_source.cpp#L231) |
| `cqmom_boundary_callback_error` | [236](../../../src/moment_source.cpp#L236) |
| `validate_request` | [241](../../../src/moment_source.cpp#L241), [267](../../../src/moment_source.cpp#L267) |
| `validate_cqmom_source_update_controls` | [277](../../../src/moment_source.cpp#L277) |
| `validate_moment_data_request` | [299](../../../src/moment_source.cpp#L299), [316](../../../src/moment_source.cpp#L316) |
| `validate_trusted_moment_data_request` | [307](../../../src/moment_source.cpp#L307), [326](../../../src/moment_source.cpp#L326) |
| `compose_source_and_boundary` | [336](../../../src/moment_source.cpp#L336) |
| `advance_moments_with_source` | [362](../../../src/moment_source.cpp#L362), [759](../../../src/moment_source.cpp#L759), [779](../../../src/moment_source.cpp#L779) |
| `compose_cqmom_source_and_boundary` | [406](../../../src/moment_source.cpp#L406) |
| `advance_validated_cqmom_moments_with_source` | [433](../../../src/moment_source.cpp#L433) |
| `QmomSourceContextAccess::create` | [497](../../../src/moment_source.cpp#L497) |
| `CqmomSourceContextAccess::create` | [503](../../../src/moment_source.cpp#L503) |
| `QmomSourceContext::QmomSourceContext` | [513](../../../src/moment_source.cpp#L513), [521](../../../src/moment_source.cpp#L521) |
| `QmomSourceContext::time` | [534](../../../src/moment_source.cpp#L534) |
| `QmomSourceContext::moment_count` | [539](../../../src/moment_source.cpp#L539) |
| `QmomSourceContext::size` | [544](../../../src/moment_source.cpp#L544) |
| `QmomSourceContext::values` | [549](../../../src/moment_source.cpp#L549) |
| `QmomSourceContext::value` | [554](../../../src/moment_source.cpp#L554) |
| `QmomSourceContext::moment` | [559](../../../src/moment_source.cpp#L559) |
| `CqmomSourceContext::CqmomSourceContext` | [586](../../../src/moment_source.cpp#L586), [597](../../../src/moment_source.cpp#L597) |
| `CqmomSourceContext::time` | [611](../../../src/moment_source.cpp#L611) |
| `CqmomSourceContext::moment_count` | [616](../../../src/moment_source.cpp#L616) |
| `CqmomSourceContext::dimension_count` | [621](../../../src/moment_source.cpp#L621) |
| `CqmomSourceContext::size` | [626](../../../src/moment_source.cpp#L626) |
| `CqmomSourceContext::indices` | [631](../../../src/moment_source.cpp#L631) |
| `CqmomSourceContext::values` | [636](../../../src/moment_source.cpp#L636) |
| `CqmomSourceContext::index` | [641](../../../src/moment_source.cpp#L641) |
| `CqmomSourceContext::value` | [646](../../../src/moment_source.cpp#L646), [661](../../../src/moment_source.cpp#L661) |
| `CqmomSourceContext::position` | [651](../../../src/moment_source.cpp#L651) |
| `CqmomSourceContext::moment` | [671](../../../src/moment_source.cpp#L671) |
| `CqmomSourceContext::generalized_moment` | [686](../../../src/moment_source.cpp#L686) |
| `CqmomSourceContext::for_each_terminal_node` | [706](../../../src/moment_source.cpp#L706) |
| `CqmomSourceContext::reconstructed_stage_quadrature` | [727](../../../src/moment_source.cpp#L727) |
| `CqmomSourceContext::hierarchy` | [754](../../../src/moment_source.cpp#L754) |
| `advance_cqmom_moments_with_source` | [795](../../../src/moment_source.cpp#L795), [810](../../../src/moment_source.cpp#L810) |
| `advance_trusted_qmom_moments_with_source` | [829](../../../src/moment_source.cpp#L829), [854](../../../src/moment_source.cpp#L854) |
| `advance_trusted_cqmom_moments_with_source` | [881](../../../src/moment_source.cpp#L881), [896](../../../src/moment_source.cpp#L896), [915](../../../src/moment_source.cpp#L915), [933](../../../src/moment_source.cpp#L933) |

### `src/node_advection.cpp`

[Source](../../../src/node_advection.cpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `all_finite` | [21](../../../src/node_advection.cpp#L21) |
| `diagnostics_from_embedded_rk` | [26](../../../src/node_advection.cpp#L26) |
| `pack_nodes` | [40](../../../src/node_advection.cpp#L40) |
| `unpack_nodes` | [62](../../../src/node_advection.cpp#L62) |
| `validate_weights` | [76](../../../src/node_advection.cpp#L76) |
| `validate_abscissas` | [88](../../../src/node_advection.cpp#L88) |
| `node_rhs_from_velocity` | [111](../../../src/node_advection.cpp#L111) |
| `add_scaled` | [148](../../../src/node_advection.cpp#L148) |
| `integrate_nodes_embedded_rk` | [162](../../../src/node_advection.cpp#L162) |
| `integrate_nodes_ssp_rk3` | [198](../../../src/node_advection.cpp#L198) |
| `integrate_nodes` | [259](../../../src/node_advection.cpp#L259) |
| `scheme_name` | [290](../../../src/node_advection.cpp#L290) |
| `advance_nodes` | [306](../../../src/node_advection.cpp#L306), [334](../../../src/node_advection.cpp#L334), [349](../../../src/node_advection.cpp#L349) |
| `format_node_advection_options` | [389](../../../src/node_advection.cpp#L389) |

### `src/output.cpp`

[Source](../../../src/output.cpp) · [Responsibilities, contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `qmom_label` | [19](../../../src/output.cpp#L19) |
| `cqmom_label` | [24](../../../src/output.cpp#L24) |
| `sequential_indices` | [40](../../../src/output.cpp#L40) |
| `qmom_labels` | [56](../../../src/output.cpp#L56) |
| `cqmom_labels` | [72](../../../src/output.cpp#L72) |
| `flatten_selected_indices` | [88](../../../src/output.cpp#L88) |
| `all_cqmom_labels` | [99](../../../src/output.cpp#L99) |
| `cqmom_abscissa_label` | [104](../../../src/output.cpp#L104) |
| `qmom_quadrature_labels` | [109](../../../src/output.cpp#L109) |
| `cqmom_quadrature_labels` | [119](../../../src/output.cpp#L119) |
| `MomentWriter::for_qmom` | [143](../../../src/output.cpp#L143), [149](../../../src/output.cpp#L149) |
| `MomentWriter::for_cqmom` | [155](../../../src/output.cpp#L155), [162](../../../src/output.cpp#L162), [169](../../../src/output.cpp#L169), [174](../../../src/output.cpp#L174) |
| `MomentWriter::MomentWriter` | [179](../../../src/output.cpp#L179), [194](../../../src/output.cpp#L194), [209](../../../src/output.cpp#L209) |
| `MomentWriter::write_header` | [223](../../../src/output.cpp#L223) |
| `MomentWriter::write_row` | [237](../../../src/output.cpp#L237) |
| `MomentWriter::write_aligned_row` | [264](../../../src/output.cpp#L264) |
| `MomentWriter::write` | [291](../../../src/output.cpp#L291), [304](../../../src/output.cpp#L304), [321](../../../src/output.cpp#L321), [334](../../../src/output.cpp#L334) |
| `QuadratureWriter::for_qmom` | [348](../../../src/output.cpp#L348) |
| `QuadratureWriter::for_cqmom` | [353](../../../src/output.cpp#L353) |
| `QuadratureWriter::QuadratureWriter` | [358](../../../src/output.cpp#L358) |
| `QuadratureWriter::write_header` | [373](../../../src/output.cpp#L373) |
| `QuadratureWriter::write_row` | [387](../../../src/output.cpp#L387) |
| `QuadratureWriter::write_qmom_row` | [414](../../../src/output.cpp#L414) |
| `QuadratureWriter::write` | [444](../../../src/output.cpp#L444), [449](../../../src/output.cpp#L449), [454](../../../src/output.cpp#L454) |

### `src/phase_space_evolution.cpp`

[Source](../../../src/phase_space_evolution.cpp) · [Responsibilities, contracts, and evidence](phase_space_evolution.md)

| Function or overload family | Definition lines |
|---|---|
| `apply_batch` | [49](../../../src/phase_space_evolution.cpp#L49) |
| `validate_time_input` | [62](../../../src/phase_space_evolution.cpp#L62) |
| `validate_phase_space_evolution_scheme` | [74](../../../src/phase_space_evolution.cpp#L74) |
| `execute_phase_space_evolution_sequence` | [87](../../../src/phase_space_evolution.cpp#L87) |
| `validate_qmom_input` | [122](../../../src/phase_space_evolution.cpp#L122) |
| `validate_boundary` | [158](../../../src/phase_space_evolution.cpp#L158) |
| `validate_cqmom_boundary` | [166](../../../src/phase_space_evolution.cpp#L166) |
| `successful_reconstruction` | [174](../../../src/phase_space_evolution.cpp#L174) |
| `record_advection_diagnostics` | [179](../../../src/phase_space_evolution.cpp#L179) |
| `record_source_update_diagnostics` | [195](../../../src/phase_space_evolution.cpp#L195) |
| `qmom_from_quadrature` | [202](../../../src/phase_space_evolution.cpp#L202) |
| `add_advection_totals` | [220](../../../src/phase_space_evolution.cpp#L220) |
| `finalize_diagnostics` | [229](../../../src/phase_space_evolution.cpp#L229) |
| `record_failure_diagnostics` | [235](../../../src/phase_space_evolution.cpp#L235) |
| `qmom_failure` | [242](../../../src/phase_space_evolution.cpp#L242) |
| `select_post_source_qmom` | [250](../../../src/phase_space_evolution.cpp#L250) |
| `validate_qmom_advection_input` | [314](../../../src/phase_space_evolution.cpp#L314) |
| `validate_qmom_moment_source_input` | [327](../../../src/phase_space_evolution.cpp#L327) |
| `advance_qmom_advection_stage` | [358](../../../src/phase_space_evolution.cpp#L358) |
| `advance_qmom_moment_source_stage` | [408](../../../src/phase_space_evolution.cpp#L408) |
| `neutral_cqmom_message` | [471](../../../src/phase_space_evolution.cpp#L471) |
| `neutral_cqmom_diagnostics` | [484](../../../src/phase_space_evolution.cpp#L484) |
| `neutral_cqmom_result` | [493](../../../src/phase_space_evolution.cpp#L493) |
| `operator<` | [507](../../../src/phase_space_evolution.cpp#L507) |
| `active_branch_key` | [520](../../../src/phase_space_evolution.cpp#L520) |
| `path_has_prefix` | [525](../../../src/phase_space_evolution.cpp#L525) |
| `has_changed_strict_ancestor` | [541](../../../src/phase_space_evolution.cpp#L541) |
| `record_cqmom_node_count_transitions` | [559](../../../src/phase_space_evolution.cpp#L559) |
| `advance_qmom_phase_space_impl` | [649](../../../src/phase_space_evolution.cpp#L649) |
| `advance_qmom_advection` | [691](../../../src/phase_space_evolution.cpp#L691) |
| `advance_qmom_moment_source_impl` | [716](../../../src/phase_space_evolution.cpp#L716) |
| `advance_qmom_moment_source` | [744](../../../src/phase_space_evolution.cpp#L744), [751](../../../src/phase_space_evolution.cpp#L751) |
| `advance_qmom_phase_space` | [758](../../../src/phase_space_evolution.cpp#L758), [765](../../../src/phase_space_evolution.cpp#L765) |
| `initialize_cqmom_diagnostics` | [775](../../../src/phase_space_evolution.cpp#L775), [782](../../../src/phase_space_evolution.cpp#L782), [789](../../../src/phase_space_evolution.cpp#L789) |
| `validate_dense_cqmom_moments_for_phase_space_evolution` | [802](../../../src/phase_space_evolution.cpp#L802) |
| `validate_lean_cqmom_moments_for_phase_space_evolution` | [809](../../../src/phase_space_evolution.cpp#L809) |
| `validate_cqmom_state_for_phase_space_evolution` | [816](../../../src/phase_space_evolution.cpp#L816) |
| `validate_cqmom_advection_input` | [831](../../../src/phase_space_evolution.cpp#L831) |
| `validate_cqmom_moment_source_input` | [841](../../../src/phase_space_evolution.cpp#L841) |
| `validate_cqmom_phase_space_evolution_input` | [853](../../../src/phase_space_evolution.cpp#L853) |
| `CqmomStageOutcome::ok` | [880](../../../src/phase_space_evolution.cpp#L880) |
| `take_cqmom_advection_active_counts` | [894](../../../src/phase_space_evolution.cpp#L894) |
| `cqmom_state_failure` | [903](../../../src/phase_space_evolution.cpp#L903) |
| `cqmom_success` | [911](../../../src/phase_space_evolution.cpp#L911) |
| `materialize_cqmom_phase_space_evolution_result` | [916](../../../src/phase_space_evolution.cpp#L916) |
| `dense_state_from_accepted_quadrature` | [937](../../../src/phase_space_evolution.cpp#L937) |
| `lean_state_from_accepted_quadrature` | [953](../../../src/phase_space_evolution.cpp#L953) |
| `prepare_cqmom_phase_space_evolution_working_state` | [964](../../../src/phase_space_evolution.cpp#L964) |
| `advance_cqmom_advection_stage` | [1008](../../../src/phase_space_evolution.cpp#L1008) |
| `advance_cqmom_moment_source_stage` | [1105](../../../src/phase_space_evolution.cpp#L1105) |
| `advance_cqmom_advection_impl` | [1244](../../../src/phase_space_evolution.cpp#L1244) |
| `advance_cqmom_moment_source_impl` | [1251](../../../src/phase_space_evolution.cpp#L1251) |
| `advance_cqmom_phase_space_state_impl` | [1258](../../../src/phase_space_evolution.cpp#L1258) |
| `advance_cqmom_advection` | [1318](../../../src/phase_space_evolution.cpp#L1318) |
| `advance_cqmom_moment_source_public_impl` | [1340](../../../src/phase_space_evolution.cpp#L1340) |
| `advance_cqmom_moment_source` | [1363](../../../src/phase_space_evolution.cpp#L1363), [1369](../../../src/phase_space_evolution.cpp#L1369) |
| `advance_cqmom_phase_space` | [1376](../../../src/phase_space_evolution.cpp#L1376), [1383](../../../src/phase_space_evolution.cpp#L1383) |
| `advance_qmom_advection_batch` | [1391](../../../src/phase_space_evolution.cpp#L1391) |
| `advance_cqmom_advection_batch` | [1399](../../../src/phase_space_evolution.cpp#L1399) |
| `advance_qmom_moment_source_batch` | [1407](../../../src/phase_space_evolution.cpp#L1407), [1424](../../../src/phase_space_evolution.cpp#L1424) |
| `advance_cqmom_moment_source_batch` | [1416](../../../src/phase_space_evolution.cpp#L1416), [1434](../../../src/phase_space_evolution.cpp#L1434) |
| `advance_qmom_phase_space_batch` | [1443](../../../src/phase_space_evolution.cpp#L1443), [1461](../../../src/phase_space_evolution.cpp#L1461) |
| `advance_cqmom_phase_space_batch` | [1452](../../../src/phase_space_evolution.cpp#L1452), [1471](../../../src/phase_space_evolution.cpp#L1471) |
| `advance_qmom_advection_batch_with_index` | [1481](../../../src/phase_space_evolution.cpp#L1481) |
| `advance_cqmom_advection_batch_with_index` | [1497](../../../src/phase_space_evolution.cpp#L1497) |
| `advance_qmom_moment_source_batch_with_index` | [1513](../../../src/phase_space_evolution.cpp#L1513), [1544](../../../src/phase_space_evolution.cpp#L1544) |
| `advance_cqmom_moment_source_batch_with_index` | [1529](../../../src/phase_space_evolution.cpp#L1529), [1566](../../../src/phase_space_evolution.cpp#L1566) |
| `advance_qmom_phase_space_batch_with_index` | [1587](../../../src/phase_space_evolution.cpp#L1587), [1632](../../../src/phase_space_evolution.cpp#L1632) |
| `advance_cqmom_phase_space_batch_with_index` | [1609](../../../src/phase_space_evolution.cpp#L1609), [1657](../../../src/phase_space_evolution.cpp#L1657) |

### `src/public_api.cpp`

[Source](../../../src/public_api.cpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `to_wheeler_options` | [27](../../../src/public_api.cpp#L27) |
| `compute_qmom_unscaled` | [34](../../../src/public_api.cpp#L34) |
| `compute_qmom_with_automatic_scaling` | [51](../../../src/public_api.cpp#L51) |
| `to_qmom_status` | [73](../../../src/public_api.cpp#L73) |
| `hierarchy` | [95](../../../src/public_api.cpp#L95) |
| `qmom_options` | [102](../../../src/public_api.cpp#L102) |
| `conditional_options` | [109](../../../src/public_api.cpp#L109) |
| `adaptive_branch_counts` | [116](../../../src/public_api.cpp#L116) |
| `diagnostic_level` | [123](../../../src/public_api.cpp#L123) |
| `vandermonde_solver` | [130](../../../src/public_api.cpp#L130) |
| `coordinate_support` | [137](../../../src/public_api.cpp#L137) |
| `moment_storage` | [144](../../../src/public_api.cpp#L144) |
| `diagnostic_level_name` | [149](../../../src/public_api.cpp#L149) |
| `vandermonde_solver_name` | [163](../../../src/public_api.cpp#L163) |
| `compute_qmom` | [182](../../../src/public_api.cpp#L182), [187](../../../src/public_api.cpp#L187), [193](../../../src/public_api.cpp#L193) |
| `QmomState::QmomState` | [198](../../../src/public_api.cpp#L198) |
| `QmomCase::initial_state` | [203](../../../src/public_api.cpp#L203) |
| `to_string` | [227](../../../src/public_api.cpp#L227) |
| `format_cqmom_settings` | [239](../../../src/public_api.cpp#L239) |
| `CqmomCase::CqmomCase` | [251](../../../src/public_api.cpp#L251), [253](../../../src/public_api.cpp#L253), [264](../../../src/public_api.cpp#L264) |
| `CqmomCase::operator=` | [255](../../../src/public_api.cpp#L255), [266](../../../src/public_api.cpp#L266) |
| `CqmomCase::hierarchy` | [270](../../../src/public_api.cpp#L270) |
| `CqmomCase::qmom_options` | [275](../../../src/public_api.cpp#L275) |
| `CqmomCase::conditional_options` | [280](../../../src/public_api.cpp#L280) |
| `CqmomCase::adaptive_branch_counts` | [285](../../../src/public_api.cpp#L285) |
| `CqmomCase::diagnostic_level` | [290](../../../src/public_api.cpp#L290) |
| `CqmomCase::vandermonde_solver` | [295](../../../src/public_api.cpp#L295) |
| `CqmomCase::coordinate_support` | [300](../../../src/public_api.cpp#L300) |
| `CqmomCase::moment_storage` | [305](../../../src/public_api.cpp#L305) |
| `CqmomCase::initial_state` | [310](../../../src/public_api.cpp#L310) |
| `CqmomState::CqmomState` | [339](../../../src/public_api.cpp#L339), [341](../../../src/public_api.cpp#L341), [352](../../../src/public_api.cpp#L352) |
| `CqmomState::operator=` | [343](../../../src/public_api.cpp#L343), [354](../../../src/public_api.cpp#L354) |
| `CqmomState::hierarchy` | [358](../../../src/public_api.cpp#L358) |
| `CqmomState::qmom_options` | [363](../../../src/public_api.cpp#L363) |
| `CqmomState::conditional_options` | [368](../../../src/public_api.cpp#L368) |
| `CqmomState::adaptive_branch_counts` | [373](../../../src/public_api.cpp#L373) |
| `CqmomState::diagnostic_level` | [378](../../../src/public_api.cpp#L378) |
| `CqmomState::vandermonde_solver` | [383](../../../src/public_api.cpp#L383) |
| `CqmomState::coordinate_support` | [388](../../../src/public_api.cpp#L388) |
| `CqmomState::moment_storage` | [393](../../../src/public_api.cpp#L393) |
| `CqmomState::active_branch_counts` | [398](../../../src/public_api.cpp#L398) |
| `CqmomState::moment_indices` | [403](../../../src/public_api.cpp#L403) |
| `CqmomState::moment_values` | [410](../../../src/public_api.cpp#L410) |
| `CqmomState::moment` | [417](../../../src/public_api.cpp#L417) |
| `validate_prepared_fixed_metadata` | [449](../../../src/public_api.cpp#L449), [461](../../../src/public_api.cpp#L461) |
| `make_prepared_fixed_metadata` | [469](../../../src/public_api.cpp#L469), [488](../../../src/public_api.cpp#L488) |
| `CqmomConstructionAccess::create_case` | [508](../../../src/public_api.cpp#L508), [514](../../../src/public_api.cpp#L514) |
| `CqmomConstructionAccess::create_state` | [520](../../../src/public_api.cpp#L520), [525](../../../src/public_api.cpp#L525), [531](../../../src/public_api.cpp#L531), [540](../../../src/public_api.cpp#L540), [547](../../../src/public_api.cpp#L547), [552](../../../src/public_api.cpp#L552), [558](../../../src/public_api.cpp#L558), [567](../../../src/public_api.cpp#L567) |
| `CqmomConstructionAccess::case_dense_case` | [574](../../../src/public_api.cpp#L574) |
| `CqmomConstructionAccess::case_lean_case` | [579](../../../src/public_api.cpp#L579) |
| `CqmomConstructionAccess::state_dense_case` | [584](../../../src/public_api.cpp#L584) |
| `CqmomConstructionAccess::state_lean_case` | [589](../../../src/public_api.cpp#L589) |
| `CqmomConstructionAccess::state_active_branch_counts` | [594](../../../src/public_api.cpp#L594) |
| `CqmomConstructionAccess::state_accepted_quadrature` | [599](../../../src/public_api.cpp#L599) |
| `CqmomConstructionAccess::case_prepared_fixed_metadata` | [604](../../../src/public_api.cpp#L604) |
| `CqmomConstructionAccess::state_prepared_fixed_metadata` | [609](../../../src/public_api.cpp#L609) |
| `CqmomConstructionAccess::create_quadrature` | [614](../../../src/public_api.cpp#L614) |
| `CqmomQuadrature::CqmomQuadrature` | [624](../../../src/public_api.cpp#L624) |
| `CqmomQuadrature::hierarchy` | [632](../../../src/public_api.cpp#L632) |
| `CqmomQuadrature::number_of_dimensions` | [637](../../../src/public_api.cpp#L637) |
| `CqmomQuadrature::total_nodes` | [642](../../../src/public_api.cpp#L642) |
| `CqmomQuadrature::weights` | [647](../../../src/public_api.cpp#L647) |
| `CqmomQuadrature::abscissas` | [652](../../../src/public_api.cpp#L652) |
| `CqmomQuadrature::adaptive_branch_counts_enabled` | [657](../../../src/public_api.cpp#L657) |
| `CqmomQuadrature::has_rectangular_storage` | [662](../../../src/public_api.cpp#L662) |
| `CqmomQuadrature::terminal_paths` | [667](../../../src/public_api.cpp#L667) |
| `CqmomQuadrature::weight` | [672](../../../src/public_api.cpp#L672) |
| `CqmomQuadrature::abscissa` | [681](../../../src/public_api.cpp#L681) |
| `CqmomQuadrature::node` | [694](../../../src/public_api.cpp#L694) |
| `to_cqmom_status` | [707](../../../src/public_api.cpp#L707) |
| `to_cqmom_diagnostics` | [723](../../../src/public_api.cpp#L723), [765](../../../src/public_api.cpp#L765) |
| `phase_space_index_to_conditioning_path` | [777](../../../src/public_api.cpp#L777) |
| `terminal_paths_for_public_quadrature` | [787](../../../src/public_api.cpp#L787) |
| `failed_phase_space_dimension` | [803](../../../src/public_api.cpp#L803) |
| `to_cqmom_result` | [813](../../../src/public_api.cpp#L813) |
| `compute_cqmom` | [830](../../../src/public_api.cpp#L830), [835](../../../src/public_api.cpp#L835), [970](../../../src/public_api.cpp#L970) |
| `compute_cqmom_internal` | [879](../../../src/public_api.cpp#L879), [936](../../../src/public_api.cpp#L936) |

### `src/realizability.cpp`

[Source](../../../src/realizability.cpp) · [Responsibilities, contracts, and evidence](realizability.md)

| Function or overload family | Definition lines |
|---|---|
| `RealizabilityResult::realizable` | [19](../../../src/realizability.cpp#L19) |
| `RealizabilityResult::on_boundary` | [24](../../../src/realizability.cpp#L24) |
| `support_name` | [59](../../../src/realizability.cpp#L59) |
| `make_result` | [73](../../../src/realizability.cpp#L73) |
| `vector_scale` | [79](../../../src/realizability.cpp#L79) |
| `symmetric_jacobi_eigen` | [89](../../../src/realizability.cpp#L89) |
| `analyze_matrix` | [122](../../../src/realizability.cpp#L122) |
| `update_limiting_matrix` | [187](../../../src/realizability.cpp#L187) |
| `incorporate_matrix` | [201](../../../src/realizability.cpp#L201) |
| `hankel_matrix` | [228](../../../src/realizability.cpp#L228) |
| `interval_quadratic_localizing_matrix` | [243](../../../src/realizability.cpp#L243) |
| `interval_linear_localizing_matrix` | [260](../../../src/realizability.cpp#L260) |
| `extension_column` | [278](../../../src/realizability.cpp#L278) |
| `projection_on_eigenvector` | [288](../../../src/realizability.cpp#L288) |
| `check_range_consistency` | [298](../../../src/realizability.cpp#L298) |
| `incorporate_range_consistency` | [326](../../../src/realizability.cpp#L326) |
| `finish_classification` | [351](../../../src/realizability.cpp#L351) |
| `validate_common_inputs` | [367](../../../src/realizability.cpp#L367) |
| `check_hamburger` | [402](../../../src/realizability.cpp#L402) |
| `check_stieltjes` | [431](../../../src/realizability.cpp#L431) |
| `check_hausdorff` | [488](../../../src/realizability.cpp#L488) |
| `with_support_context` | [543](../../../src/realizability.cpp#L543) |
| `check_realizability` | [554](../../../src/realizability.cpp#L554), [577](../../../src/realizability.cpp#L577) |

### `src/reconstruction_scaling.cpp`

[Source](../../../src/reconstruction_scaling.cpp) · [Responsibilities, contracts, and evidence](qmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `finite_positive` | [14](../../../src/reconstruction_scaling.cpp#L14) |
| `convert_to_double` | [19](../../../src/reconstruction_scaling.cpp#L19) |
| `automatic_scale` | [30](../../../src/reconstruction_scaling.cpp#L30) |
| `scale_value` | [54](../../../src/reconstruction_scaling.cpp#L54) |
| `scale_mixed_value` | [68](../../../src/reconstruction_scaling.cpp#L68) |
| `scale_bounds` | [89](../../../src/reconstruction_scaling.cpp#L89) |
| `make_qmom_scaling_context` | [102](../../../src/reconstruction_scaling.cpp#L102) |
| `make_cqmom_scaling_context` | [121](../../../src/reconstruction_scaling.cpp#L121) |
| `restore_abscissas` | [171](../../../src/reconstruction_scaling.cpp#L171) |
| `restore_cqmom_abscissas` | [186](../../../src/reconstruction_scaling.cpp#L186) |
| `restore_bounds` | [202](../../../src/reconstruction_scaling.cpp#L202) |
| `can_restore_cqmom_result` | [217](../../../src/reconstruction_scaling.cpp#L217) |

### `src/recursive_quadrature.cpp`

[Source](../../../src/recursive_quadrature.cpp) · [Responsibilities, contracts, and evidence](cqmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `IndexedMomentSet::IndexedMomentSet` | [32](../../../src/recursive_quadrature.cpp#L32) |
| `IndexedMomentSet::key_extents` | [56](../../../src/recursive_quadrature.cpp#L56) |
| `IndexedMomentSet::indices` | [61](../../../src/recursive_quadrature.cpp#L61) |
| `IndexedMomentSet::values` | [66](../../../src/recursive_quadrature.cpp#L66) |
| `IndexedMomentSet::mutable_values` | [71](../../../src/recursive_quadrature.cpp#L71) |
| `IndexedMomentSet::size` | [76](../../../src/recursive_quadrature.cpp#L76) |
| `IndexedMomentSet::contains` | [81](../../../src/recursive_quadrature.cpp#L81) |
| `IndexedMomentSet::offset` | [98](../../../src/recursive_quadrature.cpp#L98) |
| `IndexedMomentSet::moment` | [111](../../../src/recursive_quadrature.cpp#L111) |
| `CqmomMomentView::CqmomMomentView` | [116](../../../src/recursive_quadrature.cpp#L116), [118](../../../src/recursive_quadrature.cpp#L118) |
| `CqmomMomentView::moment` | [120](../../../src/recursive_quadrature.cpp#L120) |
| `create_indexed_cqmom_moment_set` | [129](../../../src/recursive_quadrature.cpp#L129) |
| `create_ordered_root_layout` | [137](../../../src/recursive_quadrature.cpp#L137) |
| `append_unique_index` | [164](../../../src/recursive_quadrature.cpp#L164) |
| `local_branch_offset` | [172](../../../src/recursive_quadrature.cpp#L172) |
| `create_child_layout_case` | [182](../../../src/recursive_quadrature.cpp#L182) |
| `create_branch_moment_layout` | [228](../../../src/recursive_quadrature.cpp#L228) |
| `PreparedFixedCqmomMetadata::PreparedFixedCqmomMetadata` | [275](../../../src/recursive_quadrature.cpp#L275) |
| `create_prepared_fixed_cqmom_metadata` | [285](../../../src/recursive_quadrature.cpp#L285) |
| `prepared_fixed_cqmom_metadata_is_compatible` | [295](../../../src/recursive_quadrature.cpp#L295) |
| `extract_wheeler_moments` | [310](../../../src/recursive_quadrature.cpp#L310) |
| `child_layout_for_parent_node_count` | [329](../../../src/recursive_quadrature.cpp#L329) |
| `make_branch_rhs_matrix` | [343](../../../src/recursive_quadrature.cpp#L343) |
| `validate_positive_extents` | [372](../../../src/recursive_quadrature.cpp#L372) |
| `checked_multiply` | [384](../../../src/recursive_quadrature.cpp#L384) |
| `default_conditioning_order` | [393](../../../src/recursive_quadrature.cpp#L393) |
| `validate_conditioning_order` | [400](../../../src/recursive_quadrature.cpp#L400) |
| `ordered_extents` | [422](../../../src/recursive_quadrature.cpp#L422) |
| `from_ordered_index` | [433](../../../src/recursive_quadrature.cpp#L433) |
| `validate_abscissas` | [443](../../../src/recursive_quadrature.cpp#L443) |
| `valid_conditional_options` | [484](../../../src/recursive_quadrature.cpp#L484) |
| `valid_diagnostic_level` | [490](../../../src/recursive_quadrature.cpp#L490) |
| `records_success_branch_diagnostics` | [502](../../../src/recursive_quadrature.cpp#L502) |
| `create_ordered_root_slice` | [507](../../../src/recursive_quadrature.cpp#L507) |
| `one_dimensional_moments` | [522](../../../src/recursive_quadrature.cpp#L522) |
| `reconstruct_child_slices` | [527](../../../src/recursive_quadrature.cpp#L527) |
| `create_recursive_result` | [618](../../../src/recursive_quadrature.cpp#L618) |
| `record_failure` | [623](../../../src/recursive_quadrature.cpp#L623) |
| `validate_recursive_request_preflight` | [632](../../../src/recursive_quadrature.cpp#L632) |
| `resolve_cqmom_wheeler_options` | [656](../../../src/recursive_quadrature.cpp#L656) |
| `operator<` | [672](../../../src/recursive_quadrature.cpp#L672) |
| `path_has_prefix` | [697](../../../src/recursive_quadrature.cpp#L697) |
| `branch_limit_invalidated` | [713](../../../src/recursive_quadrature.cpp#L713) |
| `validate_active_branch_topology` | [725](../../../src/recursive_quadrature.cpp#L725) |
| `recursive_branch_node_count_limit` | [797](../../../src/recursive_quadrature.cpp#L797) |
| `compute_branch_quadrature` | [815](../../../src/recursive_quadrature.cpp#L815) |
| `validate_coordinate_support` | [820](../../../src/recursive_quadrature.cpp#L820) |
| `support_interval_for_dimension` | [847](../../../src/recursive_quadrature.cpp#L847) |
| `record_bounded_candidate_diagnostics` | [857](../../../src/recursive_quadrature.cpp#L857) |
| `branch_growth_trial_allowed` | [879](../../../src/recursive_quadrature.cpp#L879) |
| `product_extents` | [889](../../../src/recursive_quadrature.cpp#L889) |
| `flatten_multi_index` | [901](../../../src/recursive_quadrature.cpp#L901) |
| `unflatten_multi_index` | [921](../../../src/recursive_quadrature.cpp#L921) |
| `enumerate_multi_indices` | [940](../../../src/recursive_quadrature.cpp#L940) |
| `PhaseSpaceHierarchy::PhaseSpaceHierarchy` | [952](../../../src/recursive_quadrature.cpp#L952), [959](../../../src/recursive_quadrature.cpp#L959) |
| `PhaseSpaceHierarchy::number_of_dimensions` | [965](../../../src/recursive_quadrature.cpp#L965) |
| `PhaseSpaceHierarchy::node_counts` | [970](../../../src/recursive_quadrature.cpp#L970) |
| `PhaseSpaceHierarchy::conditioning_order` | [975](../../../src/recursive_quadrature.cpp#L975) |
| `PhaseSpaceHierarchy::node_count` | [980](../../../src/recursive_quadrature.cpp#L980) |
| `PhaseSpaceHierarchy::total_nodes` | [989](../../../src/recursive_quadrature.cpp#L989) |
| `PhaseSpaceHierarchy::moment_extents` | [994](../../../src/recursive_quadrature.cpp#L994) |
| `PhaseSpaceHierarchy::flatten_node_index` | [1005](../../../src/recursive_quadrature.cpp#L1005) |
| `PhaseSpaceHierarchy::unflatten_node_index` | [1014](../../../src/recursive_quadrature.cpp#L1014) |
| `create_cqmom_moment_layout` | [1019](../../../src/recursive_quadrature.cpp#L1019) |
| `CqmomMomentData::CqmomMomentData` | [1063](../../../src/recursive_quadrature.cpp#L1063) |
| `CqmomMomentData::contains_required_index` | [1068](../../../src/recursive_quadrature.cpp#L1068) |
| `CqmomMomentData::moment` | [1073](../../../src/recursive_quadrature.cpp#L1073) |
| `extract_cqmom_moments` | [1089](../../../src/recursive_quadrature.cpp#L1089) |
| `create_cqmom_moment_data` | [1116](../../../src/recursive_quadrature.cpp#L1116) |
| `enumerate_quadrature_node_multi_indices` | [1138](../../../src/recursive_quadrature.cpp#L1138) |
| `phase_space_node_index_to_path` | [1153](../../../src/recursive_quadrature.cpp#L1153) |
| `rectangular_terminal_paths` | [1164](../../../src/recursive_quadrature.cpp#L1164) |
| `RecursiveQuadrature::RecursiveQuadrature` | [1180](../../../src/recursive_quadrature.cpp#L1180), [1187](../../../src/recursive_quadrature.cpp#L1187) |
| `RecursiveQuadrature::hierarchy` | [1198](../../../src/recursive_quadrature.cpp#L1198) |
| `RecursiveQuadrature::number_of_dimensions` | [1203](../../../src/recursive_quadrature.cpp#L1203) |
| `RecursiveQuadrature::total_nodes` | [1208](../../../src/recursive_quadrature.cpp#L1208) |
| `RecursiveQuadrature::weights` | [1213](../../../src/recursive_quadrature.cpp#L1213) |
| `RecursiveQuadrature::abscissas` | [1218](../../../src/recursive_quadrature.cpp#L1218) |
| `RecursiveQuadrature::adaptive_branch_counts_enabled` | [1223](../../../src/recursive_quadrature.cpp#L1223) |
| `RecursiveQuadrature::has_rectangular_storage` | [1228](../../../src/recursive_quadrature.cpp#L1228) |
| `RecursiveQuadrature::tree_nodes` | [1233](../../../src/recursive_quadrature.cpp#L1233) |
| `RecursiveQuadrature::terminal_node_ids` | [1238](../../../src/recursive_quadrature.cpp#L1238) |
| `RecursiveQuadrature::terminal_paths` | [1243](../../../src/recursive_quadrature.cpp#L1243) |
| `RecursiveQuadrature::weight` | [1248](../../../src/recursive_quadrature.cpp#L1248) |
| `RecursiveQuadrature::abscissa` | [1257](../../../src/recursive_quadrature.cpp#L1257) |
| `RecursiveQuadrature::node` | [1270](../../../src/recursive_quadrature.cpp#L1270) |
| `RecursiveQuadrature::set_weight` | [1280](../../../src/recursive_quadrature.cpp#L1280) |
| `RecursiveQuadrature::set_abscissa` | [1289](../../../src/recursive_quadrature.cpp#L1289) |
| `RecursiveQuadrature::set_adaptive_branch_counts_enabled` | [1302](../../../src/recursive_quadrature.cpp#L1302) |
| `RecursiveQuadrature::add_tree_node` | [1315](../../../src/recursive_quadrature.cpp#L1315) |
| `RecursiveQuadrature::set_terminal_quadrature` | [1327](../../../src/recursive_quadrature.cpp#L1327) |
| `compute_recursive_quadrature_from_moment_view` | [1353](../../../src/recursive_quadrature.cpp#L1353) |
| `compute_recursive_quadrature_from_validated_moment_view` | [1676](../../../src/recursive_quadrature.cpp#L1676) |
| `compute_recursive_quadrature_from_moment_data` | [1686](../../../src/recursive_quadrature.cpp#L1686) |
| `compute_recursive_quadrature_from_trusted_moment_data` | [1721](../../../src/recursive_quadrature.cpp#L1721) |
| `compute_recursive_quadrature` | [1754](../../../src/recursive_quadrature.cpp#L1754) |

### `src/runge_kutta.cpp`

[Source](../../../src/runge_kutta.cpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `require_same_size` | [21](../../../src/runge_kutta.cpp#L21) |
| `euclidean_norm` | [29](../../../src/runge_kutta.cpp#L29) |
| `rms_norm` | [39](../../../src/runge_kutta.cpp#L39) |
| `rhs_checked` | [48](../../../src/runge_kutta.cpp#L48) |
| `add_scaled` | [59](../../../src/runge_kutta.cpp#L59) |
| `sign` | [70](../../../src/runge_kutta.cpp#L70) |
| `validate_or_default_max_step` | [75](../../../src/runge_kutta.cpp#L75) |
| `minimum_step` | [84](../../../src/runge_kutta.cpp#L84) |
| `error_norm` | [100](../../../src/runge_kutta.cpp#L100) |
| `invalid_embedded_rk_result` | [112](../../../src/runge_kutta.cpp#L112) |
| `stage_accepted` | [117](../../../src/runge_kutta.cpp#L117) |
| `try_ssp_rk3_step` | [122](../../../src/runge_kutta.cpp#L122) |
| `ssp_rk3_step` | [179](../../../src/runge_kutta.cpp#L179) |
| `integrate_embedded_rk` | [224](../../../src/runge_kutta.cpp#L224) |

### `src/wheeler.cpp`

[Source](../../../src/wheeler.cpp) · [Responsibilities, contracts, and evidence](qmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `MomentSet::MomentSet` | [23](../../../src/wheeler.cpp#L23), [25](../../../src/wheeler.cpp#L25) |
| `MomentSet::size` | [27](../../../src/wheeler.cpp#L27) |
| `MomentSet::empty` | [32](../../../src/wheeler.cpp#L32) |
| `MomentSet::values` | [37](../../../src/wheeler.cpp#L37) |
| `MomentSet::operator[]` | [42](../../../src/wheeler.cpp#L42) |
| `one_node_result` | [70](../../../src/wheeler.cpp#L70) |
| `negative_mass_result` | [83](../../../src/wheeler.cpp#L83) |
| `compute_recurrence` | [91](../../../src/wheeler.cpp#L91) |
| `make_jacobi_matrix` | [121](../../../src/wheeler.cpp#L121) |
| `symmetric_jacobi_eigen` | [142](../../../src/wheeler.cpp#L142) |
| `recurrence_is_finite` | [151](../../../src/wheeler.cpp#L151) |
| `has_negative_recurrence_b` | [156](../../../src/wheeler.cpp#L156) |
| `adaptive_sigma_tolerance` | [161](../../../src/wheeler.cpp#L161) |
| `compute_weights` | [172](../../../src/wheeler.cpp#L172) |
| `adaptive_accepts` | [184](../../../src/wheeler.cpp#L184) |
| `required_moments_are_finite` | [225](../../../src/wheeler.cpp#L225) |
| `attempt_quadrature_count` | [237](../../../src/wheeler.cpp#L237) |
| `compute_quadrature` | [320](../../../src/wheeler.cpp#L320) |
| `select_bounded_quadrature` | [392](../../../src/wheeler.cpp#L392) |

## Declaration-only headers and aliases

These headers introduce declarations, aliases, option/result structures, or umbrella
includes rather than additional explicit executable bodies in the inspected build.
Their callable definitions are listed above. Aggregate initialization is covered
with the owning representation or options contract.

| Header | Functional owner |
|---|---|
| [include/qmomkit/abscissa_bounds.hpp](../../../include/qmomkit/abscissa_bounds.hpp) | [qmom reconstruction](qmom_reconstruction.md) |
| [include/qmomkit/conditioning_path.hpp](../../../include/qmomkit/conditioning_path.hpp) | [state and representations](state_and_representations.md) |
| [include/qmomkit/moment_reconstruction.hpp](../../../include/qmomkit/moment_reconstruction.hpp) | [moment evaluation](moment_evaluation.md) |
| [include/qmomkit/moment_set.hpp](../../../include/qmomkit/moment_set.hpp) | [state and representations](state_and_representations.md) |
| [include/qmomkit/multi_index.hpp](../../../include/qmomkit/multi_index.hpp) | [state and representations](state_and_representations.md) |
| [include/qmomkit/output.hpp](../../../include/qmomkit/output.hpp) | [utilities](utilities.md) |
| [include/qmomkit/phase_space_hierarchy.hpp](../../../include/qmomkit/phase_space_hierarchy.hpp) | [state and representations](state_and_representations.md) |
| [include/qmomkit/qmomkit.hpp](../../../include/qmomkit/qmomkit.hpp) | [README](README.md) |
| [include/qmomkit/realizability.hpp](../../../include/qmomkit/realizability.hpp) | [realizability](realizability.md) |
| [include/qmomkit/vandermonde_solver.hpp](../../../include/qmomkit/vandermonde_solver.hpp) | [conditional systems](conditional_systems.md) |
| [include/qmomkit/wheeler.hpp](../../../include/qmomkit/wheeler.hpp) | [qmom reconstruction](qmom_reconstruction.md) |
| [src/internal/bjorck_pereyra.hpp](../../../src/internal/bjorck_pereyra.hpp) | [conditional systems](conditional_systems.md) |
| [src/internal/bounded_quadrature.hpp](../../../src/internal/bounded_quadrature.hpp) | [qmom reconstruction](qmom_reconstruction.md) |
| [src/internal/cqmom_case.hpp](../../../src/internal/cqmom_case.hpp) | [state and representations](state_and_representations.md) |
| [src/internal/embedded_runge_kutta.hpp](../../../src/internal/embedded_runge_kutta.hpp) | [time integration](time_integration.md) |
| [src/internal/householder_qr.hpp](../../../src/internal/householder_qr.hpp) | [conditional systems](conditional_systems.md) |
| [src/internal/jacobi_svd.hpp](../../../src/internal/jacobi_svd.hpp) | [conditional systems](conditional_systems.md) |
| [src/internal/moment_reconstruction.hpp](../../../src/internal/moment_reconstruction.hpp) | [moment evaluation](moment_evaluation.md) |
| [src/internal/moment_source.hpp](../../../src/internal/moment_source.hpp) | [moment sources](moment_sources.md) |
| [src/internal/reconstruction_scaling.hpp](../../../src/internal/reconstruction_scaling.hpp) | [qmom reconstruction](qmom_reconstruction.md) |
| [src/internal/vandermonde.hpp](../../../src/internal/vandermonde.hpp) | [conditional systems](conditional_systems.md) |
