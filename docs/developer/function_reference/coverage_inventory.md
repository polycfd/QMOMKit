# Definition coverage and navigation inventory

This inventory routes core implementation definitions to functional explanations.
It complements the algorithm sections and responsibility tables; a symbol listing
does not by itself establish intent or test coverage.

## Scope and method

The inventory covers the 23 core `.cpp` files and installed/private headers.
It lists explicit functions, methods, constructors, and operators by defining
source range. Defaulted or implicit special members are omitted; consequential
lambdas are documented with their owning function. Instrumentation definitions
are included for its enabled preprocessor branch, and debug-only contract checks
for theirs. Line links are navigation aids, not permanent API identifiers: they
drift whenever a file is edited above a listed definition, so treat the file and
the name as authoritative and the line as a hint.

### Checking and refreshing the line numbers

Because the line numbers drift, verify them rather than trusting them. This
reports every link whose line no longer contains the named definition:

```bash
python3 - <<'EOF'
import re, os
inv = open('docs/developer/function_reference/coverage_inventory.md').read()
for m in re.finditer(r'\| `([^`]+)` \| ((?:\[\d+\]\([^)]+\),?\s*)+)\|', inv):
    name = m.group(1).split('::')[-1]
    for line, path in re.findall(r'\[(\d+)\]\(\.\./\.\./\.\./([^)#]+)#L\d+\)', m.group(2)):
        src = open(path).read().split('\n') if os.path.exists(path) else []
        if int(line) > len(src) or name not in src[int(line)-1]:
            print(f'{path}:{line} no longer defines {name}')
EOF
```

A clean run prints nothing. When it reports drift, the nearest definition of the
same name in the same file is almost always the new location; when it reports a
name the file no longer contains at all, the definition moved or was removed and
the row needs a decision rather than a new number.

The inventory includes **943 definition locations**, grouped by named function
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
| `ConditionalMomentResult::ok` | [70](../../../include/qmomkit/conditional_moments.hpp#L70) |

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
| `CqmomQuadrature::number_of_dimensions` | [48](../../../include/qmomkit/cqmom.hpp#L48) |
| `CqmomQuadrature::total_nodes` | [52](../../../include/qmomkit/cqmom.hpp#L52) |
| `CqmomQuadrature::weights` | [56](../../../include/qmomkit/cqmom.hpp#L56) |
| `CqmomQuadrature::abscissas` | [60](../../../include/qmomkit/cqmom.hpp#L60) |
| `CqmomQuadrature::adaptive_branch_counts_enabled` | [64](../../../include/qmomkit/cqmom.hpp#L64) |
| `CqmomQuadrature::has_rectangular_storage` | [68](../../../include/qmomkit/cqmom.hpp#L68) |
| `CqmomQuadrature::terminal_paths` | [72](../../../include/qmomkit/cqmom.hpp#L72) |
| `CqmomResult::ok` | [335](../../../include/qmomkit/cqmom.hpp#L335) |

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
| `MomentSourceUpdateResult::ok` | [259](../../../include/qmomkit/moment_source.hpp#L259) |
| `CqmomMomentSourceUpdateResult::ok` | [259](../../../include/qmomkit/moment_source.hpp#L259) |

### `include/qmomkit/node_advection.hpp`

[Source](../../../include/qmomkit/node_advection.hpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomNodeAdvectionResult::ok` | [101](../../../include/qmomkit/node_advection.hpp#L101) |
| `CqmomNodeAdvectionResult::ok` | [101](../../../include/qmomkit/node_advection.hpp#L101) |

### `include/qmomkit/phase_space_evolution.hpp`

[Source](../../../include/qmomkit/phase_space_evolution.hpp) · [Responsibilities, contracts, and evidence](phase_space_evolution.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomPhaseSpaceEvolutionResult::ok` | [143](../../../include/qmomkit/phase_space_evolution.hpp#L143) |
| `CqmomPhaseSpaceEvolutionResult::ok` | [158](../../../include/qmomkit/phase_space_evolution.hpp#L158) |

### `include/qmomkit/phase_space_hierarchy.hpp`

[Source](../../../include/qmomkit/phase_space_hierarchy.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `PhaseSpaceHierarchy::node_counts` | [32](../../../include/qmomkit/phase_space_hierarchy.hpp#L32) |
| `PhaseSpaceHierarchy::conditioning_order` | [36](../../../include/qmomkit/phase_space_hierarchy.hpp#L36) |

### `include/qmomkit/qmom.hpp`

[Source](../../../include/qmomkit/qmom.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomResult::ok` | [40](../../../include/qmomkit/qmom.hpp#L40) |
| `QmomState::moments` | [70](../../../include/qmomkit/qmom.hpp#L70) |
| `QmomState::quadrature` | [75](../../../include/qmomkit/qmom.hpp#L75) |
| `QmomState::active_node_count` | [80](../../../include/qmomkit/qmom.hpp#L80) |
| `QmomState::maximum_node_count` | [85](../../../include/qmomkit/qmom.hpp#L85) |
| `QmomState::qmom_options` | [90](../../../include/qmomkit/qmom.hpp#L90) |
| `QmomState::abscissa_bounds` | [95](../../../include/qmomkit/qmom.hpp#L95) |

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
| `SspRk3Result::ok` | [84](../../../include/qmomkit/runge_kutta.hpp#L84) |

### `src/conditional_moments.cpp`

[Source](../../../src/conditional_moments.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `invalid_result` | [25](../../../src/conditional_moments.cpp#L25) |
| `valid_options` | [40](../../../src/conditional_moments.cpp#L40) |
| `extract_rhs` | [47](../../../src/conditional_moments.cpp#L47) |
| `ConditionalMomentResult::moment` | [72](../../../src/conditional_moments.cpp#L72) |
| `reconstruct_conditional_moments` | [81](../../../src/conditional_moments.cpp#L81) |

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
| `ConditioningPath::phase_space_node_index` | [68](../../../src/conditioning_path.cpp#L68) |
| `ConditioningPath::to_string` | [73](../../../src/conditioning_path.cpp#L73) |
| `ConditioningPath::describe` | [89](../../../src/conditioning_path.cpp#L89) |
| `ConditioningPathView::phase_space_node_index` | [47](../../../src/conditioning_path.cpp#L47) |
| `ConditioningPathView::node_index_at_level` | [109](../../../src/conditioning_path.cpp#L109) |
| `ConditioningPathTable::ConditioningPathTable` | [128](../../../src/conditioning_path.cpp#L128), [145](../../../src/conditioning_path.cpp#L145) |
| `ConditioningPathTable::at` | [164](../../../src/conditioning_path.cpp#L164) |
| `ConditioningPathTable::push_back` | [173](../../../src/conditioning_path.cpp#L173), [178](../../../src/conditioning_path.cpp#L178) |
| `ConditioningPathTable::clear` | [191](../../../src/conditioning_path.cpp#L191) |
| `ConditioningPathTable::reserve` | [196](../../../src/conditioning_path.cpp#L196) |
| `ConditioningPathTable::to_conditioning_paths` | [201](../../../src/conditioning_path.cpp#L201) |
| `operator==` | [118](../../../src/conditioning_path.cpp#L118), [212](../../../src/conditioning_path.cpp#L212), [222](../../../src/conditioning_path.cpp#L222) |
| `operator!=` | [123](../../../src/conditioning_path.cpp#L123), [217](../../../src/conditioning_path.cpp#L217), [227](../../../src/conditioning_path.cpp#L227) |
| `operator<` | [232](../../../src/conditioning_path.cpp#L232) |
| `recursive_path_to_node_index` | [238](../../../src/conditioning_path.cpp#L238) |

### `src/configuration.cpp`

[Source](../../../src/configuration.cpp) · [Responsibilities, contracts, and evidence](initialization_and_configuration.md)

| Function or overload family | Definition lines |
|---|---|
| `trim` | [55](../../../src/configuration.cpp#L55) |
| `add_error` | [66](../../../src/configuration.cpp#L66) |
| `parse_size` | [72](../../../src/configuration.cpp#L72) |
| `parse_coordinate_realizability_section` | [84](../../../src/configuration.cpp#L84) |
| `is_coordinate_realizability_key` | [93](../../../src/configuration.cpp#L93) |
| `is_configured_coordinate_realizability_key` | [108](../../../src/configuration.cpp#L108) |
| `parse_double` | [118](../../../src/configuration.cpp#L118) |
| `parse_size_list` | [132](../../../src/configuration.cpp#L132) |
| `parse_double_list` | [148](../../../src/configuration.cpp#L148) |
| `parse_optional_double_list` | [164](../../../src/configuration.cpp#L164) |
| `validate_names` | [186](../../../src/configuration.cpp#L186) |
| `parse_moment_key` | [207](../../../src/configuration.cpp#L207) |
| `mixture_component_prefix` | [228](../../../src/configuration.cpp#L228) |
| `is_mixture_component_section` | [233](../../../src/configuration.cpp#L233) |
| `is_independent_product_dimension_section` | [244](../../../src/configuration.cpp#L244) |
| `parse_file` | [254](../../../src/configuration.cpp#L254) |
| `require_entry` | [331](../../../src/configuration.cpp#L331) |
| `find_entry` | [342](../../../src/configuration.cpp#L342) |
| `reject_unknown_keys` | [352](../../../src/configuration.cpp#L352) |
| `reject_unknown_sections` | [370](../../../src/configuration.cpp#L370) |
| `validate_format_version` | [389](../../../src/configuration.cpp#L389) |
| `parse_phase_space_node_counts` | [398](../../../src/configuration.cpp#L398) |
| `validate_phase_space_names` | [410](../../../src/configuration.cpp#L410) |
| `validate_qmom_phase_space` | [420](../../../src/configuration.cpp#L420) |
| `validate_cqmom_phase_space` | [433](../../../src/configuration.cpp#L433) |
| `validate_optional_qmom_phase_space` | [453](../../../src/configuration.cpp#L453) |
| `validate_optional_cqmom_phase_space` | [469](../../../src/configuration.cpp#L469) |
| `parse_required_double` | [485](../../../src/configuration.cpp#L485) |
| `parse_optional_double` | [500](../../../src/configuration.cpp#L500) |
| `parse_optional_bool` | [515](../../../src/configuration.cpp#L515) |
| `parse_optional_size` | [536](../../../src/configuration.cpp#L536) |
| `parse_optional_positive_double` | [551](../../../src/configuration.cpp#L551) |
| `infer_qmom_initial_maximum_order` | [568](../../../src/configuration.cpp#L568) |
| `parse_optional_positive_size` | [595](../../../src/configuration.cpp#L595) |
| `parse_distribution` | [612](../../../src/configuration.cpp#L612), [699](../../../src/configuration.cpp#L699) |
| `component_key_allowed` | [704](../../../src/configuration.cpp#L704) |
| `marginal_key_allowed` | [716](../../../src/configuration.cpp#L716) |
| `has_explicit_qmom_moments` | [732](../../../src/configuration.cpp#L732) |
| `has_explicit_cqmom_moments` | [737](../../../src/configuration.cpp#L737) |
| `qmom_allowed_keys` | [742](../../../src/configuration.cpp#L742) |
| `append_time_integration_allowed_keys` | [767](../../../src/configuration.cpp#L767) |
| `append_output_allowed_keys` | [786](../../../src/configuration.cpp#L786) |
| `with_runtime_allowed_keys` | [795](../../../src/configuration.cpp#L795) |
| `time_integration_sections` | [802](../../../src/configuration.cpp#L802) |
| `reject_unknown_qmom_sections` | [807](../../../src/configuration.cpp#L807) |
| `reject_unknown_qmom_keys` | [820](../../../src/configuration.cpp#L820) |
| `parse_mixture` | [847](../../../src/configuration.cpp#L847) |
| `has_mixture_component_sections` | [916](../../../src/configuration.cpp#L916) |
| `cqmom_explicit_allowed_keys` | [928](../../../src/configuration.cpp#L928) |
| `cqmom_initial_allowed_keys` | [936](../../../src/configuration.cpp#L936) |
| `reject_unknown_cqmom_sections` | [945](../../../src/configuration.cpp#L945) |
| `reject_unknown_cqmom_keys` | [958](../../../src/configuration.cpp#L958) |
| `reject_duplicate_indexed_section` | [838](../../../src/configuration.cpp#L838) |
| `parse_independent_product` | [968](../../../src/configuration.cpp#L968) |
| `checked_product` | [1061](../../../src/configuration.cpp#L1061) |
| `parse_explicit_tensor` | [1072](../../../src/configuration.cpp#L1072) |
| `format_configuration_error` | [1134](../../../src/configuration.cpp#L1134) |
| `format_configuration_errors` | [1150](../../../src/configuration.cpp#L1150) |
| `failure_result` | [1159](../../../src/configuration.cpp#L1159), [1164](../../../src/configuration.cpp#L1164) |
| `file_failure` | [1169](../../../src/configuration.cpp#L1169) |
| `validate_replacement_tensor_structure` | [1174](../../../src/configuration.cpp#L1174) |
| `format_moment_index` | [1197](../../../src/configuration.cpp#L1197) |
| `read_configuration_file` | [1213](../../../src/configuration.cpp#L1213) |
| `has_realizability_section` | [1239](../../../src/configuration.cpp#L1239) |
| `has_unnumbered_realizability_section` | [1249](../../../src/configuration.cpp#L1249) |
| `coordinate_realizability_sections` | [1259](../../../src/configuration.cpp#L1259) |
| `is_realizability_key` | [1273](../../../src/configuration.cpp#L1273) |
| `realizability_policy_fields` | [1278](../../../src/configuration.cpp#L1278) |
| `is_node_advection_key` | [1284](../../../src/configuration.cpp#L1284) |
| `node_advection_allowed_keys` | [1289](../../../src/configuration.cpp#L1289) |
| `parse_qmom_options_from_parsed` | [1298](../../../src/configuration.cpp#L1298) |
| `parse_automatic_scaling_from_parsed` | [1326](../../../src/configuration.cpp#L1326) |
| `reject_unknown_qmom_option_keys` | [1345](../../../src/configuration.cpp#L1345) |
| `parse_realizability_check_mode_value` | [1363](../../../src/configuration.cpp#L1363) |
| `parse_realizability_check_mode` | [1382](../../../src/configuration.cpp#L1382) |
| `parse_moment_support` | [1393](../../../src/configuration.cpp#L1393) |
| `parse_phase_space_evolution_scheme_value` | [1416](../../../src/configuration.cpp#L1416) |
| `parse_node_advection_scheme_value` | [1435](../../../src/configuration.cpp#L1435) |
| `parse_node_advection_options_from_parsed` | [1454](../../../src/configuration.cpp#L1454) |
| `parse_realizability_policy_from_parsed` | [1502](../../../src/configuration.cpp#L1502) |
| `parse_coordinate_realizability_policies_from_parsed` | [1583](../../../src/configuration.cpp#L1583) |
| `parse_cqmom_diagnostic_level_value` | [1611](../../../src/configuration.cpp#L1611) |
| `parse_cqmom_diagnostic_level` | [1629](../../../src/configuration.cpp#L1629) |
| `parse_cqmom_moment_storage_value` | [1639](../../../src/configuration.cpp#L1639) |
| `parse_cqmom_moment_storage` | [1653](../../../src/configuration.cpp#L1653) |
| `parse_vandermonde_solver_value` | [1663](../../../src/configuration.cpp#L1663) |
| `parse_vandermonde_solver` | [1712](../../../src/configuration.cpp#L1712) |
| `has_output_section` | [1722](../../../src/configuration.cpp#L1722) |
| `is_output_key` | [1732](../../../src/configuration.cpp#L1732) |
| `validate_optional_phase_space_evolution_scheme` | [1737](../../../src/configuration.cpp#L1737) |
| `validate_optional_node_advection_options` | [1746](../../../src/configuration.cpp#L1746) |
| `validate_optional_realizability_policy` | [1751](../../../src/configuration.cpp#L1751) |
| `validate_cqmom_coordinate_realizability_ranks` | [1763](../../../src/configuration.cpp#L1763) |
| `reject_unknown_phase_space_evolution_sections` | [1777](../../../src/configuration.cpp#L1777) |
| `reject_unknown_phase_space_evolution_keys` | [1790](../../../src/configuration.cpp#L1790) |
| `parse_phase_space_evolution_scheme_from_parsed` | [1809](../../../src/configuration.cpp#L1809) |
| `validate_optional_cqmom_diagnostic_level` | [1829](../../../src/configuration.cpp#L1829) |
| `validate_optional_cqmom_moment_storage` | [1834](../../../src/configuration.cpp#L1834) |
| `validate_optional_vandermonde_solver` | [1839](../../../src/configuration.cpp#L1839) |
| `parse_coordinate_support` | [1845](../../../src/configuration.cpp#L1845) |
| `parse_cqmom_case_settings` | [1907](../../../src/configuration.cpp#L1907) |
| `parse_qmom_initial_condition_from_parsed` | [2027](../../../src/configuration.cpp#L2027) |
| `parse_cqmom_initial_condition_from_parsed` | [2122](../../../src/configuration.cpp#L2122) |
| `initialize_moments` | [2236](../../../src/configuration.cpp#L2236) |
| `initialize_cqmom_dense_moments` | [2271](../../../src/configuration.cpp#L2271) |
| `initialize_cqmom_moment_data` | [2301](../../../src/configuration.cpp#L2301) |
| `parse_qmom_initial_condition` | [2348](../../../src/configuration.cpp#L2348) |
| `parse_cqmom_initial_condition` | [2354](../../../src/configuration.cpp#L2354) |
| `create_qmom_case` | [2360](../../../src/configuration.cpp#L2360) |
| `create_cqmom_dense_case` | [2387](../../../src/configuration.cpp#L2387) |
| `create_cqmom_lean_case` | [2444](../../../src/configuration.cpp#L2444) |
| `with_initial_moments` | [2463](../../../src/configuration.cpp#L2463) |
| `parse_qmom` | [2519](../../../src/configuration.cpp#L2519) |
| `parse_cqmom_dense` | [2574](../../../src/configuration.cpp#L2574) |
| `parse_cqmom_lean` | [2605](../../../src/configuration.cpp#L2605) |
| `parse_cqmom` | [2636](../../../src/configuration.cpp#L2636) |
| `parse_realizability_policy` | [2679](../../../src/configuration.cpp#L2679) |
| `parse_coordinate_realizability_policies` | [2708](../../../src/configuration.cpp#L2708) |
| `parse_phase_space_evolution_scheme` | [2728](../../../src/configuration.cpp#L2728) |
| `parse_node_advection_options` | [2747](../../../src/configuration.cpp#L2747) |
| `parse_qmom_phase_space_evolution_options` | [2766](../../../src/configuration.cpp#L2766) |
| `parse_cqmom_phase_space_evolution_options` | [2796](../../../src/configuration.cpp#L2796) |
| `parse_output_options` | [2822](../../../src/configuration.cpp#L2822) |
| `load_configuration` | [2928](../../../src/configuration.cpp#L2928) |
| `load_qmom` | [2940](../../../src/configuration.cpp#L2940) |
| `load_cqmom_dense` | [2945](../../../src/configuration.cpp#L2945) |
| `load_cqmom_lean` | [2950](../../../src/configuration.cpp#L2950) |
| `load_cqmom` | [2955](../../../src/configuration.cpp#L2955) |
| `load_output_options` | [2960](../../../src/configuration.cpp#L2960) |
| `load_node_advection_options` | [2965](../../../src/configuration.cpp#L2965) |
| `load_realizability_policy` | [2970](../../../src/configuration.cpp#L2970) |
| `load_qmom_phase_space_evolution_options` | [2975](../../../src/configuration.cpp#L2975) |
| `load_cqmom_phase_space_evolution_options` | [2980](../../../src/configuration.cpp#L2980) |

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
| `compute_normal_selected` | [349](../../../src/distributions.cpp#L349) |
| `compute_log_normal_selected` | [383](../../../src/distributions.cpp#L383) |
| `compute_beta_selected` | [416](../../../src/distributions.cpp#L416) |
| `compute_gamma_selected` | [443](../../../src/distributions.cpp#L443) |
| `compute_chi_selected` | [471](../../../src/distributions.cpp#L471) |
| `compute_chi_squared_selected` | [499](../../../src/distributions.cpp#L499) |
| `compute_weibull_selected` | [527](../../../src/distributions.cpp#L527) |
| `compute_normal` | [554](../../../src/distributions.cpp#L554) |
| `compute_log_normal` | [583](../../../src/distributions.cpp#L583) |
| `compute_beta` | [601](../../../src/distributions.cpp#L601) |
| `compute_gamma` | [613](../../../src/distributions.cpp#L613) |
| `compute_chi` | [626](../../../src/distributions.cpp#L626) |
| `compute_chi_squared` | [639](../../../src/distributions.cpp#L639) |
| `compute_weibull` | [653](../../../src/distributions.cpp#L653) |
| `RawMomentVisitor::operator()` | [670](../../../src/distributions.cpp#L670), [675](../../../src/distributions.cpp#L675), [680](../../../src/distributions.cpp#L680), [685](../../../src/distributions.cpp#L685), [690](../../../src/distributions.cpp#L690), [695](../../../src/distributions.cpp#L695), [700](../../../src/distributions.cpp#L700) |
| `SelectedRawMomentVisitor::operator()` | [713](../../../src/distributions.cpp#L713), [718](../../../src/distributions.cpp#L718), [723](../../../src/distributions.cpp#L723), [728](../../../src/distributions.cpp#L728), [733](../../../src/distributions.cpp#L733), [738](../../../src/distributions.cpp#L738), [743](../../../src/distributions.cpp#L743) |
| `add_compensated` | [749](../../../src/distributions.cpp#L749) |
| `add_selected_compensated` | [780](../../../src/distributions.cpp#L780) |
| `contains_order` | [814](../../../src/distributions.cpp#L814) |
| `order_position` | [826](../../../src/distributions.cpp#L826) |
| `compute_independent_product_marginal_moments` | [838](../../../src/distributions.cpp#L838), [848](../../../src/distributions.cpp#L848) |
| `first_requested_index_with_order` | [858](../../../src/distributions.cpp#L858) |
| `make_selected_raw_moment_result` | [872](../../../src/distributions.cpp#L872) |
| `compute_raw_moments` | [916](../../../src/distributions.cpp#L916), [932](../../../src/distributions.cpp#L932), [942](../../../src/distributions.cpp#L942), [1012](../../../src/distributions.cpp#L1012), [1093](../../../src/distributions.cpp#L1093), [1200](../../../src/distributions.cpp#L1200) |

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
| `raw_moment_tensor_value_count` | [30](../../../src/internal/cqmom_validation.hpp#L30) |
| `same_cqmom_layout` | [66](../../../src/internal/cqmom_validation.hpp#L66) |
| `multi_index_within_extents` | [72](../../../src/internal/cqmom_validation.hpp#L72) |
| `validate_cqmom_moment_layout_matches_expected` | [88](../../../src/internal/cqmom_validation.hpp#L88) |
| `validate_fixed_cqmom_moment_layout` | [99](../../../src/internal/cqmom_validation.hpp#L99) |
| `validate_cqmom_moment_values_match_layout` | [105](../../../src/internal/cqmom_validation.hpp#L105) |
| `validate_fixed_cqmom_moment_data` | [117](../../../src/internal/cqmom_validation.hpp#L117) |
| `validate_cqmom_moment_data_matches_hierarchy` | [123](../../../src/internal/cqmom_validation.hpp#L123) |
| `validate_dense_cqmom_moment_tensor` | [131](../../../src/internal/cqmom_validation.hpp#L131) |
| `validate_dense_cqmom_reconstruction_tensor` | [148](../../../src/internal/cqmom_validation.hpp#L148) |
| `validate_cqmom_quadrature_storage` | [163](../../../src/internal/cqmom_validation.hpp#L163) |
| `validate_cqmom_moment_rate` | [193](../../../src/internal/cqmom_validation.hpp#L193) |

### `src/internal/embedded_runge_kutta.cpp`

[Source](../../../src/internal/embedded_runge_kutta.cpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `embedded_rk_tableau` | [47](../../../src/internal/embedded_runge_kutta.cpp#L47) |
| `take_embedded_rk_step` | [59](../../../src/internal/embedded_runge_kutta.cpp#L59) |
| `embedded_rk_controller_exponent` | [112](../../../src/internal/embedded_runge_kutta.cpp#L112) |

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
| `solve_qr` | [281](../../../src/internal/householder_qr.cpp#L281), [348](../../../src/internal/householder_qr.cpp#L348) |

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
| `solve_svd` | [444](../../../src/internal/jacobi_svd.cpp#L444), [514](../../../src/internal/jacobi_svd.cpp#L514) |
| `solve_svd_pseudoinverse` | [524](../../../src/internal/jacobi_svd.cpp#L524) |

### `src/internal/linear_algebra.cpp`

[Source](../../../src/internal/linear_algebra.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_infinity_norm` | [33](../../../src/internal/linear_algebra.cpp#L33) |
| `solve_dense` | [48](../../../src/internal/linear_algebra.cpp#L48) |
| `make_identity` | [124](../../../src/internal/linear_algebra.cpp#L124) |
| `residual_norm` | [134](../../../src/internal/linear_algebra.cpp#L134) |
| `condition_estimate` | [154](../../../src/internal/linear_algebra.cpp#L154) |
| `rotate` | [168](../../../src/internal/linear_algebra.cpp#L168) |
| `symmetric_jacobi_eigen` | [213](../../../src/internal/linear_algebra.cpp#L213) |

### `src/internal/linear_algebra.hpp`

[Source](../../../src/internal/linear_algebra.hpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `DenseMatrix::DenseMatrix` | [18](../../../src/internal/linear_algebra.hpp#L18) |
| `DenseMatrix::rows` | [20](../../../src/internal/linear_algebra.hpp#L20) |
| `DenseMatrix::cols` | [25](../../../src/internal/linear_algebra.hpp#L25) |
| `DenseMatrix::operator()` | [30](../../../src/internal/linear_algebra.hpp#L30), [39](../../../src/internal/linear_algebra.hpp#L39) |

### `src/internal/numeric_utils.hpp`

[Source](../../../src/internal/numeric_utils.hpp) · [Responsibilities, contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `all_finite` | [12](../../../src/internal/numeric_utils.hpp#L12) |
| `output_buffer_poison` | [30](../../../src/internal/numeric_utils.hpp#L30) |
| `is_output_buffer_poison` | [37](../../../src/internal/numeric_utils.hpp#L37) |
| `poison_output_buffer` | [47](../../../src/internal/numeric_utils.hpp#L47) |

### `src/internal/ode_rhs_validation.hpp`

[Source](../../../src/internal/ode_rhs_validation.hpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `ode_rhs_unfilled_message` | [16](../../../src/internal/ode_rhs_validation.hpp#L16) |
| `rhs_checked` | [35](../../../src/internal/ode_rhs_validation.hpp#L35) |

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

### `src/internal/phase_space_velocity_validation.hpp`

[Source](../../../src/internal/phase_space_velocity_validation.hpp) · [Responsibilities, contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `format_diagnostic_vector` | [17](../../../src/internal/phase_space_velocity_validation.hpp#L17) |
| `phase_space_velocity_wrong_dimension_message` | [33](../../../src/internal/phase_space_velocity_validation.hpp#L33) |
| `phase_space_node_slice` | [42](../../../src/internal/phase_space_velocity_validation.hpp#L42) |
| `phase_space_velocity_non_finite_message` | [49](../../../src/internal/phase_space_velocity_validation.hpp#L49) |
| `validate_phase_space_velocity_result` | [62](../../../src/internal/phase_space_velocity_validation.hpp#L62) |
| `phase_space_velocity_unfilled_message` | [85](../../../src/internal/phase_space_velocity_validation.hpp#L85) |
| `validate_packed_phase_space_velocity_result` | [102](../../../src/internal/phase_space_velocity_validation.hpp#L102) |

### `src/internal/qmom_state_access.hpp`

[Source](../../../src/internal/qmom_state_access.hpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomStateAccess::from_quadrature` | [20](../../../src/internal/qmom_state_access.hpp#L20) |

### `src/internal/recursive_quadrature.hpp`

[Source](../../../src/internal/recursive_quadrature.hpp) · [Responsibilities, contracts, and evidence](cqmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `RecursiveQuadratureResult::ok` | [303](../../../src/internal/recursive_quadrature.hpp#L303) |

### `src/internal/vandermonde.cpp`

[Source](../../../src/internal/vandermonde.cpp) · [Responsibilities, contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_all_finite` | [19](../../../src/internal/vandermonde.cpp#L19) |
| `failure_result` | [34](../../../src/internal/vandermonde.cpp#L34) |
| `select_automatic_solver` | [39](../../../src/internal/vandermonde.cpp#L39) |
| `resolve_solver` | [44](../../../src/internal/vandermonde.cpp#L44) |
| `vandermonde_result` | [60](../../../src/internal/vandermonde.cpp#L60) |
| `qr_linear_solve_result` | [67](../../../src/internal/vandermonde.cpp#L67) |
| `svd_linear_solve_result` | [84](../../../src/internal/vandermonde.cpp#L84) |
| `residual_is_accepted` | [101](../../../src/internal/vandermonde.cpp#L101) |
| `make_conditional_vandermonde_matrix` | [114](../../../src/internal/vandermonde.cpp#L114) |
| `solve_conditional_vandermonde_system` | [128](../../../src/internal/vandermonde.cpp#L128) |

### `src/moment_reconstruction.cpp`

[Source](../../../src/moment_reconstruction.cpp) · [Responsibilities, contracts, and evidence](moment_evaluation.md)

| Function or overload family | Definition lines |
|---|---|
| `integer_power` | [20](../../../src/moment_reconstruction.cpp#L20) |
| `validate_reconstructed_value` | [40](../../../src/moment_reconstruction.cpp#L40) |
| `is_integer_exponent` | [48](../../../src/moment_reconstruction.cpp#L48) |
| `validate_generalized_value` | [53](../../../src/moment_reconstruction.cpp#L53) |
| `validate_generalized_exponent` | [61](../../../src/moment_reconstruction.cpp#L61) |
| `generalized_power` | [69](../../../src/moment_reconstruction.cpp#L69) |
| `cqmom_reconstruction_quadrature_messages` | [81](../../../src/moment_reconstruction.cpp#L81) |
| `reconstruct_cqmom_moment_value` | [88](../../../src/moment_reconstruction.cpp#L88) |
| `reconstruct_cqmom_moments_impl` | [107](../../../src/moment_reconstruction.cpp#L107) |
| `reconstruct_raw_moments` | [150](../../../src/moment_reconstruction.cpp#L150), [216](../../../src/moment_reconstruction.cpp#L216), [233](../../../src/moment_reconstruction.cpp#L233) |
| `evaluate_generalized_moment` | [183](../../../src/moment_reconstruction.cpp#L183), [210](../../../src/moment_reconstruction.cpp#L210), [252](../../../src/moment_reconstruction.cpp#L252) |
| `reconstruct_cqmom_moments` | [279](../../../src/moment_reconstruction.cpp#L279) |
| `reconstruct_trusted_cqmom_moments` | [287](../../../src/moment_reconstruction.cpp#L287) |

### `src/moment_set.cpp`

[Source](../../../src/moment_set.cpp) · [Responsibilities, contracts, and evidence](qmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `MomentSet::MomentSet` | [9](../../../src/moment_set.cpp#L9), [11](../../../src/moment_set.cpp#L11) |
| `MomentSet::size` | [13](../../../src/moment_set.cpp#L13) |
| `MomentSet::empty` | [18](../../../src/moment_set.cpp#L18) |
| `MomentSet::values` | [23](../../../src/moment_set.cpp#L23) |
| `MomentSet::operator[]` | [28](../../../src/moment_set.cpp#L28) |

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
| `validate_source_generalized_order` | [102](../../../src/moment_source.cpp#L102) |
| `validate_fixed_cqmom_moment_data` | [114](../../../src/moment_source.cpp#L114) |
| `validate_rate` | [121](../../../src/moment_source.cpp#L121) |
| `validate_source_rate` | [128](../../../src/moment_source.cpp#L128) |
| `validate_boundary_rate` | [140](../../../src/moment_source.cpp#L140) |
| `validate_cqmom_boundary_rate` | [152](../../../src/moment_source.cpp#L152) |
| `validate_cqmom_stage_state` | [159](../../../src/moment_source.cpp#L159) |
| `validate_coordinate_realizability_policies` | [171](../../../src/moment_source.cpp#L171) |
| `has_coordinate_realizability_policy` | [183](../../../src/moment_source.cpp#L183) |
| `extract_coordinate_marginal` | [190](../../../src/moment_source.cpp#L190) |
| `accepts_coordinate_realizability_policies` | [222](../../../src/moment_source.cpp#L222) |
| `boundary_callback_error` | [243](../../../src/moment_source.cpp#L243) |
| `cqmom_boundary_callback_error` | [248](../../../src/moment_source.cpp#L248) |
| `validate_request` | [253](../../../src/moment_source.cpp#L253), [284](../../../src/moment_source.cpp#L284) |
| `validate_cqmom_source_update_controls` | [294](../../../src/moment_source.cpp#L294) |
| `validate_moment_data_request` | [316](../../../src/moment_source.cpp#L316), [333](../../../src/moment_source.cpp#L333) |
| `validate_trusted_moment_data_request` | [324](../../../src/moment_source.cpp#L324), [343](../../../src/moment_source.cpp#L343) |
| `compose_source_and_boundary` | [353](../../../src/moment_source.cpp#L353) |
| `advance_moments_with_source` | [398](../../../src/moment_source.cpp#L398), [811](../../../src/moment_source.cpp#L811), [828](../../../src/moment_source.cpp#L828) |
| `compose_cqmom_source_and_boundary` | [450](../../../src/moment_source.cpp#L450) |
| `advance_validated_cqmom_moments_with_source` | [477](../../../src/moment_source.cpp#L477) |
| `QmomSourceContextAccess::create` | [495](../../../src/moment_source.cpp#L495) |
| `CqmomSourceContextAccess::create` | [541](../../../src/moment_source.cpp#L541) |
| `QmomSourceContext::QmomSourceContext` | [557](../../../src/moment_source.cpp#L557), [565](../../../src/moment_source.cpp#L565) |
| `QmomSourceContext::time` | [578](../../../src/moment_source.cpp#L578) |
| `QmomSourceContext::moment_count` | [583](../../../src/moment_source.cpp#L583) |
| `QmomSourceContext::size` | [588](../../../src/moment_source.cpp#L588) |
| `QmomSourceContext::values` | [593](../../../src/moment_source.cpp#L593) |
| `QmomSourceContext::value` | [598](../../../src/moment_source.cpp#L598) |
| `QmomSourceContext::moment` | [603](../../../src/moment_source.cpp#L603) |
| `QmomSourceContext::generalized_moment` | [620](../../../src/moment_source.cpp#L620) |
| `CqmomSourceContext::CqmomSourceContext` | [644](../../../src/moment_source.cpp#L644), [655](../../../src/moment_source.cpp#L655) |
| `CqmomSourceContext::time` | [669](../../../src/moment_source.cpp#L669) |
| `CqmomSourceContext::moment_count` | [674](../../../src/moment_source.cpp#L674) |
| `CqmomSourceContext::dimension_count` | [679](../../../src/moment_source.cpp#L679) |
| `CqmomSourceContext::size` | [684](../../../src/moment_source.cpp#L684) |
| `CqmomSourceContext::indices` | [689](../../../src/moment_source.cpp#L689) |
| `CqmomSourceContext::values` | [694](../../../src/moment_source.cpp#L694) |
| `CqmomSourceContext::index` | [699](../../../src/moment_source.cpp#L699) |
| `CqmomSourceContext::value` | [704](../../../src/moment_source.cpp#L704), [719](../../../src/moment_source.cpp#L719) |
| `CqmomSourceContext::position` | [709](../../../src/moment_source.cpp#L709) |
| `CqmomSourceContext::moment` | [729](../../../src/moment_source.cpp#L729) |
| `CqmomSourceContext::generalized_moment` | [745](../../../src/moment_source.cpp#L745) |
| `CqmomSourceContext::for_each_terminal_node` | [758](../../../src/moment_source.cpp#L758) |
| `CqmomSourceContext::reconstructed_stage_quadrature` | [765](../../../src/moment_source.cpp#L765) |
| `CqmomSourceContext::hierarchy` | [806](../../../src/moment_source.cpp#L806) |
| `advance_cqmom_moments_with_source` | [847](../../../src/moment_source.cpp#L847), [862](../../../src/moment_source.cpp#L862) |
| `advance_trusted_qmom_moments_with_source` | [881](../../../src/moment_source.cpp#L881), [906](../../../src/moment_source.cpp#L906) |
| `advance_trusted_cqmom_moments_with_source` | [933](../../../src/moment_source.cpp#L933), [948](../../../src/moment_source.cpp#L948), [967](../../../src/moment_source.cpp#L967), [985](../../../src/moment_source.cpp#L985) |

### `src/node_advection.cpp`

[Source](../../../src/node_advection.cpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `all_finite` | [24](../../../src/node_advection.cpp#L24) |
| `diagnostics_from_embedded_rk` | [29](../../../src/node_advection.cpp#L29) |
| `pack_nodes` | [43](../../../src/node_advection.cpp#L43) |
| `unpack_nodes` | [65](../../../src/node_advection.cpp#L65) |
| `validate_weights` | [79](../../../src/node_advection.cpp#L79) |
| `validate_abscissas` | [91](../../../src/node_advection.cpp#L91) |
| `node_rhs_from_velocity` | [118](../../../src/node_advection.cpp#L118) |
| `add_scaled` | [149](../../../src/node_advection.cpp#L149) |
| `integrate_nodes_embedded_rk` | [163](../../../src/node_advection.cpp#L163) |
| `integrate_nodes_ssp_rk3` | [200](../../../src/node_advection.cpp#L200) |
| `integrate_nodes` | [265](../../../src/node_advection.cpp#L265) |
| `scheme_name` | [296](../../../src/node_advection.cpp#L296) |
| `make_phase_space_velocity` | [395](../../../src/node_advection.cpp#L395) |
| `advance_nodes` | [312](../../../src/node_advection.cpp#L312), [340](../../../src/node_advection.cpp#L340), [351](../../../src/node_advection.cpp#L351) |
| `format_node_advection_options` | [439](../../../src/node_advection.cpp#L439) |

### `src/output.cpp`

[Source](../../../src/output.cpp) · [Responsibilities, contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `qmom_label` | [19](../../../src/output.cpp#L19) |
| `write_column_header` | [36](../../../src/output.cpp#L36) |
| `cqmom_label` | [52](../../../src/output.cpp#L52) |
| `sequential_indices` | [105](../../../src/output.cpp#L105) |
| `qmom_labels` | [121](../../../src/output.cpp#L121) |
| `cqmom_labels` | [137](../../../src/output.cpp#L137) |
| `flatten_selected_indices` | [153](../../../src/output.cpp#L153) |
| `all_cqmom_labels` | [164](../../../src/output.cpp#L164) |
| `cqmom_abscissa_label` | [169](../../../src/output.cpp#L169) |
| `qmom_quadrature_labels` | [174](../../../src/output.cpp#L174) |
| `cqmom_quadrature_labels` | [184](../../../src/output.cpp#L184) |
| `MomentWriter::for_qmom` | [208](../../../src/output.cpp#L208), [214](../../../src/output.cpp#L214), [220](../../../src/output.cpp#L220), [230](../../../src/output.cpp#L230) |
| `MomentWriter::for_cqmom` | [241](../../../src/output.cpp#L241), [248](../../../src/output.cpp#L248), [251](../../../src/output.cpp#L251), [254](../../../src/output.cpp#L254), [265](../../../src/output.cpp#L265), [267](../../../src/output.cpp#L267), [270](../../../src/output.cpp#L270) |
| `MomentWriter::MomentWriter` | [276](../../../src/output.cpp#L276), [294](../../../src/output.cpp#L294), [311](../../../src/output.cpp#L311), [328](../../../src/output.cpp#L328) |
| `MomentWriter::write_header` | [359](../../../src/output.cpp#L359) |
| `MomentWriter::write_row` | [364](../../../src/output.cpp#L364) |
| `MomentWriter::write_aligned_row` | [391](../../../src/output.cpp#L391) |
| `MomentWriter::write` | [422](../../../src/output.cpp#L422), [439](../../../src/output.cpp#L439), [462](../../../src/output.cpp#L462), [478](../../../src/output.cpp#L478), [500](../../../src/output.cpp#L500), [521](../../../src/output.cpp#L521), [548](../../../src/output.cpp#L548), [565](../../../src/output.cpp#L565), [578](../../../src/output.cpp#L578), [597](../../../src/output.cpp#L597) |
| `QuadratureWriter::for_qmom` | [611](../../../src/output.cpp#L611) |
| `QuadratureWriter::for_cqmom` | [616](../../../src/output.cpp#L616) |
| `QuadratureWriter::QuadratureWriter` | [613](../../../src/output.cpp#L613) |
| `QuadratureWriter::write_header` | [637](../../../src/output.cpp#L637) |
| `QuadratureWriter::write_row` | [642](../../../src/output.cpp#L642) |
| `QuadratureWriter::write_qmom_row` | [669](../../../src/output.cpp#L669) |
| `QuadratureWriter::write` | [699](../../../src/output.cpp#L699), [704](../../../src/output.cpp#L704), [709](../../../src/output.cpp#L709) |

### `src/phase_space_evolution.cpp`

[Source](../../../src/phase_space_evolution.cpp) · [Responsibilities, contracts, and evidence](phase_space_evolution.md)

| Function or overload family | Definition lines |
|---|---|
| `apply_batch` | [50](../../../src/phase_space_evolution.cpp#L50) |
| `validate_time_input` | [63](../../../src/phase_space_evolution.cpp#L63) |
| `validate_phase_space_evolution_scheme` | [75](../../../src/phase_space_evolution.cpp#L75) |
| `execute_phase_space_evolution_sequence` | [88](../../../src/phase_space_evolution.cpp#L88) |
| `validate_qmom_input` | [123](../../../src/phase_space_evolution.cpp#L123) |
| `validate_boundary` | [163](../../../src/phase_space_evolution.cpp#L163) |
| `validate_cqmom_boundary` | [171](../../../src/phase_space_evolution.cpp#L171) |
| `successful_reconstruction` | [179](../../../src/phase_space_evolution.cpp#L179) |
| `record_advection_diagnostics` | [184](../../../src/phase_space_evolution.cpp#L184) |
| `record_source_update_diagnostics` | [200](../../../src/phase_space_evolution.cpp#L200) |
| `qmom_from_quadrature` | [207](../../../src/phase_space_evolution.cpp#L207) |
| `add_advection_totals` | [225](../../../src/phase_space_evolution.cpp#L225) |
| `finalize_diagnostics` | [234](../../../src/phase_space_evolution.cpp#L234) |
| `record_failure_diagnostics` | [240](../../../src/phase_space_evolution.cpp#L240) |
| `qmom_failure` | [247](../../../src/phase_space_evolution.cpp#L247) |
| `select_post_source_qmom` | [255](../../../src/phase_space_evolution.cpp#L255) |
| `validate_qmom_advection_input` | [319](../../../src/phase_space_evolution.cpp#L319) |
| `validate_qmom_moment_source_input` | [332](../../../src/phase_space_evolution.cpp#L332) |
| `advance_qmom_advection_stage` | [368](../../../src/phase_space_evolution.cpp#L368) |
| `advance_qmom_moment_source_stage` | [418](../../../src/phase_space_evolution.cpp#L418) |
| `neutral_cqmom_message` | [481](../../../src/phase_space_evolution.cpp#L481) |
| `neutral_cqmom_diagnostics` | [494](../../../src/phase_space_evolution.cpp#L494) |
| `neutral_cqmom_result` | [503](../../../src/phase_space_evolution.cpp#L503) |
| `operator<` | [517](../../../src/phase_space_evolution.cpp#L517) |
| `active_branch_key` | [530](../../../src/phase_space_evolution.cpp#L530) |
| `path_has_prefix` | [535](../../../src/phase_space_evolution.cpp#L535) |
| `has_changed_strict_ancestor` | [551](../../../src/phase_space_evolution.cpp#L551) |
| `record_cqmom_node_count_transitions` | [569](../../../src/phase_space_evolution.cpp#L569) |
| `advance_qmom_phase_space_impl` | [659](../../../src/phase_space_evolution.cpp#L659) |
| `advance_qmom_advection` | [701](../../../src/phase_space_evolution.cpp#L701) |
| `advance_qmom_moment_source_impl` | [726](../../../src/phase_space_evolution.cpp#L726) |
| `advance_qmom_moment_source` | [752](../../../src/phase_space_evolution.cpp#L752), [758](../../../src/phase_space_evolution.cpp#L758) |
| `advance_qmom_phase_space` | [764](../../../src/phase_space_evolution.cpp#L764), [771](../../../src/phase_space_evolution.cpp#L771) |
| `initialize_cqmom_diagnostics` | [781](../../../src/phase_space_evolution.cpp#L781), [788](../../../src/phase_space_evolution.cpp#L788), [795](../../../src/phase_space_evolution.cpp#L795) |
| `validate_dense_cqmom_moments_for_phase_space_evolution` | [808](../../../src/phase_space_evolution.cpp#L808) |
| `validate_lean_cqmom_moments_for_phase_space_evolution` | [815](../../../src/phase_space_evolution.cpp#L815) |
| `validate_cqmom_state_for_phase_space_evolution` | [822](../../../src/phase_space_evolution.cpp#L822) |
| `validate_cqmom_advection_input` | [837](../../../src/phase_space_evolution.cpp#L837) |
| `validate_cqmom_moment_source_input` | [847](../../../src/phase_space_evolution.cpp#L847) |
| `validate_cqmom_phase_space_evolution_input` | [859](../../../src/phase_space_evolution.cpp#L859) |
| `CqmomStageOutcome::ok` | [886](../../../src/phase_space_evolution.cpp#L886) |
| `take_cqmom_advection_active_counts` | [900](../../../src/phase_space_evolution.cpp#L900) |
| `cqmom_state_failure` | [909](../../../src/phase_space_evolution.cpp#L909) |
| `cqmom_success` | [917](../../../src/phase_space_evolution.cpp#L917) |
| `materialize_cqmom_phase_space_evolution_result` | [922](../../../src/phase_space_evolution.cpp#L922) |
| `dense_state_from_accepted_quadrature` | [943](../../../src/phase_space_evolution.cpp#L943) |
| `lean_state_from_accepted_quadrature` | [959](../../../src/phase_space_evolution.cpp#L959) |
| `prepare_cqmom_phase_space_evolution_working_state` | [970](../../../src/phase_space_evolution.cpp#L970) |
| `advance_cqmom_advection_stage` | [1014](../../../src/phase_space_evolution.cpp#L1014) |
| `advance_cqmom_moment_source_stage` | [1111](../../../src/phase_space_evolution.cpp#L1111) |
| `advance_cqmom_advection_impl` | [1250](../../../src/phase_space_evolution.cpp#L1250) |
| `advance_cqmom_moment_source_impl` | [1257](../../../src/phase_space_evolution.cpp#L1257) |
| `advance_cqmom_phase_space_state_impl` | [1264](../../../src/phase_space_evolution.cpp#L1264) |
| `advance_cqmom_advection` | [1324](../../../src/phase_space_evolution.cpp#L1324) |
| `advance_cqmom_moment_source_public_impl` | [1346](../../../src/phase_space_evolution.cpp#L1346) |
| `advance_cqmom_moment_source` | [1369](../../../src/phase_space_evolution.cpp#L1369), [1375](../../../src/phase_space_evolution.cpp#L1375) |
| `advance_cqmom_phase_space` | [1382](../../../src/phase_space_evolution.cpp#L1382), [1389](../../../src/phase_space_evolution.cpp#L1389) |
| `advance_qmom_advection_batch` | [1397](../../../src/phase_space_evolution.cpp#L1397) |
| `advance_cqmom_advection_batch` | [1405](../../../src/phase_space_evolution.cpp#L1405) |
| `advance_qmom_moment_source_batch` | [1413](../../../src/phase_space_evolution.cpp#L1413), [1429](../../../src/phase_space_evolution.cpp#L1429) |
| `advance_cqmom_moment_source_batch` | [1421](../../../src/phase_space_evolution.cpp#L1421), [1438](../../../src/phase_space_evolution.cpp#L1438) |
| `advance_qmom_phase_space_batch` | [1447](../../../src/phase_space_evolution.cpp#L1447), [1464](../../../src/phase_space_evolution.cpp#L1464) |
| `advance_cqmom_phase_space_batch` | [1455](../../../src/phase_space_evolution.cpp#L1455), [1473](../../../src/phase_space_evolution.cpp#L1473) |
| `advance_qmom_advection_batch_with_index` | [1483](../../../src/phase_space_evolution.cpp#L1483) |
| `advance_cqmom_advection_batch_with_index` | [1500](../../../src/phase_space_evolution.cpp#L1500) |
| `advance_qmom_moment_source_batch_with_index` | [1517](../../../src/phase_space_evolution.cpp#L1517), [1547](../../../src/phase_space_evolution.cpp#L1547) |
| `advance_cqmom_moment_source_batch_with_index` | [1532](../../../src/phase_space_evolution.cpp#L1532), [1568](../../../src/phase_space_evolution.cpp#L1568) |
| `advance_qmom_phase_space_batch_with_index` | [1589](../../../src/phase_space_evolution.cpp#L1589), [1634](../../../src/phase_space_evolution.cpp#L1634) |
| `advance_cqmom_phase_space_batch_with_index` | [1611](../../../src/phase_space_evolution.cpp#L1611), [1662](../../../src/phase_space_evolution.cpp#L1662) |

### `src/public_api.cpp`

[Source](../../../src/public_api.cpp) · [Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `to_wheeler_options` | [28](../../../src/public_api.cpp#L28) |
| `compute_qmom_unscaled` | [35](../../../src/public_api.cpp#L35) |
| `compute_qmom_with_automatic_scaling` | [52](../../../src/public_api.cpp#L52) |
| `to_qmom_status` | [74](../../../src/public_api.cpp#L74) |
| `hierarchy` | [96](../../../src/public_api.cpp#L96) |
| `qmom_options` | [103](../../../src/public_api.cpp#L103) |
| `conditional_options` | [110](../../../src/public_api.cpp#L110) |
| `adaptive_branch_counts` | [117](../../../src/public_api.cpp#L117) |
| `diagnostic_level` | [124](../../../src/public_api.cpp#L124) |
| `vandermonde_solver` | [131](../../../src/public_api.cpp#L131) |
| `coordinate_support` | [138](../../../src/public_api.cpp#L138) |
| `moment_storage` | [145](../../../src/public_api.cpp#L145) |
| `diagnostic_level_name` | [150](../../../src/public_api.cpp#L150) |
| `vandermonde_solver_name` | [164](../../../src/public_api.cpp#L164) |
| `compute_qmom` | [185](../../../src/public_api.cpp#L185), [190](../../../src/public_api.cpp#L190), [196](../../../src/public_api.cpp#L196) |
| `QmomState::QmomState` | [201](../../../src/public_api.cpp#L201) |
| `QmomCase::initial_state` | [208](../../../src/public_api.cpp#L208) |
| `to_string` | [232](../../../src/public_api.cpp#L232) |
| `format_cqmom_settings` | [244](../../../src/public_api.cpp#L244) |
| `CqmomCase::CqmomCase` | [256](../../../src/public_api.cpp#L256), [258](../../../src/public_api.cpp#L258), [260](../../../src/public_api.cpp#L260), [262](../../../src/public_api.cpp#L262), [273](../../../src/public_api.cpp#L273) |
| `CqmomCase::operator=` | [264](../../../src/public_api.cpp#L264), [275](../../../src/public_api.cpp#L275) |
| `CqmomCase::hierarchy` | [279](../../../src/public_api.cpp#L279) |
| `CqmomCase::qmom_options` | [284](../../../src/public_api.cpp#L284) |
| `CqmomCase::conditional_options` | [289](../../../src/public_api.cpp#L289) |
| `CqmomCase::adaptive_branch_counts` | [294](../../../src/public_api.cpp#L294) |
| `CqmomCase::diagnostic_level` | [299](../../../src/public_api.cpp#L299) |
| `CqmomCase::vandermonde_solver` | [304](../../../src/public_api.cpp#L304) |
| `CqmomCase::coordinate_support` | [309](../../../src/public_api.cpp#L309) |
| `CqmomCase::moment_storage` | [314](../../../src/public_api.cpp#L314) |
| `CqmomCase::initial_state` | [319](../../../src/public_api.cpp#L319) |
| `CqmomState::CqmomState` | [364](../../../src/public_api.cpp#L364), [366](../../../src/public_api.cpp#L366), [377](../../../src/public_api.cpp#L377) |
| `CqmomState::operator=` | [368](../../../src/public_api.cpp#L368), [379](../../../src/public_api.cpp#L379) |
| `CqmomState::hierarchy` | [383](../../../src/public_api.cpp#L383) |
| `CqmomState::quadrature` | [388](../../../src/public_api.cpp#L388) |
| `CqmomState::qmom_options` | [397](../../../src/public_api.cpp#L397) |
| `CqmomState::conditional_options` | [402](../../../src/public_api.cpp#L402) |
| `CqmomState::adaptive_branch_counts` | [407](../../../src/public_api.cpp#L407) |
| `CqmomState::diagnostic_level` | [412](../../../src/public_api.cpp#L412) |
| `CqmomState::vandermonde_solver` | [417](../../../src/public_api.cpp#L417) |
| `CqmomState::coordinate_support` | [422](../../../src/public_api.cpp#L422) |
| `CqmomState::moment_storage` | [427](../../../src/public_api.cpp#L427) |
| `CqmomState::active_branch_counts` | [432](../../../src/public_api.cpp#L432) |
| `CqmomState::moment_indices` | [437](../../../src/public_api.cpp#L437) |
| `CqmomState::moment_values` | [444](../../../src/public_api.cpp#L444) |
| `CqmomState::moment` | [451](../../../src/public_api.cpp#L451) |
| `validate_prepared_fixed_metadata` | [483](../../../src/public_api.cpp#L483), [495](../../../src/public_api.cpp#L495) |
| `make_prepared_fixed_metadata` | [503](../../../src/public_api.cpp#L503), [522](../../../src/public_api.cpp#L522) |
| `CqmomConstructionAccess::create_case` | [542](../../../src/public_api.cpp#L542), [548](../../../src/public_api.cpp#L548) |
| `CqmomConstructionAccess::create_state` | [554](../../../src/public_api.cpp#L554), [556](../../../src/public_api.cpp#L556), [562](../../../src/public_api.cpp#L562), [574](../../../src/public_api.cpp#L574), [577](../../../src/public_api.cpp#L577), [583](../../../src/public_api.cpp#L583), [589](../../../src/public_api.cpp#L589), [601](../../../src/public_api.cpp#L601) |
| `CqmomConstructionAccess::case_dense_case` | [608](../../../src/public_api.cpp#L608) |
| `CqmomConstructionAccess::case_lean_case` | [613](../../../src/public_api.cpp#L613) |
| `CqmomConstructionAccess::state_dense_case` | [618](../../../src/public_api.cpp#L618) |
| `CqmomConstructionAccess::state_lean_case` | [623](../../../src/public_api.cpp#L623) |
| `CqmomConstructionAccess::state_active_branch_counts` | [628](../../../src/public_api.cpp#L628) |
| `CqmomConstructionAccess::state_accepted_quadrature` | [633](../../../src/public_api.cpp#L633) |
| `CqmomConstructionAccess::case_prepared_fixed_metadata` | [638](../../../src/public_api.cpp#L638) |
| `CqmomConstructionAccess::state_prepared_fixed_metadata` | [643](../../../src/public_api.cpp#L643) |
| `CqmomConstructionAccess::create_quadrature` | [648](../../../src/public_api.cpp#L648) |
| `CqmomQuadrature::CqmomQuadrature` | [658](../../../src/public_api.cpp#L658) |
| `CqmomQuadrature::hierarchy` | [385](../../../src/public_api.cpp#L385) |
| `CqmomQuadrature::weight` | [665](../../../src/public_api.cpp#L665) |
| `CqmomQuadrature::abscissa` | [674](../../../src/public_api.cpp#L674) |
| `CqmomQuadrature::node` | [687](../../../src/public_api.cpp#L687) |
| `to_cqmom_status` | [700](../../../src/public_api.cpp#L700) |
| `to_cqmom_diagnostics` | [758](../../../src/public_api.cpp#L758), [814](../../../src/public_api.cpp#L814) |
| `terminal_paths_for_public_quadrature` | [813](../../../src/public_api.cpp#L813) |
| `failed_phase_space_dimension` | [815](../../../src/public_api.cpp#L815) |
| `to_cqmom_result` | [859](../../../src/public_api.cpp#L859) |
| `compute_cqmom` | [828](../../../src/public_api.cpp#L828), [833](../../../src/public_api.cpp#L833), [980](../../../src/public_api.cpp#L980) |
| `compute_cqmom_internal` | [942](../../../src/public_api.cpp#L942), [985](../../../src/public_api.cpp#L985) |

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
| `hankel_matrix` | [229](../../../src/realizability.cpp#L229) |
| `interval_quadratic_localizing_matrix` | [244](../../../src/realizability.cpp#L244) |
| `interval_linear_localizing_matrix` | [261](../../../src/realizability.cpp#L261) |
| `extension_column` | [279](../../../src/realizability.cpp#L279) |
| `projection_on_eigenvector` | [289](../../../src/realizability.cpp#L289) |
| `check_range_consistency` | [299](../../../src/realizability.cpp#L299) |
| `incorporate_range_consistency` | [327](../../../src/realizability.cpp#L327) |
| `finish_classification` | [358](../../../src/realizability.cpp#L358) |
| `validate_common_inputs` | [374](../../../src/realizability.cpp#L374) |
| `check_hamburger` | [409](../../../src/realizability.cpp#L409) |
| `check_stieltjes` | [454](../../../src/realizability.cpp#L454) |
| `check_hausdorff` | [511](../../../src/realizability.cpp#L511) |
| `with_support_context` | [566](../../../src/realizability.cpp#L566) |
| `check_realizability` | [577](../../../src/realizability.cpp#L577), [600](../../../src/realizability.cpp#L600) |

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
| `IndexedMomentSet::IndexedMomentSet` | [33](../../../src/recursive_quadrature.cpp#L33) |
| `IndexedMomentSet::key_extents` | [70](../../../src/recursive_quadrature.cpp#L70) |
| `IndexedMomentSet::indices` | [75](../../../src/recursive_quadrature.cpp#L75) |
| `IndexedMomentSet::values` | [80](../../../src/recursive_quadrature.cpp#L80) |
| `IndexedMomentSet::mutable_values` | [85](../../../src/recursive_quadrature.cpp#L85) |
| `IndexedMomentSet::size` | [90](../../../src/recursive_quadrature.cpp#L90) |
| `IndexedMomentSet::contains` | [95](../../../src/recursive_quadrature.cpp#L95) |
| `IndexedMomentSet::offset` | [112](../../../src/recursive_quadrature.cpp#L112) |
| `IndexedMomentSet::moment` | [123](../../../src/recursive_quadrature.cpp#L123) |
| `CqmomMomentView::CqmomMomentView` | [128](../../../src/recursive_quadrature.cpp#L128), [130](../../../src/recursive_quadrature.cpp#L130) |
| `CqmomMomentView::moment` | [123](../../../src/recursive_quadrature.cpp#L123) |
| `create_indexed_cqmom_moment_set` | [141](../../../src/recursive_quadrature.cpp#L141) |
| `create_ordered_root_layout` | [149](../../../src/recursive_quadrature.cpp#L149) |
| `append_unique_index` | [176](../../../src/recursive_quadrature.cpp#L176) |
| `local_branch_offset` | [184](../../../src/recursive_quadrature.cpp#L184) |
| `create_child_layout_case` | [194](../../../src/recursive_quadrature.cpp#L194) |
| `create_branch_moment_layout` | [234](../../../src/recursive_quadrature.cpp#L234) |
| `PreparedFixedCqmomMetadata::PreparedFixedCqmomMetadata` | [287](../../../src/recursive_quadrature.cpp#L287) |
| `create_prepared_fixed_cqmom_metadata` | [297](../../../src/recursive_quadrature.cpp#L297) |
| `prepared_fixed_cqmom_metadata_is_compatible` | [307](../../../src/recursive_quadrature.cpp#L307) |
| `extract_wheeler_moments` | [322](../../../src/recursive_quadrature.cpp#L322) |
| `child_layout_for_parent_node_count` | [341](../../../src/recursive_quadrature.cpp#L341) |
| `make_branch_rhs_matrix` | [355](../../../src/recursive_quadrature.cpp#L355) |
| `validate_positive_extents` | [384](../../../src/recursive_quadrature.cpp#L384) |
| `checked_multiply` | [396](../../../src/recursive_quadrature.cpp#L396) |
| `default_conditioning_order` | [406](../../../src/recursive_quadrature.cpp#L406) |
| `validate_conditioning_order` | [413](../../../src/recursive_quadrature.cpp#L413) |
| `ordered_extents` | [435](../../../src/recursive_quadrature.cpp#L435) |
| `from_ordered_index` | [446](../../../src/recursive_quadrature.cpp#L446) |
| `validate_abscissas` | [456](../../../src/recursive_quadrature.cpp#L456) |
| `valid_conditional_options` | [497](../../../src/recursive_quadrature.cpp#L497) |
| `valid_diagnostic_level` | [504](../../../src/recursive_quadrature.cpp#L504) |
| `records_success_branch_diagnostics` | [516](../../../src/recursive_quadrature.cpp#L516) |
| `create_ordered_root_slice` | [521](../../../src/recursive_quadrature.cpp#L521) |
| `one_dimensional_moments` | [536](../../../src/recursive_quadrature.cpp#L536) |
| `reconstruct_child_slices` | [541](../../../src/recursive_quadrature.cpp#L541) |
| `create_recursive_result` | [633](../../../src/recursive_quadrature.cpp#L633) |
| `record_failure` | [638](../../../src/recursive_quadrature.cpp#L638) |
| `validate_recursive_request_preflight` | [647](../../../src/recursive_quadrature.cpp#L647) |
| `resolve_cqmom_wheeler_options` | [671](../../../src/recursive_quadrature.cpp#L671) |
| `operator<` | [687](../../../src/recursive_quadrature.cpp#L687) |
| `path_has_prefix` | [712](../../../src/recursive_quadrature.cpp#L712) |
| `branch_limit_invalidated` | [728](../../../src/recursive_quadrature.cpp#L728) |
| `validate_active_branch_topology` | [740](../../../src/recursive_quadrature.cpp#L740) |
| `recursive_branch_node_count_limit` | [812](../../../src/recursive_quadrature.cpp#L812) |
| `compute_branch_quadrature` | [830](../../../src/recursive_quadrature.cpp#L830) |
| `validate_coordinate_support` | [835](../../../src/recursive_quadrature.cpp#L835) |
| `support_interval_for_dimension` | [862](../../../src/recursive_quadrature.cpp#L862) |
| `record_bounded_candidate_diagnostics` | [872](../../../src/recursive_quadrature.cpp#L872) |
| `branch_growth_trial_allowed` | [890](../../../src/recursive_quadrature.cpp#L890) |
| `product_extents` | [904](../../../src/recursive_quadrature.cpp#L904) |
| `flatten_multi_index` | [916](../../../src/recursive_quadrature.cpp#L916) |
| `unflatten_multi_index` | [936](../../../src/recursive_quadrature.cpp#L936) |
| `enumerate_multi_indices` | [1086](../../../src/recursive_quadrature.cpp#L1086) |
| `PhaseSpaceHierarchy::PhaseSpaceHierarchy` | [1098](../../../src/recursive_quadrature.cpp#L1098), [1105](../../../src/recursive_quadrature.cpp#L1105) |
| `PhaseSpaceHierarchy::number_of_dimensions` | [1133](../../../src/recursive_quadrature.cpp#L1133) |
| `PhaseSpaceHierarchy::node_count` | [1111](../../../src/recursive_quadrature.cpp#L1111) |
| `PhaseSpaceHierarchy::total_nodes` | [1342](../../../src/recursive_quadrature.cpp#L1342) |
| `PhaseSpaceHierarchy::moment_extents` | [1120](../../../src/recursive_quadrature.cpp#L1120) |
| `PhaseSpaceHierarchy::flatten_node_index` | [1131](../../../src/recursive_quadrature.cpp#L1131) |
| `PhaseSpaceHierarchy::unflatten_node_index` | [1140](../../../src/recursive_quadrature.cpp#L1140) |
| `create_cqmom_moment_layout` | [1145](../../../src/recursive_quadrature.cpp#L1145) |
| `CqmomMomentData::CqmomMomentData` | [1189](../../../src/recursive_quadrature.cpp#L1189) |
| `CqmomMomentData::contains_required_index` | [1194](../../../src/recursive_quadrature.cpp#L1194) |
| `CqmomMomentData::moment` | [1199](../../../src/recursive_quadrature.cpp#L1199) |
| `extract_cqmom_moments` | [1215](../../../src/recursive_quadrature.cpp#L1215) |
| `create_cqmom_moment_data` | [1242](../../../src/recursive_quadrature.cpp#L1242) |
| `enumerate_quadrature_node_multi_indices` | [1264](../../../src/recursive_quadrature.cpp#L1264) |
| `rectangular_terminal_paths` | [1279](../../../src/recursive_quadrature.cpp#L1279) |
| `RecursiveQuadrature::RecursiveQuadrature` | [1314](../../../src/recursive_quadrature.cpp#L1314), [1321](../../../src/recursive_quadrature.cpp#L1321) |
| `RecursiveQuadrature::hierarchy` | [1332](../../../src/recursive_quadrature.cpp#L1332) |
| `RecursiveQuadrature::number_of_dimensions` | [1133](../../../src/recursive_quadrature.cpp#L1133) |
| `RecursiveQuadrature::total_nodes` | [1342](../../../src/recursive_quadrature.cpp#L1342) |
| `RecursiveQuadrature::weights` | [1347](../../../src/recursive_quadrature.cpp#L1347) |
| `RecursiveQuadrature::abscissas` | [1352](../../../src/recursive_quadrature.cpp#L1352) |
| `RecursiveQuadrature::adaptive_branch_counts_enabled` | [1357](../../../src/recursive_quadrature.cpp#L1357) |
| `RecursiveQuadrature::has_rectangular_storage` | [1362](../../../src/recursive_quadrature.cpp#L1362) |
| `RecursiveQuadrature::tree_nodes` | [1367](../../../src/recursive_quadrature.cpp#L1367) |
| `RecursiveQuadrature::terminal_node_ids` | [1372](../../../src/recursive_quadrature.cpp#L1372) |
| `RecursiveQuadrature::terminal_paths` | [1377](../../../src/recursive_quadrature.cpp#L1377) |
| `RecursiveQuadrature::weight` | [1382](../../../src/recursive_quadrature.cpp#L1382) |
| `RecursiveQuadrature::abscissa` | [1391](../../../src/recursive_quadrature.cpp#L1391) |
| `RecursiveQuadrature::node` | [1404](../../../src/recursive_quadrature.cpp#L1404) |
| `RecursiveQuadrature::set_weight` | [1414](../../../src/recursive_quadrature.cpp#L1414) |
| `RecursiveQuadrature::set_abscissa` | [1423](../../../src/recursive_quadrature.cpp#L1423) |
| `RecursiveQuadrature::set_adaptive_branch_counts_enabled` | [1436](../../../src/recursive_quadrature.cpp#L1436) |
| `RecursiveQuadrature::add_tree_node` | [1449](../../../src/recursive_quadrature.cpp#L1449) |
| `RecursiveQuadrature::set_terminal_quadrature` | [1461](../../../src/recursive_quadrature.cpp#L1461) |
| `compute_recursive_quadrature_from_moment_view` | [1487](../../../src/recursive_quadrature.cpp#L1487) |
| `compute_recursive_quadrature_from_validated_moment_view` | [1812](../../../src/recursive_quadrature.cpp#L1812) |
| `compute_recursive_quadrature_from_moment_data` | [1821](../../../src/recursive_quadrature.cpp#L1821) |
| `compute_recursive_quadrature_from_trusted_moment_data` | [1856](../../../src/recursive_quadrature.cpp#L1856) |
| `compute_recursive_quadrature` | [1889](../../../src/recursive_quadrature.cpp#L1889) |

### `src/runge_kutta.cpp`

[Source](../../../src/runge_kutta.cpp) · [Responsibilities, contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `require_same_size` | [22](../../../src/runge_kutta.cpp#L22) |
| `euclidean_norm` | [30](../../../src/runge_kutta.cpp#L30) |
| `add_scaled` | [40](../../../src/runge_kutta.cpp#L40) |
| `sign` | [51](../../../src/runge_kutta.cpp#L51) |
| `validate_or_default_max_step` | [56](../../../src/runge_kutta.cpp#L56) |
| `minimum_step` | [65](../../../src/runge_kutta.cpp#L65) |
| `error_norm` | [81](../../../src/runge_kutta.cpp#L81) |
| `invalid_embedded_rk_result` | [101](../../../src/runge_kutta.cpp#L101) |
| `stage_accepted` | [106](../../../src/runge_kutta.cpp#L106) |
| `try_ssp_rk3_step` | [111](../../../src/runge_kutta.cpp#L111) |
| `ssp_rk3_step` | [173](../../../src/runge_kutta.cpp#L173) |
| `integrate_embedded_rk` | [219](../../../src/runge_kutta.cpp#L219) |

### `src/wheeler.cpp`

[Source](../../../src/wheeler.cpp) · [Responsibilities, contracts, and evidence](qmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `one_node_result` | [42](../../../src/wheeler.cpp#L42) |
| `negative_mass_result` | [55](../../../src/wheeler.cpp#L55) |
| `compute_recurrence` | [63](../../../src/wheeler.cpp#L63) |
| `make_jacobi_matrix` | [93](../../../src/wheeler.cpp#L93) |
| `symmetric_jacobi_eigen` | [120](../../../src/wheeler.cpp#L120) |
| `recurrence_is_finite` | [123](../../../src/wheeler.cpp#L123) |
| `has_negative_recurrence_b` | [128](../../../src/wheeler.cpp#L128) |
| `adaptive_sigma_tolerance` | [133](../../../src/wheeler.cpp#L133) |
| `compute_weights` | [144](../../../src/wheeler.cpp#L144) |
| `adaptive_accepts` | [156](../../../src/wheeler.cpp#L156) |
| `required_moments_are_finite` | [209](../../../src/wheeler.cpp#L209) |
| `attempt_quadrature_count` | [221](../../../src/wheeler.cpp#L221) |
| `compute_quadrature` | [309](../../../src/wheeler.cpp#L309) |
| `select_bounded_quadrature` | [382](../../../src/wheeler.cpp#L382) |

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
