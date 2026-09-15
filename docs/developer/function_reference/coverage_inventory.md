# Definition coverage and navigation inventory

This inventory routes core implementation definitions to functional
explanations. It complements the algorithm sections and responsibility tables;
a symbol listing does not by itself establish intent or test coverage.

## Scope and method

The inventory covers the 23 core `.cpp` files and installed/private headers. It
lists explicit functions, methods, constructors, and operators by defining
source range. Defaulted or implicit special members are omitted; consequential
lambdas are documented with their owning function. Instrumentation definitions
are included for its enabled preprocessor branch, and debug-only contract
checks for theirs. Line links are navigation aids, not permanent API
identifiers: they drift whenever a file is edited above a listed definition, so
treat the file and the name as authoritative and the line as a hint.

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

A clean run prints nothing. When it reports drift, the nearest definition of
the same name in the same file is almost always the new location; when it
reports a name the file no longer contains at all, the definition moved or was
removed and the row needs a decision rather than a new number.

That check only confirms that existing rows point somewhere valid. It cannot
see a definition that has no row at all, which is how a whole type can be added
without the ledger noticing — the two flat index-table headers were missed
exactly that way and had no section until they were audited.

The rule that closes it: **when a change adds a file under `src/` or
`include/qmomkit/`, add its section here in the same change**, and when it adds
a definition to a file that already has a section, add its row. A section that
lists a file is a claim that it lists all of that file's definitions.

The inventory includes **1015 definition locations**, grouped by named function
within each defining file. Overloads and separately compiled instrumentation
branches share a row with distinct line links. The count includes trivial
helpers and is not a count of numerical algorithms. Detailed pages cover the
important algorithm families; entry, accessor, diagnostic, parsing, formatting,
and forwarding functions are grouped by responsibility. No coverage percentage
or exhaustive branch-coverage claim is made.

The core numerical backbone receives Level A treatment: Wheeler and scaling;
support-specific realizability; weighted conditional systems and numerical
kernels; required/conditional layouts and recursive branch policy;
index/topology transformations; raw/generalized moments and rates; analytical
initialization; stage-local source closure; SSP/embedded RK; and accepted
evolution transitions. Configuration policy resolution and public workflows
combine Levels A/B; ordinary parsing, ownership, getters, output, and counters
receive lighter B/C coverage. Specific tests and remaining gaps are attached to
the functional explanations, not inferred from proximity of filenames.

Tests, examples, configuration fixtures, package metadata, and build files are
linked where they clarify a contract, but this inventory is limited to the core
library. The optional `tools/monte_carlo/` companion, plotting scripts,
benchmark harnesses, and application physics in examples remain separate
consumer/package responsibilities. The inventory does not claim to document
every function in those consumers.

## Definition ledger

Each file heading links to its definitions and owning functional reference.
Local scope means an anonymous namespace; public/internal ownership and class
contracts are explained on the destination page. For out-of-class methods, the
owner is retained when present in the definition signature. Trivial defaulted
ownership operations are covered by the state/result helper tables, even when
excluded from the executable-body count.

### `include/qmomkit/conditional_moments.hpp`

[Source](../../../include/qmomkit/conditional_moments.hpp) · [Responsibilities,
contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `ConditionalMomentResult::ok` | [70](../../../include/qmomkit/conditional_moments.hpp#L70) |

### `include/qmomkit/conditioning_path.hpp`

[Source](../../../include/qmomkit/conditioning_path.hpp) · [Responsibilities,
contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `ConditioningPathView::ConditioningPathView` | [60](../../../include/qmomkit/conditioning_path.hpp#L60) |
| `ConditioningPathView::empty` | [65](../../../include/qmomkit/conditioning_path.hpp#L65) |
| `ConditioningPathView::size` | [69](../../../include/qmomkit/conditioning_path.hpp#L69) |
| `ConditioningPathView::begin` | [75](../../../include/qmomkit/conditioning_path.hpp#L75) |
| `ConditioningPathView::end` | [79](../../../include/qmomkit/conditioning_path.hpp#L79) |
| `ConditioningPathView::data` | [83](../../../include/qmomkit/conditioning_path.hpp#L83) |
| `ConditioningPathView::to_conditioning_path` | [89](../../../include/qmomkit/conditioning_path.hpp#L89) |
| `ConditioningPathTable::const_iterator::operator*` | [132](../../../include/qmomkit/conditioning_path.hpp#L132) |
| `ConditioningPathTable::const_iterator::operator++` | [136](../../../include/qmomkit/conditioning_path.hpp#L136), [141](../../../include/qmomkit/conditioning_path.hpp#L141) |
| `ConditioningPathTable::const_iterator::operator==` | [147](../../../include/qmomkit/conditioning_path.hpp#L147) |
| `ConditioningPathTable::const_iterator::operator!=` | [151](../../../include/qmomkit/conditioning_path.hpp#L151) |
| `ConditioningPathTable::size` | [172](../../../include/qmomkit/conditioning_path.hpp#L172) |
| `ConditioningPathTable::levels` | [177](../../../include/qmomkit/conditioning_path.hpp#L177) |
| `ConditioningPathTable::empty` | [181](../../../include/qmomkit/conditioning_path.hpp#L181) |
| `ConditioningPathTable::operator[]` | [185](../../../include/qmomkit/conditioning_path.hpp#L185) |
| `ConditioningPathTable::begin` | [191](../../../include/qmomkit/conditioning_path.hpp#L191) |
| `ConditioningPathTable::end` | [195](../../../include/qmomkit/conditioning_path.hpp#L195) |
| `ConditioningPathTable::node_indices` | [201](../../../include/qmomkit/conditioning_path.hpp#L201) |

### `include/qmomkit/configuration.hpp`

[Source](../../../include/qmomkit/configuration.hpp) · [Responsibilities,
contracts, and evidence](initialization_and_configuration.md)

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

[Source](../../../include/qmomkit/cqmom.hpp) · [Responsibilities, contracts,
and evidence](state_and_representations.md)

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
| `CqmomQuadrature::hierarchy` | [44](../../../include/qmomkit/cqmom.hpp#L44) |

### `include/qmomkit/distributions.hpp`

[Source](../../../include/qmomkit/distributions.hpp) · [Responsibilities,
contracts, and evidence](initialization_and_configuration.md)

| Function or overload family | Definition lines |
|---|---|
| `RawMomentResult::ok` | [118](../../../include/qmomkit/distributions.hpp#L118) |
| `SelectedRawMomentResult::ok` | [154](../../../include/qmomkit/distributions.hpp#L154) |
| `TensorMomentResult::ok` | [168](../../../include/qmomkit/distributions.hpp#L168) |

### `include/qmomkit/moment_source.hpp`

[Source](../../../include/qmomkit/moment_source.hpp) · [Responsibilities,
contracts, and evidence](moment_sources.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomMomentSourceUpdateResult::ok` | [261](../../../include/qmomkit/moment_source.hpp#L261) |
| `CqmomMomentSourceUpdateResult::ok` | [261](../../../include/qmomkit/moment_source.hpp#L261) |

### `include/qmomkit/multi_index.hpp`

[Source](../../../include/qmomkit/multi_index.hpp) · [Responsibilities,
contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `MultiIndexView::size` | [45](../../../include/qmomkit/multi_index.hpp#L45) |
| `MultiIndexView::empty` | [49](../../../include/qmomkit/multi_index.hpp#L49) |
| `MultiIndexView::operator[]` | [53](../../../include/qmomkit/multi_index.hpp#L53) |
| `MultiIndexView::begin` | [57](../../../include/qmomkit/multi_index.hpp#L57) |
| `MultiIndexView::end` | [61](../../../include/qmomkit/multi_index.hpp#L61) |
| `MultiIndexView::data` | [65](../../../include/qmomkit/multi_index.hpp#L65) |
| `MultiIndexView::to_multi_index` | [71](../../../include/qmomkit/multi_index.hpp#L71) |
| `MultiIndexTable::const_iterator::operator*` | [127](../../../include/qmomkit/multi_index.hpp#L127) |
| `MultiIndexTable::const_iterator::operator++` | [131](../../../include/qmomkit/multi_index.hpp#L131), [136](../../../include/qmomkit/multi_index.hpp#L136) |
| `MultiIndexTable::const_iterator::operator==` | [142](../../../include/qmomkit/multi_index.hpp#L142) |
| `MultiIndexTable::const_iterator::operator!=` | [146](../../../include/qmomkit/multi_index.hpp#L146) |
| `MultiIndexTable::const_iterator::operator-` | [150](../../../include/qmomkit/multi_index.hpp#L150) |
| `MultiIndexTable::size` | [168](../../../include/qmomkit/multi_index.hpp#L168) |
| `MultiIndexTable::rank` | [172](../../../include/qmomkit/multi_index.hpp#L172) |
| `MultiIndexTable::empty` | [176](../../../include/qmomkit/multi_index.hpp#L176) |
| `MultiIndexTable::operator[]` | [180](../../../include/qmomkit/multi_index.hpp#L180) |
| `MultiIndexTable::begin` | [186](../../../include/qmomkit/multi_index.hpp#L186) |
| `MultiIndexTable::end` | [190](../../../include/qmomkit/multi_index.hpp#L190) |
| `MultiIndexTable::values` | [196](../../../include/qmomkit/multi_index.hpp#L196) |

### `include/qmomkit/node_advection.hpp`

[Source](../../../include/qmomkit/node_advection.hpp) · [Responsibilities,
contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomNodeAdvectionResult::ok` | [101](../../../include/qmomkit/node_advection.hpp#L101) |
| `CqmomNodeAdvectionResult::ok` | [101](../../../include/qmomkit/node_advection.hpp#L101) |

### `include/qmomkit/phase_space_evolution.hpp`

[Source](../../../include/qmomkit/phase_space_evolution.hpp) ·
[Responsibilities, contracts, and evidence](phase_space_evolution.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomPhaseSpaceEvolutionResult::ok` | [143](../../../include/qmomkit/phase_space_evolution.hpp#L143) |
| `CqmomPhaseSpaceEvolutionResult::ok` | [158](../../../include/qmomkit/phase_space_evolution.hpp#L158) |

### `include/qmomkit/phase_space_hierarchy.hpp`

[Source](../../../include/qmomkit/phase_space_hierarchy.hpp) ·
[Responsibilities, contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `PhaseSpaceHierarchy::node_counts` | [32](../../../include/qmomkit/phase_space_hierarchy.hpp#L32) |
| `PhaseSpaceHierarchy::conditioning_order` | [36](../../../include/qmomkit/phase_space_hierarchy.hpp#L36) |
| `PhaseSpaceHierarchy::number_of_dimensions` | [28](../../../include/qmomkit/phase_space_hierarchy.hpp#L28) |
| `PhaseSpaceHierarchy::total_nodes` | [40](../../../include/qmomkit/phase_space_hierarchy.hpp#L40) |

### `include/qmomkit/qmom.hpp`

[Source](../../../include/qmomkit/qmom.hpp) · [Responsibilities, contracts, and
evidence](state_and_representations.md)

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

[Source](../../../include/qmomkit/quadrature.hpp) · [Responsibilities,
contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `Quadrature1D::ok` | [30](../../../include/qmomkit/quadrature.hpp#L30) |

### `include/qmomkit/realizability_policy.hpp`

[Source](../../../include/qmomkit/realizability_policy.hpp) ·
[Responsibilities, contracts, and evidence](realizability.md)

| Function or overload family | Definition lines |
|---|---|
| `accepts_realizability_policy` | [56](../../../include/qmomkit/realizability_policy.hpp#L56) |
| `make_realizability_acceptance_predicate` | [95](../../../include/qmomkit/realizability_policy.hpp#L95) |
| `make_realizability_stage_acceptance_predicate` | [113](../../../include/qmomkit/realizability_policy.hpp#L113) |

### `include/qmomkit/runge_kutta.hpp`

[Source](../../../include/qmomkit/runge_kutta.hpp) · [Responsibilities,
contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `SspRk3Result::ok` | [84](../../../include/qmomkit/runge_kutta.hpp#L84) |

### `src/conditional_moments.cpp`

[Source](../../../src/conditional_moments.cpp) · [Responsibilities, contracts,
and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `invalid_result` | [25](../../../src/conditional_moments.cpp#L25) |
| `valid_options` | [40](../../../src/conditional_moments.cpp#L40) |
| `extract_rhs` | [47](../../../src/conditional_moments.cpp#L47) |
| `ConditionalMomentResult::moment` | [72](../../../src/conditional_moments.cpp#L72) |
| `reconstruct_conditional_moments` | [81](../../../src/conditional_moments.cpp#L81) |

### `src/conditioning_path.cpp`

[Source](../../../src/conditioning_path.cpp) · [Responsibilities, contracts,
and evidence](state_and_representations.md)

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

[Source](../../../src/configuration.cpp) · [Responsibilities, contracts, and
evidence](initialization_and_configuration.md)

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
| `parse_vandermonde_rcond` | [1694](../../../src/configuration.cpp#L1694) |

### `src/distributions.cpp`

[Source](../../../src/distributions.cpp) · [Responsibilities, contracts, and
evidence](initialization_and_configuration.md)

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

[Source](../../../src/internal/bjorck_pereyra.cpp) · [Responsibilities,
contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_all_finite` | [16](../../../src/internal/bjorck_pereyra.cpp#L16) |
| `failure_result` | [31](../../../src/internal/bjorck_pereyra.cpp#L31) |
| `node_scale` | [36](../../../src/internal/bjorck_pereyra.cpp#L36) |
| `denominator_is_singular` | [41](../../../src/internal/bjorck_pereyra.cpp#L41) |
| `solve_bjorck_pereyra` | [48](../../../src/internal/bjorck_pereyra.cpp#L48) |

### `src/internal/cqmom_active_counts.hpp`

[Source](../../../src/internal/cqmom_active_counts.hpp) · [Responsibilities,
contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `parent_path_at_level` | [14](../../../src/internal/cqmom_active_counts.hpp#L14) |
| `find_active_branch_accumulator` | [33](../../../src/internal/cqmom_active_counts.hpp#L33) |
| `contains_node_index` | [49](../../../src/internal/cqmom_active_counts.hpp#L49) |
| `validate_contiguous_child_indices` | [54](../../../src/internal/cqmom_active_counts.hpp#L54) |
| `active_branch_counts_from_quadrature` | [67](../../../src/internal/cqmom_active_counts.hpp#L67) |

### `src/internal/cqmom_validation.hpp`

[Source](../../../src/internal/cqmom_validation.hpp) · [Responsibilities,
contracts, and evidence](utilities.md)

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
| `cqmom_required_moment_extents_representable` | [18](../../../src/internal/cqmom_validation.hpp#L18) |

### `src/internal/embedded_runge_kutta.cpp`

[Source](../../../src/internal/embedded_runge_kutta.cpp) · [Responsibilities,
contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `embedded_rk_tableau` | [47](../../../src/internal/embedded_runge_kutta.cpp#L47) |
| `take_embedded_rk_step` | [59](../../../src/internal/embedded_runge_kutta.cpp#L59) |
| `embedded_rk_controller_exponent` | [112](../../../src/internal/embedded_runge_kutta.cpp#L112) |

### `src/internal/householder_qr.cpp`

[Source](../../../src/internal/householder_qr.cpp) · [Responsibilities,
contracts, and evidence](conditional_systems.md)

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

[Source](../../../src/internal/jacobi_svd.cpp) · [Responsibilities, contracts,
and evidence](conditional_systems.md)

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

[Source](../../../src/internal/linear_algebra.cpp) · [Responsibilities,
contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_infinity_norm` | [33](../../../src/internal/linear_algebra.cpp#L33) |
| `solve_dense` | [48](../../../src/internal/linear_algebra.cpp#L48) |
| `make_identity` | [124](../../../src/internal/linear_algebra.cpp#L124) |
| `residual_norm` | [134](../../../src/internal/linear_algebra.cpp#L134) |
| `condition_estimate` | [154](../../../src/internal/linear_algebra.cpp#L154) |
| `rotate` | [168](../../../src/internal/linear_algebra.cpp#L168) |
| `symmetric_jacobi_eigen` | [213](../../../src/internal/linear_algebra.cpp#L213) |
| `matrix_all_finite` | [16](../../../src/internal/linear_algebra.cpp#L16) |

### `src/internal/linear_algebra.hpp`

[Source](../../../src/internal/linear_algebra.hpp) · [Responsibilities,
contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `DenseMatrix::DenseMatrix` | [18](../../../src/internal/linear_algebra.hpp#L18) |
| `DenseMatrix::rows` | [20](../../../src/internal/linear_algebra.hpp#L20) |
| `DenseMatrix::cols` | [25](../../../src/internal/linear_algebra.hpp#L25) |
| `DenseMatrix::operator()` | [30](../../../src/internal/linear_algebra.hpp#L30), [39](../../../src/internal/linear_algebra.hpp#L39) |
| `DenseMatrix::storage_size` | [49](../../../src/internal/linear_algebra.hpp#L49) |

### `src/internal/numeric_utils.hpp`

[Source](../../../src/internal/numeric_utils.hpp) · [Responsibilities,
contracts, and evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `all_finite` | [12](../../../src/internal/numeric_utils.hpp#L12) |
| `output_buffer_poison` | [30](../../../src/internal/numeric_utils.hpp#L30) |
| `is_output_buffer_poison` | [37](../../../src/internal/numeric_utils.hpp#L37) |
| `poison_output_buffer` | [47](../../../src/internal/numeric_utils.hpp#L47) |

### `src/internal/ode_rhs_validation.hpp`

[Source](../../../src/internal/ode_rhs_validation.hpp) · [Responsibilities,
contracts, and evidence](time_integration.md)

| Function or overload family | Definition lines |
|---|---|
| `ode_rhs_unfilled_message` | [16](../../../src/internal/ode_rhs_validation.hpp#L16) |
| `rhs_checked` | [35](../../../src/internal/ode_rhs_validation.hpp#L35) |

### `src/internal/performance.hpp`

[Source](../../../src/internal/performance.hpp) · [Responsibilities, contracts,
and evidence](utilities.md)

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

[Source](../../../src/internal/phase_space_velocity_validation.hpp) ·
[Responsibilities, contracts, and evidence](utilities.md)

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

[Source](../../../src/internal/qmom_state_access.hpp) · [Responsibilities,
contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomStateAccess::from_quadrature` | [20](../../../src/internal/qmom_state_access.hpp#L20) |

### `src/internal/recursive_quadrature.hpp`

[Source](../../../src/internal/recursive_quadrature.hpp) · [Responsibilities,
contracts, and evidence](cqmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `RecursiveQuadratureResult::ok` | [318](../../../src/internal/recursive_quadrature.hpp#L318) |

### `src/internal/vandermonde.cpp`

[Source](../../../src/internal/vandermonde.cpp) · [Responsibilities, contracts,
and evidence](conditional_systems.md)

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

[Source](../../../src/moment_reconstruction.cpp) · [Responsibilities,
contracts, and evidence](moment_evaluation.md)

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

[Source](../../../src/moment_set.cpp) · [Responsibilities, contracts, and
evidence](qmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `MomentSet::MomentSet` | [9](../../../src/moment_set.cpp#L9), [11](../../../src/moment_set.cpp#L11) |
| `MomentSet::size` | [13](../../../src/moment_set.cpp#L13) |
| `MomentSet::empty` | [18](../../../src/moment_set.cpp#L18) |
| `MomentSet::values` | [23](../../../src/moment_set.cpp#L23) |
| `MomentSet::operator[]` | [28](../../../src/moment_set.cpp#L28) |

### `src/moment_source.cpp`

[Source](../../../src/moment_source.cpp) · [Responsibilities, contracts, and
evidence](moment_sources.md)

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
| `advance_qmom_moments_with_source` | [399](../../../src/moment_source.cpp#L399), [813](../../../src/moment_source.cpp#L813), [833](../../../src/moment_source.cpp#L833) |
| `compose_cqmom_source_and_boundary` | [452](../../../src/moment_source.cpp#L452) |
| `advance_validated_cqmom_moments_with_source` | [479](../../../src/moment_source.cpp#L479) |
| `QmomSourceContextAccess::create` | [543](../../../src/moment_source.cpp#L543) |
| `CqmomSourceContextAccess::create` | [543](../../../src/moment_source.cpp#L543) |
| `QmomSourceContext::QmomSourceContext` | [559](../../../src/moment_source.cpp#L559), [567](../../../src/moment_source.cpp#L567) |
| `QmomSourceContext::time` | [580](../../../src/moment_source.cpp#L580) |
| `QmomSourceContext::moment_count` | [585](../../../src/moment_source.cpp#L585) |
| `QmomSourceContext::size` | [590](../../../src/moment_source.cpp#L590) |
| `QmomSourceContext::values` | [595](../../../src/moment_source.cpp#L595) |
| `QmomSourceContext::value` | [600](../../../src/moment_source.cpp#L600) |
| `QmomSourceContext::moment` | [605](../../../src/moment_source.cpp#L605) |
| `QmomSourceContext::generalized_moment` | [622](../../../src/moment_source.cpp#L622) |
| `CqmomSourceContext::CqmomSourceContext` | [646](../../../src/moment_source.cpp#L646), [657](../../../src/moment_source.cpp#L657) |
| `CqmomSourceContext::time` | [671](../../../src/moment_source.cpp#L671) |
| `CqmomSourceContext::moment_count` | [676](../../../src/moment_source.cpp#L676) |
| `CqmomSourceContext::dimension_count` | [681](../../../src/moment_source.cpp#L681) |
| `CqmomSourceContext::size` | [686](../../../src/moment_source.cpp#L686) |
| `CqmomSourceContext::indices` | [691](../../../src/moment_source.cpp#L691) |
| `CqmomSourceContext::values` | [696](../../../src/moment_source.cpp#L696) |
| `CqmomSourceContext::index` | [701](../../../src/moment_source.cpp#L701) |
| `CqmomSourceContext::value` | [706](../../../src/moment_source.cpp#L706), [721](../../../src/moment_source.cpp#L721) |
| `CqmomSourceContext::position` | [711](../../../src/moment_source.cpp#L711) |
| `CqmomSourceContext::moment` | [731](../../../src/moment_source.cpp#L731) |
| `CqmomSourceContext::generalized_moment` | [747](../../../src/moment_source.cpp#L747) |
| `CqmomSourceContext::for_each_terminal_node` | [760](../../../src/moment_source.cpp#L760) |
| `CqmomSourceContext::reconstructed_stage_quadrature` | [781](../../../src/moment_source.cpp#L781) |
| `CqmomSourceContext::hierarchy` | [808](../../../src/moment_source.cpp#L808) |
| `advance_cqmom_moments_with_source` | [849](../../../src/moment_source.cpp#L849), [864](../../../src/moment_source.cpp#L864) |
| `advance_trusted_qmom_moments_with_source` | [883](../../../src/moment_source.cpp#L883), [908](../../../src/moment_source.cpp#L908) |
| `advance_trusted_cqmom_moments_with_source` | [935](../../../src/moment_source.cpp#L935), [950](../../../src/moment_source.cpp#L950), [969](../../../src/moment_source.cpp#L969), [987](../../../src/moment_source.cpp#L987) |
| `assign_source_rates` | [394](../../../src/moment_source.cpp#L394) |

### `src/node_advection.cpp`

[Source](../../../src/node_advection.cpp) · [Responsibilities, contracts, and
evidence](time_integration.md)

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

[Source](../../../src/output.cpp) · [Responsibilities, contracts, and
evidence](utilities.md)

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
| `centered_column_label` | [24](../../../src/output.cpp#L24) |
| `generalized_moment_label` | [68](../../../src/output.cpp#L68) |
| `append_generalized_labels` | [84](../../../src/output.cpp#L84) |
| `validate_generalized_moment_orders` | [93](../../../src/output.cpp#L93) |

### `src/phase_space_evolution.cpp`

[Source](../../../src/phase_space_evolution.cpp) · [Responsibilities,
contracts, and evidence](phase_space_evolution.md)

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
| `advance_qmom_moment_source` | [753](../../../src/phase_space_evolution.cpp#L753), [759](../../../src/phase_space_evolution.cpp#L759) |
| `advance_qmom_phase_space` | [766](../../../src/phase_space_evolution.cpp#L766), [773](../../../src/phase_space_evolution.cpp#L773) |
| `initialize_cqmom_diagnostics` | [783](../../../src/phase_space_evolution.cpp#L783), [790](../../../src/phase_space_evolution.cpp#L790), [797](../../../src/phase_space_evolution.cpp#L797) |
| `validate_dense_cqmom_moments_for_phase_space_evolution` | [810](../../../src/phase_space_evolution.cpp#L810) |
| `validate_lean_cqmom_moments_for_phase_space_evolution` | [817](../../../src/phase_space_evolution.cpp#L817) |
| `validate_cqmom_state_for_phase_space_evolution` | [824](../../../src/phase_space_evolution.cpp#L824) |
| `validate_cqmom_advection_input` | [839](../../../src/phase_space_evolution.cpp#L839) |
| `validate_cqmom_moment_source_input` | [849](../../../src/phase_space_evolution.cpp#L849) |
| `validate_cqmom_phase_space_evolution_input` | [861](../../../src/phase_space_evolution.cpp#L861) |
| `CqmomStageOutcome::ok` | [888](../../../src/phase_space_evolution.cpp#L888) |
| `take_cqmom_advection_active_counts` | [902](../../../src/phase_space_evolution.cpp#L902) |
| `cqmom_state_failure` | [911](../../../src/phase_space_evolution.cpp#L911) |
| `cqmom_success` | [919](../../../src/phase_space_evolution.cpp#L919) |
| `materialize_cqmom_phase_space_evolution_result` | [924](../../../src/phase_space_evolution.cpp#L924) |
| `dense_state_from_accepted_quadrature` | [945](../../../src/phase_space_evolution.cpp#L945) |
| `lean_state_from_accepted_quadrature` | [961](../../../src/phase_space_evolution.cpp#L961) |
| `prepare_cqmom_phase_space_evolution_working_state` | [972](../../../src/phase_space_evolution.cpp#L972) |
| `advance_cqmom_advection_stage` | [1015](../../../src/phase_space_evolution.cpp#L1015) |
| `advance_cqmom_moment_source_stage` | [1111](../../../src/phase_space_evolution.cpp#L1111) |
| `advance_cqmom_advection_impl` | [1248](../../../src/phase_space_evolution.cpp#L1248) |
| `advance_cqmom_moment_source_impl` | [1255](../../../src/phase_space_evolution.cpp#L1255) |
| `advance_cqmom_phase_space_state_impl` | [1262](../../../src/phase_space_evolution.cpp#L1262) |
| `advance_cqmom_advection` | [1322](../../../src/phase_space_evolution.cpp#L1322) |
| `advance_cqmom_moment_source_public_impl` | [1344](../../../src/phase_space_evolution.cpp#L1344) |
| `advance_cqmom_moment_source` | [1367](../../../src/phase_space_evolution.cpp#L1367), [1373](../../../src/phase_space_evolution.cpp#L1373) |
| `advance_cqmom_phase_space` | [1380](../../../src/phase_space_evolution.cpp#L1380), [1387](../../../src/phase_space_evolution.cpp#L1387) |
| `advance_qmom_advection_batch` | [1395](../../../src/phase_space_evolution.cpp#L1395) |
| `advance_cqmom_advection_batch` | [1403](../../../src/phase_space_evolution.cpp#L1403) |
| `advance_qmom_moment_source_batch` | [1411](../../../src/phase_space_evolution.cpp#L1411), [1427](../../../src/phase_space_evolution.cpp#L1427) |
| `advance_cqmom_moment_source_batch` | [1419](../../../src/phase_space_evolution.cpp#L1419), [1436](../../../src/phase_space_evolution.cpp#L1436) |
| `advance_qmom_phase_space_batch` | [1445](../../../src/phase_space_evolution.cpp#L1445), [1463](../../../src/phase_space_evolution.cpp#L1463) |
| `advance_cqmom_phase_space_batch` | [1454](../../../src/phase_space_evolution.cpp#L1454), [1473](../../../src/phase_space_evolution.cpp#L1473) |
| `advance_qmom_advection_batch_with_index` | [1483](../../../src/phase_space_evolution.cpp#L1483) |
| `advance_cqmom_advection_batch_with_index` | [1500](../../../src/phase_space_evolution.cpp#L1500) |
| `advance_qmom_moment_source_batch_with_index` | [1517](../../../src/phase_space_evolution.cpp#L1517), [1547](../../../src/phase_space_evolution.cpp#L1547) |
| `advance_cqmom_moment_source_batch_with_index` | [1532](../../../src/phase_space_evolution.cpp#L1532), [1568](../../../src/phase_space_evolution.cpp#L1568) |
| `advance_qmom_phase_space_batch_with_index` | [1589](../../../src/phase_space_evolution.cpp#L1589), [1634](../../../src/phase_space_evolution.cpp#L1634) |
| `advance_cqmom_phase_space_batch_with_index` | [1611](../../../src/phase_space_evolution.cpp#L1611), [1662](../../../src/phase_space_evolution.cpp#L1662) |
| `validate_cqmom_moment_values_match_layout` | [819](../../../src/phase_space_evolution.cpp#L819) |

### `src/public_api.cpp`

[Source](../../../src/public_api.cpp) · [Responsibilities, contracts, and
evidence](state_and_representations.md)

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
| `CqmomState::CqmomState` | [365](../../../src/public_api.cpp#L365), [367](../../../src/public_api.cpp#L367), [378](../../../src/public_api.cpp#L378) |
| `CqmomState::operator=` | [369](../../../src/public_api.cpp#L369), [380](../../../src/public_api.cpp#L380) |
| `CqmomState::hierarchy` | [384](../../../src/public_api.cpp#L384) |
| `CqmomState::quadrature` | [389](../../../src/public_api.cpp#L389) |
| `CqmomState::qmom_options` | [398](../../../src/public_api.cpp#L398) |
| `CqmomState::conditional_options` | [403](../../../src/public_api.cpp#L403) |
| `CqmomState::adaptive_branch_counts` | [408](../../../src/public_api.cpp#L408) |
| `CqmomState::diagnostic_level` | [413](../../../src/public_api.cpp#L413) |
| `CqmomState::vandermonde_solver` | [418](../../../src/public_api.cpp#L418) |
| `CqmomState::coordinate_support` | [423](../../../src/public_api.cpp#L423) |
| `CqmomState::moment_storage` | [428](../../../src/public_api.cpp#L428) |
| `CqmomState::active_branch_counts` | [433](../../../src/public_api.cpp#L433) |
| `CqmomState::moment_indices` | [438](../../../src/public_api.cpp#L438) |
| `CqmomState::moment_values` | [445](../../../src/public_api.cpp#L445) |
| `CqmomState::moment` | [452](../../../src/public_api.cpp#L452) |
| `validate_prepared_metadata` | [484](../../../src/public_api.cpp#L484), [496](../../../src/public_api.cpp#L496) |
| `make_prepared_metadata` | [505](../../../src/public_api.cpp#L505), [515](../../../src/public_api.cpp#L515) |
| `CqmomConstructionAccess::create_case` | [532](../../../src/public_api.cpp#L532), [526](../../../src/public_api.cpp#L526) |
| `CqmomConstructionAccess::create_state` | [557](../../../src/public_api.cpp#L557), [549](../../../src/public_api.cpp#L549), [564](../../../src/public_api.cpp#L564), [575](../../../src/public_api.cpp#L575), [569](../../../src/public_api.cpp#L569), [583](../../../src/public_api.cpp#L583), [543](../../../src/public_api.cpp#L543), [538](../../../src/public_api.cpp#L538) |
| `CqmomConstructionAccess::case_dense_case` | [590](../../../src/public_api.cpp#L590) |
| `CqmomConstructionAccess::case_lean_case` | [595](../../../src/public_api.cpp#L595) |
| `CqmomConstructionAccess::state_dense_case` | [600](../../../src/public_api.cpp#L600) |
| `CqmomConstructionAccess::state_lean_case` | [605](../../../src/public_api.cpp#L605) |
| `CqmomConstructionAccess::state_active_branch_counts` | [610](../../../src/public_api.cpp#L610) |
| `CqmomConstructionAccess::state_accepted_quadrature` | [615](../../../src/public_api.cpp#L615) |
| `CqmomConstructionAccess::case_prepared_metadata` | [620](../../../src/public_api.cpp#L620) |
| `CqmomConstructionAccess::state_prepared_metadata` | [625](../../../src/public_api.cpp#L625) |
| `CqmomConstructionAccess::create_quadrature` | [630](../../../src/public_api.cpp#L630) |
| `CqmomQuadrature::CqmomQuadrature` | [656](../../../src/public_api.cpp#L656) |
| `CqmomQuadrature::hierarchy` | [384](../../../src/public_api.cpp#L384) |
| `CqmomQuadrature::weight` | [647](../../../src/public_api.cpp#L647) |
| `CqmomQuadrature::abscissa` | [656](../../../src/public_api.cpp#L656) |
| `CqmomQuadrature::node` | [669](../../../src/public_api.cpp#L669) |
| `to_cqmom_status` | [682](../../../src/public_api.cpp#L682) |
| `to_cqmom_diagnostics` | [740](../../../src/public_api.cpp#L740), [698](../../../src/public_api.cpp#L698) |
| `terminal_paths_for_public_quadrature` | [752](../../../src/public_api.cpp#L752) |
| `failed_phase_space_dimension` | [774](../../../src/public_api.cpp#L774) |
| `to_cqmom_result` | [784](../../../src/public_api.cpp#L784) |
| `compute_cqmom` | [815](../../../src/public_api.cpp#L815), [810](../../../src/public_api.cpp#L810), [962](../../../src/public_api.cpp#L962) |
| `compute_cqmom_internal` | [924](../../../src/public_api.cpp#L924), [863](../../../src/public_api.cpp#L863) |
| `invalid_cqmom_hierarchy_result` | [800](../../../src/public_api.cpp#L800) |

### `src/realizability.cpp`

[Source](../../../src/realizability.cpp) · [Responsibilities, contracts, and
evidence](realizability.md)

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
| `numerical_rank` | [352](../../../src/realizability.cpp#L352) |

### `src/reconstruction_scaling.cpp`

[Source](../../../src/reconstruction_scaling.cpp) · [Responsibilities,
contracts, and evidence](qmom_reconstruction.md)

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

[Source](../../../src/recursive_quadrature.cpp) · [Responsibilities, contracts,
and evidence](cqmom_reconstruction.md)

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
| `append_unique_index` | [223](../../../src/recursive_quadrature.cpp#L223) |
| `local_branch_offset` | [243](../../../src/recursive_quadrature.cpp#L243) |
| `create_child_layout_case` | [272](../../../src/recursive_quadrature.cpp#L272) |
| `create_branch_moment_layout` | [320](../../../src/recursive_quadrature.cpp#L320) |
| `PreparedCqmomMetadata::PreparedCqmomMetadata` | [383](../../../src/recursive_quadrature.cpp#L383) |
| `create_prepared_cqmom_metadata` | [395](../../../src/recursive_quadrature.cpp#L395) |
| `prepared_cqmom_metadata_is_compatible` | [410](../../../src/recursive_quadrature.cpp#L410) |
| `extract_wheeler_moments` | [425](../../../src/recursive_quadrature.cpp#L425) |
| `child_layout_for_parent_node_count` | [444](../../../src/recursive_quadrature.cpp#L444) |
| `make_branch_rhs_matrix` | [458](../../../src/recursive_quadrature.cpp#L458) |
| `validate_positive_extents` | [487](../../../src/recursive_quadrature.cpp#L487) |
| `checked_multiply` | [499](../../../src/recursive_quadrature.cpp#L499) |
| `default_conditioning_order` | [509](../../../src/recursive_quadrature.cpp#L509) |
| `validate_conditioning_order` | [516](../../../src/recursive_quadrature.cpp#L516) |
| `ordered_extents` | [538](../../../src/recursive_quadrature.cpp#L538) |
| `from_ordered_index` | [549](../../../src/recursive_quadrature.cpp#L549) |
| `validate_abscissas` | [559](../../../src/recursive_quadrature.cpp#L559) |
| `valid_conditional_options` | [600](../../../src/recursive_quadrature.cpp#L600) |
| `valid_diagnostic_level` | [607](../../../src/recursive_quadrature.cpp#L607) |
| `records_success_branch_diagnostics` | [619](../../../src/recursive_quadrature.cpp#L619) |
| `create_ordered_root_slice` | [624](../../../src/recursive_quadrature.cpp#L624) |
| `one_dimensional_moments` | [639](../../../src/recursive_quadrature.cpp#L639) |
| `reconstruct_child_slices` | [644](../../../src/recursive_quadrature.cpp#L644) |
| `create_recursive_result` | [736](../../../src/recursive_quadrature.cpp#L736) |
| `record_failure` | [741](../../../src/recursive_quadrature.cpp#L741) |
| `validate_recursive_request_preflight` | [750](../../../src/recursive_quadrature.cpp#L750) |
| `resolve_cqmom_wheeler_options` | [774](../../../src/recursive_quadrature.cpp#L774) |
| `operator<` | [790](../../../src/recursive_quadrature.cpp#L790) |
| `path_has_prefix` | [815](../../../src/recursive_quadrature.cpp#L815) |
| `branch_limit_invalidated` | [831](../../../src/recursive_quadrature.cpp#L831) |
| `validate_active_branch_topology` | [843](../../../src/recursive_quadrature.cpp#L843) |
| `recursive_branch_node_count_limit` | [915](../../../src/recursive_quadrature.cpp#L915) |
| `compute_branch_quadrature` | [933](../../../src/recursive_quadrature.cpp#L933) |
| `validate_coordinate_support` | [938](../../../src/recursive_quadrature.cpp#L938) |
| `support_interval_for_dimension` | [965](../../../src/recursive_quadrature.cpp#L965) |
| `record_bounded_candidate_diagnostics` | [975](../../../src/recursive_quadrature.cpp#L975) |
| `branch_growth_trial_allowed` | [993](../../../src/recursive_quadrature.cpp#L993) |
| `product_extents` | [1007](../../../src/recursive_quadrature.cpp#L1007) |
| `flatten_multi_index` | [1019](../../../src/recursive_quadrature.cpp#L1019) |
| `unflatten_multi_index` | [1039](../../../src/recursive_quadrature.cpp#L1039) |
| `enumerate_multi_indices` | [1189](../../../src/recursive_quadrature.cpp#L1189) |
| `PhaseSpaceHierarchy::PhaseSpaceHierarchy` | [1201](../../../src/recursive_quadrature.cpp#L1201), [1208](../../../src/recursive_quadrature.cpp#L1208) |
| `PhaseSpaceHierarchy::number_of_dimensions` | [1458](../../../src/recursive_quadrature.cpp#L1458) |
| `PhaseSpaceHierarchy::node_count` | [1201](../../../src/recursive_quadrature.cpp#L1201) |
| `PhaseSpaceHierarchy::total_nodes` | [1463](../../../src/recursive_quadrature.cpp#L1463) |
| `PhaseSpaceHierarchy::moment_extents` | [1223](../../../src/recursive_quadrature.cpp#L1223) |
| `PhaseSpaceHierarchy::flatten_node_index` | [1234](../../../src/recursive_quadrature.cpp#L1234) |
| `PhaseSpaceHierarchy::unflatten_node_index` | [1243](../../../src/recursive_quadrature.cpp#L1243) |
| `create_cqmom_moment_layout` | [1248](../../../src/recursive_quadrature.cpp#L1248) |
| `CqmomMomentData::CqmomMomentData` | [1310](../../../src/recursive_quadrature.cpp#L1310) |
| `CqmomMomentData::contains_required_index` | [1315](../../../src/recursive_quadrature.cpp#L1315) |
| `CqmomMomentData::moment` | [1223](../../../src/recursive_quadrature.cpp#L1223) |
| `extract_cqmom_moments` | [1336](../../../src/recursive_quadrature.cpp#L1336) |
| `create_cqmom_moment_data` | [1363](../../../src/recursive_quadrature.cpp#L1363) |
| `enumerate_quadrature_node_multi_indices` | [1385](../../../src/recursive_quadrature.cpp#L1385) |
| `rectangular_terminal_paths` | [1400](../../../src/recursive_quadrature.cpp#L1400) |
| `RecursiveQuadrature::RecursiveQuadrature` | [1435](../../../src/recursive_quadrature.cpp#L1435), [1442](../../../src/recursive_quadrature.cpp#L1442) |
| `RecursiveQuadrature::hierarchy` | [1385](../../../src/recursive_quadrature.cpp#L1385) |
| `RecursiveQuadrature::number_of_dimensions` | [1458](../../../src/recursive_quadrature.cpp#L1458) |
| `RecursiveQuadrature::total_nodes` | [1463](../../../src/recursive_quadrature.cpp#L1463) |
| `RecursiveQuadrature::weights` | [1442](../../../src/recursive_quadrature.cpp#L1442) |
| `RecursiveQuadrature::abscissas` | [1442](../../../src/recursive_quadrature.cpp#L1442) |
| `RecursiveQuadrature::adaptive_branch_counts_enabled` | [1478](../../../src/recursive_quadrature.cpp#L1478) |
| `RecursiveQuadrature::has_rectangular_storage` | [1483](../../../src/recursive_quadrature.cpp#L1483) |
| `RecursiveQuadrature::tree_nodes` | [1488](../../../src/recursive_quadrature.cpp#L1488) |
| `RecursiveQuadrature::terminal_node_ids` | [1493](../../../src/recursive_quadrature.cpp#L1493) |
| `RecursiveQuadrature::terminal_paths` | [1400](../../../src/recursive_quadrature.cpp#L1400) |
| `RecursiveQuadrature::weight` | [1442](../../../src/recursive_quadrature.cpp#L1442) |
| `RecursiveQuadrature::abscissa` | [1442](../../../src/recursive_quadrature.cpp#L1442) |
| `RecursiveQuadrature::node` | [1404](../../../src/recursive_quadrature.cpp#L1404) |
| `RecursiveQuadrature::set_weight` | [1535](../../../src/recursive_quadrature.cpp#L1535) |
| `RecursiveQuadrature::set_abscissa` | [1544](../../../src/recursive_quadrature.cpp#L1544) |
| `RecursiveQuadrature::set_adaptive_branch_counts_enabled` | [1557](../../../src/recursive_quadrature.cpp#L1557) |
| `RecursiveQuadrature::add_tree_node` | [1570](../../../src/recursive_quadrature.cpp#L1570) |
| `RecursiveQuadrature::set_terminal_quadrature` | [1582](../../../src/recursive_quadrature.cpp#L1582) |
| `compute_recursive_quadrature_from_moment_view` | [1608](../../../src/recursive_quadrature.cpp#L1608) |
| `compute_recursive_quadrature_from_validated_moment_view` | [1948](../../../src/recursive_quadrature.cpp#L1948) |
| `compute_recursive_quadrature_from_moment_data` | [1957](../../../src/recursive_quadrature.cpp#L1957) |
| `compute_recursive_quadrature_from_trusted_moment_data` | [1991](../../../src/recursive_quadrature.cpp#L1991) |
| `compute_recursive_quadrature` | [1948](../../../src/recursive_quadrature.cpp#L1948) |
| `operator==` | [1058](../../../src/recursive_quadrature.cpp#L1058), [1068](../../../src/recursive_quadrature.cpp#L1068), [1073](../../../src/recursive_quadrature.cpp#L1073), [1179](../../../src/recursive_quadrature.cpp#L1179) |
| `operator!=` | [1063](../../../src/recursive_quadrature.cpp#L1063), [1078](../../../src/recursive_quadrature.cpp#L1078), [1083](../../../src/recursive_quadrature.cpp#L1083), [1184](../../../src/recursive_quadrature.cpp#L1184) |
| `MultiIndexTable::MultiIndexTable` | [1088](../../../src/recursive_quadrature.cpp#L1088), [1104](../../../src/recursive_quadrature.cpp#L1104) |
| `MultiIndexTable::at` | [1019](../../../src/recursive_quadrature.cpp#L1019) |
| `MultiIndexTable::push_back` | [1131](../../../src/recursive_quadrature.cpp#L1131), [1136](../../../src/recursive_quadrature.cpp#L1136) |
| `MultiIndexTable::pop_back` | [1149](../../../src/recursive_quadrature.cpp#L1149) |
| `MultiIndexTable::reserve` | [1158](../../../src/recursive_quadrature.cpp#L1158) |
| `MultiIndexTable::clear` | [1163](../../../src/recursive_quadrature.cpp#L1163) |
| `MultiIndexTable::to_multi_indices` | [1168](../../../src/recursive_quadrature.cpp#L1168) |
| `validate_cqmom_moment_data_matches_hierarchy` | [1972](../../../src/recursive_quadrature.cpp#L1972) |

### `src/runge_kutta.cpp`

[Source](../../../src/runge_kutta.cpp) · [Responsibilities, contracts, and
evidence](time_integration.md)

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

[Source](../../../src/wheeler.cpp) · [Responsibilities, contracts, and
evidence](qmom_reconstruction.md)

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

These headers introduce declarations, aliases, option/result structures, or
umbrella includes rather than additional explicit executable bodies in the
inspected build. Their callable definitions are listed above. Aggregate
initialization is covered with the owning representation or options contract.

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
