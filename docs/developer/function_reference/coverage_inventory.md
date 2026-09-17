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

The inventory includes **1016 definition locations**, grouped by named function
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
| `CqmomQuadrature::number_of_dimensions` | [49](../../../include/qmomkit/cqmom.hpp#L49) |
| `CqmomQuadrature::total_nodes` | [53](../../../include/qmomkit/cqmom.hpp#L53) |
| `CqmomQuadrature::weights` | [57](../../../include/qmomkit/cqmom.hpp#L57) |
| `CqmomQuadrature::abscissas` | [61](../../../include/qmomkit/cqmom.hpp#L61) |
| `CqmomQuadrature::adaptive_branch_counts_enabled` | [65](../../../include/qmomkit/cqmom.hpp#L65) |
| `CqmomQuadrature::has_rectangular_storage` | [69](../../../include/qmomkit/cqmom.hpp#L69) |
| `CqmomQuadrature::terminal_paths` | [73](../../../include/qmomkit/cqmom.hpp#L73) |
| `CqmomResult::ok` | [336](../../../include/qmomkit/cqmom.hpp#L336) |
| `CqmomQuadrature::hierarchy` | [45](../../../include/qmomkit/cqmom.hpp#L45) |

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
| `QmomMomentSourceUpdateResult::ok` | [269](../../../include/qmomkit/moment_source.hpp#L269) |
| `CqmomMomentSourceUpdateResult::ok` | [269](../../../include/qmomkit/moment_source.hpp#L269) |

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

### `src/internal/cqmom_case.hpp`

[Source](../../../src/internal/cqmom_case.hpp) · [Responsibilities,
contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `with_moments` | [30](../../../src/internal/cqmom_case.hpp#L30), [37](../../../src/internal/cqmom_case.hpp#L37) |
| `lean_case_with_moments` | [46](../../../src/internal/cqmom_case.hpp#L46) |

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
| `RecursiveQuadratureResult::ok` | [316](../../../src/internal/recursive_quadrature.hpp#L316) |

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
| `status_from_ssp_rk3` | [24](../../../src/moment_source.cpp#L24) |
| `message_from_status` | [38](../../../src/moment_source.cpp#L38) |
| `diagnostics_from_ssp_rk3` | [56](../../../src/moment_source.cpp#L56) |
| `invalid_result` | [67](../../../src/moment_source.cpp#L67) |
| `invalid_moment_data_result` | [74](../../../src/moment_source.cpp#L74) |
| `failed_moment_data_result` | [81](../../../src/moment_source.cpp#L81) |
| `validate_source_context` | [93](../../../src/moment_source.cpp#L93) |
| `validate_source_generalized_order` | [108](../../../src/moment_source.cpp#L108) |
| `validate_fixed_cqmom_moment_data` | [120](../../../src/moment_source.cpp#L120) |
| `validate_rate` | [127](../../../src/moment_source.cpp#L127) |
| `validate_source_rate` | [134](../../../src/moment_source.cpp#L134) |
| `validate_boundary_rate` | [146](../../../src/moment_source.cpp#L146) |
| `validate_cqmom_boundary_rate` | [158](../../../src/moment_source.cpp#L158) |
| `validate_cqmom_stage_state` | [165](../../../src/moment_source.cpp#L165) |
| `validate_coordinate_realizability_policies` | [177](../../../src/moment_source.cpp#L177) |
| `has_coordinate_realizability_policy` | [189](../../../src/moment_source.cpp#L189) |
| `extract_coordinate_marginal` | [196](../../../src/moment_source.cpp#L196) |
| `accepts_coordinate_realizability_policies` | [228](../../../src/moment_source.cpp#L228) |
| `boundary_callback_error` | [249](../../../src/moment_source.cpp#L249) |
| `cqmom_boundary_callback_error` | [254](../../../src/moment_source.cpp#L254) |
| `validate_request` | [259](../../../src/moment_source.cpp#L259), [290](../../../src/moment_source.cpp#L290) |
| `validate_cqmom_source_update_controls` | [300](../../../src/moment_source.cpp#L300) |
| `validate_moment_data_request` | [322](../../../src/moment_source.cpp#L322), [339](../../../src/moment_source.cpp#L339) |
| `validate_trusted_moment_data_request` | [330](../../../src/moment_source.cpp#L330), [349](../../../src/moment_source.cpp#L349) |
| `compose_source_and_boundary` | [359](../../../src/moment_source.cpp#L359) |
| `advance_qmom_moments_with_source` | [405](../../../src/moment_source.cpp#L405), [859](../../../src/moment_source.cpp#L859), [879](../../../src/moment_source.cpp#L879) |
| `compose_cqmom_source_and_boundary` | [458](../../../src/moment_source.cpp#L458) |
| `advance_validated_cqmom_moments_with_source` | [485](../../../src/moment_source.cpp#L485) |
| `QmomSourceContextAccess::create` | [549](../../../src/moment_source.cpp#L549) |
| `CqmomSourceContextAccess::create` | [549](../../../src/moment_source.cpp#L549) |
| `QmomSourceContext::QmomSourceContext` | [565](../../../src/moment_source.cpp#L565), [573](../../../src/moment_source.cpp#L573) |
| `QmomSourceContext::time` | [586](../../../src/moment_source.cpp#L586) |
| `QmomSourceContext::moment_count` | [591](../../../src/moment_source.cpp#L591) |
| `QmomSourceContext::size` | [596](../../../src/moment_source.cpp#L596) |
| `QmomSourceContext::values` | [601](../../../src/moment_source.cpp#L601) |
| `QmomSourceContext::value` | [606](../../../src/moment_source.cpp#L606) |
| `QmomSourceContext::moment` | [611](../../../src/moment_source.cpp#L611) |
| `QmomSourceContext::generalized_moment` | [628](../../../src/moment_source.cpp#L628) |
| `CqmomSourceContext::CqmomSourceContext` | [652](../../../src/moment_source.cpp#L652), [663](../../../src/moment_source.cpp#L663) |
| `CqmomSourceContext::time` | [677](../../../src/moment_source.cpp#L677) |
| `CqmomSourceContext::moment_count` | [682](../../../src/moment_source.cpp#L682) |
| `CqmomSourceContext::dimension_count` | [687](../../../src/moment_source.cpp#L687) |
| `CqmomSourceContext::size` | [692](../../../src/moment_source.cpp#L692) |
| `CqmomSourceContext::indices` | [697](../../../src/moment_source.cpp#L697) |
| `CqmomSourceContext::values` | [702](../../../src/moment_source.cpp#L702) |
| `CqmomSourceContext::index` | [707](../../../src/moment_source.cpp#L707) |
| `CqmomSourceContext::value` | [712](../../../src/moment_source.cpp#L712), [767](../../../src/moment_source.cpp#L767) |
| `CqmomSourceContext::position` | [754](../../../src/moment_source.cpp#L754) |
| `CqmomSourceContext::moment` | [777](../../../src/moment_source.cpp#L777) |
| `CqmomSourceContext::generalized_moment` | [793](../../../src/moment_source.cpp#L793) |
| `CqmomSourceContext::for_each_terminal_node` | [806](../../../src/moment_source.cpp#L806) |
| `CqmomSourceContext::reconstructed_stage_quadrature` | [827](../../../src/moment_source.cpp#L827) |
| `CqmomSourceContext::hierarchy` | [854](../../../src/moment_source.cpp#L854) |
| `advance_cqmom_moments_with_source` | [895](../../../src/moment_source.cpp#L895), [910](../../../src/moment_source.cpp#L910) |
| `advance_trusted_qmom_moments_with_source` | [929](../../../src/moment_source.cpp#L929), [954](../../../src/moment_source.cpp#L954) |
| `advance_trusted_cqmom_moments_with_source` | [981](../../../src/moment_source.cpp#L981), [996](../../../src/moment_source.cpp#L996), [1015](../../../src/moment_source.cpp#L1015), [1033](../../../src/moment_source.cpp#L1033) |
| `assign_source_rates` | [400](../../../src/moment_source.cpp#L400) |

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
| `PooledCoordinateBuffer::PooledCoordinateBuffer` | [326](../../../src/node_advection.cpp#L326) |
| `PooledCoordinateBuffer::~PooledCoordinateBuffer` | [347](../../../src/node_advection.cpp#L347) |
| `PooledCoordinateBuffer::free_list` | [359](../../../src/node_advection.cpp#L359) |
| `make_phase_space_velocity` | [453](../../../src/node_advection.cpp#L453) |
| `advance_nodes` | [370](../../../src/node_advection.cpp#L370), [398](../../../src/node_advection.cpp#L398), [409](../../../src/node_advection.cpp#L409) |
| `format_node_advection_options` | [493](../../../src/node_advection.cpp#L493) |

### `src/output.cpp`

[Source](../../../src/output.cpp) · [Responsibilities, contracts, and
evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `qmom_label` | [19](../../../src/output.cpp#L19) |
| `column_width_for` | [36](../../../src/output.cpp#L36) |
| `write_column_header` | [41](../../../src/output.cpp#L41) |
| `cqmom_label` | [55](../../../src/output.cpp#L55) |
| `sequential_indices` | [108](../../../src/output.cpp#L108) |
| `qmom_labels` | [124](../../../src/output.cpp#L124) |
| `cqmom_labels` | [140](../../../src/output.cpp#L140) |
| `flatten_selected_indices` | [156](../../../src/output.cpp#L156) |
| `all_cqmom_labels` | [167](../../../src/output.cpp#L167) |
| `cqmom_abscissa_label` | [172](../../../src/output.cpp#L172) |
| `qmom_quadrature_labels` | [177](../../../src/output.cpp#L177) |
| `cqmom_quadrature_labels` | [187](../../../src/output.cpp#L187) |
| `MomentWriter::for_qmom` | [211](../../../src/output.cpp#L211), [217](../../../src/output.cpp#L217), [223](../../../src/output.cpp#L223), [233](../../../src/output.cpp#L233) |
| `MomentWriter::for_cqmom` | [244](../../../src/output.cpp#L244), [251](../../../src/output.cpp#L251), [257](../../../src/output.cpp#L257), [268](../../../src/output.cpp#L268), [273](../../../src/output.cpp#L273), [282](../../../src/output.cpp#L282), [292](../../../src/output.cpp#L292) |
| `MomentWriter::MomentWriter` | [297](../../../src/output.cpp#L297), [315](../../../src/output.cpp#L315), [332](../../../src/output.cpp#L332), [348](../../../src/output.cpp#L348) |
| `MomentWriter::layout_mismatch_message` | [364](../../../src/output.cpp#L364) |
| `MomentWriter::write_header` | [390](../../../src/output.cpp#L390) |
| `MomentWriter::write_row` | [395](../../../src/output.cpp#L395) |
| `MomentWriter::write_aligned_row` | [422](../../../src/output.cpp#L422) |
| `MomentWriter::write` | [453](../../../src/output.cpp#L453), [472](../../../src/output.cpp#L472), [495](../../../src/output.cpp#L495), [514](../../../src/output.cpp#L514), [533](../../../src/output.cpp#L533), [556](../../../src/output.cpp#L556), [583](../../../src/output.cpp#L583), [602](../../../src/output.cpp#L602), [616](../../../src/output.cpp#L616), [635](../../../src/output.cpp#L635) |
| `QuadratureWriter::for_qmom` | [649](../../../src/output.cpp#L649) |
| `QuadratureWriter::for_cqmom` | [654](../../../src/output.cpp#L654) |
| `QuadratureWriter::QuadratureWriter` | [659](../../../src/output.cpp#L659) |
| `QuadratureWriter::write_header` | [675](../../../src/output.cpp#L675) |
| `QuadratureWriter::write_row` | [680](../../../src/output.cpp#L680) |
| `QuadratureWriter::write_qmom_row` | [707](../../../src/output.cpp#L707) |
| `QuadratureWriter::write` | [737](../../../src/output.cpp#L737), [742](../../../src/output.cpp#L742), [747](../../../src/output.cpp#L747) |
| `centered_column_label` | [24](../../../src/output.cpp#L24) |
| `generalized_moment_label` | [71](../../../src/output.cpp#L71) |
| `append_generalized_labels` | [87](../../../src/output.cpp#L87) |
| `validate_generalized_moment_orders` | [96](../../../src/output.cpp#L96) |

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
| `lean_state_from_accepted_quadrature` | [954](../../../src/phase_space_evolution.cpp#L954) |
| `prepare_cqmom_phase_space_evolution_working_state` | [963](../../../src/phase_space_evolution.cpp#L963) |
| `advance_cqmom_advection_stage` | [1006](../../../src/phase_space_evolution.cpp#L1006) |
| `advance_cqmom_moment_source_stage` | [1102](../../../src/phase_space_evolution.cpp#L1102) |
| `advance_cqmom_advection_impl` | [1235](../../../src/phase_space_evolution.cpp#L1235) |
| `advance_cqmom_moment_source_impl` | [1242](../../../src/phase_space_evolution.cpp#L1242) |
| `advance_cqmom_phase_space_state_impl` | [1249](../../../src/phase_space_evolution.cpp#L1249) |
| `advance_cqmom_advection` | [1309](../../../src/phase_space_evolution.cpp#L1309) |
| `advance_cqmom_moment_source_public_impl` | [1331](../../../src/phase_space_evolution.cpp#L1331) |
| `advance_cqmom_moment_source` | [1354](../../../src/phase_space_evolution.cpp#L1354), [1360](../../../src/phase_space_evolution.cpp#L1360) |
| `advance_cqmom_phase_space` | [1367](../../../src/phase_space_evolution.cpp#L1367), [1374](../../../src/phase_space_evolution.cpp#L1374) |
| `advance_qmom_advection_batch` | [1382](../../../src/phase_space_evolution.cpp#L1382) |
| `advance_cqmom_advection_batch` | [1390](../../../src/phase_space_evolution.cpp#L1390) |
| `advance_qmom_moment_source_batch` | [1398](../../../src/phase_space_evolution.cpp#L1398), [1414](../../../src/phase_space_evolution.cpp#L1414) |
| `advance_cqmom_moment_source_batch` | [1406](../../../src/phase_space_evolution.cpp#L1406), [1423](../../../src/phase_space_evolution.cpp#L1423) |
| `advance_qmom_phase_space_batch` | [1432](../../../src/phase_space_evolution.cpp#L1432), [1450](../../../src/phase_space_evolution.cpp#L1450) |
| `advance_cqmom_phase_space_batch` | [1441](../../../src/phase_space_evolution.cpp#L1441), [1460](../../../src/phase_space_evolution.cpp#L1460) |
| `advance_qmom_advection_batch_with_index` | [1470](../../../src/phase_space_evolution.cpp#L1470) |
| `advance_cqmom_advection_batch_with_index` | [1487](../../../src/phase_space_evolution.cpp#L1487) |
| `advance_qmom_moment_source_batch_with_index` | [1504](../../../src/phase_space_evolution.cpp#L1504), [1534](../../../src/phase_space_evolution.cpp#L1534) |
| `advance_cqmom_moment_source_batch_with_index` | [1519](../../../src/phase_space_evolution.cpp#L1519), [1555](../../../src/phase_space_evolution.cpp#L1555) |
| `advance_qmom_phase_space_batch_with_index` | [1576](../../../src/phase_space_evolution.cpp#L1576), [1621](../../../src/phase_space_evolution.cpp#L1621) |
| `advance_cqmom_phase_space_batch_with_index` | [1598](../../../src/phase_space_evolution.cpp#L1598), [1649](../../../src/phase_space_evolution.cpp#L1649) |
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
| `visit_case` | [101](../../../src/public_api.cpp#L101) |
| `hierarchy` | [106](../../../src/public_api.cpp#L106) |
| `qmom_options` | [113](../../../src/public_api.cpp#L113) |
| `conditional_options` | [118](../../../src/public_api.cpp#L118) |
| `adaptive_branch_counts` | [123](../../../src/public_api.cpp#L123) |
| `diagnostic_level` | [128](../../../src/public_api.cpp#L128) |
| `vandermonde_solver` | [133](../../../src/public_api.cpp#L133) |
| `coordinate_support` | [138](../../../src/public_api.cpp#L138) |
| `moment_storage` | [143](../../../src/public_api.cpp#L143) |
| `diagnostic_level_name` | [148](../../../src/public_api.cpp#L148) |
| `vandermonde_solver_name` | [162](../../../src/public_api.cpp#L162) |
| `compute_qmom` | [183](../../../src/public_api.cpp#L183), [188](../../../src/public_api.cpp#L188), [194](../../../src/public_api.cpp#L194) |
| `QmomState::QmomState` | [199](../../../src/public_api.cpp#L199) |
| `QmomCase::initial_state` | [206](../../../src/public_api.cpp#L206) |
| `to_string` | [230](../../../src/public_api.cpp#L230) |
| `format_cqmom_settings` | [242](../../../src/public_api.cpp#L242) |
| `CqmomCase::CqmomCase` | [254](../../../src/public_api.cpp#L254), [256](../../../src/public_api.cpp#L256), [258](../../../src/public_api.cpp#L258), [260](../../../src/public_api.cpp#L260), [271](../../../src/public_api.cpp#L271) |
| `CqmomCase::operator=` | [262](../../../src/public_api.cpp#L262), [273](../../../src/public_api.cpp#L273) |
| `CqmomCase::hierarchy` | [277](../../../src/public_api.cpp#L277) |
| `CqmomCase::qmom_options` | [282](../../../src/public_api.cpp#L282) |
| `CqmomCase::conditional_options` | [287](../../../src/public_api.cpp#L287) |
| `CqmomCase::adaptive_branch_counts` | [292](../../../src/public_api.cpp#L292) |
| `CqmomCase::diagnostic_level` | [297](../../../src/public_api.cpp#L297) |
| `CqmomCase::vandermonde_solver` | [302](../../../src/public_api.cpp#L302) |
| `CqmomCase::coordinate_support` | [307](../../../src/public_api.cpp#L307) |
| `CqmomCase::moment_storage` | [312](../../../src/public_api.cpp#L312) |
| `CqmomCase::initial_state` | [317](../../../src/public_api.cpp#L317) |
| `CqmomState::CqmomState` | [351](../../../src/public_api.cpp#L351), [353](../../../src/public_api.cpp#L353), [364](../../../src/public_api.cpp#L364) |
| `CqmomState::operator=` | [355](../../../src/public_api.cpp#L355), [366](../../../src/public_api.cpp#L366) |
| `CqmomState::hierarchy` | [370](../../../src/public_api.cpp#L370) |
| `CqmomState::quadrature` | [375](../../../src/public_api.cpp#L375) |
| `CqmomState::qmom_options` | [384](../../../src/public_api.cpp#L384) |
| `CqmomState::conditional_options` | [389](../../../src/public_api.cpp#L389) |
| `CqmomState::adaptive_branch_counts` | [394](../../../src/public_api.cpp#L394) |
| `CqmomState::diagnostic_level` | [399](../../../src/public_api.cpp#L399) |
| `CqmomState::vandermonde_solver` | [404](../../../src/public_api.cpp#L404) |
| `CqmomState::coordinate_support` | [409](../../../src/public_api.cpp#L409) |
| `CqmomState::moment_storage` | [414](../../../src/public_api.cpp#L414) |
| `CqmomState::active_branch_counts` | [419](../../../src/public_api.cpp#L419) |
| `CqmomState::moment_indices` | [424](../../../src/public_api.cpp#L424) |
| `CqmomState::moment_values` | [431](../../../src/public_api.cpp#L431) |
| `CqmomState::moment` | [438](../../../src/public_api.cpp#L438) |
| `validate_prepared_metadata` | [470](../../../src/public_api.cpp#L470), [482](../../../src/public_api.cpp#L482) |
| `make_prepared_metadata` | [491](../../../src/public_api.cpp#L491), [501](../../../src/public_api.cpp#L501) |
| `CqmomConstructionAccess::create_case` | [518](../../../src/public_api.cpp#L518), [512](../../../src/public_api.cpp#L512) |
| `CqmomConstructionAccess::create_state` | [543](../../../src/public_api.cpp#L543), [535](../../../src/public_api.cpp#L535), [550](../../../src/public_api.cpp#L550), [561](../../../src/public_api.cpp#L561), [555](../../../src/public_api.cpp#L555), [569](../../../src/public_api.cpp#L569), [529](../../../src/public_api.cpp#L529), [524](../../../src/public_api.cpp#L524) |
| `CqmomConstructionAccess::case_dense_case` | [576](../../../src/public_api.cpp#L576) |
| `CqmomConstructionAccess::case_lean_case` | [581](../../../src/public_api.cpp#L581) |
| `CqmomConstructionAccess::state_dense_case` | [586](../../../src/public_api.cpp#L586) |
| `CqmomConstructionAccess::state_lean_case` | [591](../../../src/public_api.cpp#L591) |
| `CqmomConstructionAccess::state_active_branch_counts` | [596](../../../src/public_api.cpp#L596) |
| `CqmomConstructionAccess::state_accepted_quadrature` | [601](../../../src/public_api.cpp#L601) |
| `CqmomConstructionAccess::case_prepared_metadata` | [606](../../../src/public_api.cpp#L606) |
| `CqmomConstructionAccess::state_prepared_metadata` | [611](../../../src/public_api.cpp#L611) |
| `CqmomConstructionAccess::create_quadrature` | [616](../../../src/public_api.cpp#L616) |
| `CqmomQuadrature::CqmomQuadrature` | [626](../../../src/public_api.cpp#L626) |
| `CqmomQuadrature::hierarchy` | [370](../../../src/public_api.cpp#L370) |
| `CqmomQuadrature::weight` | [633](../../../src/public_api.cpp#L633) |
| `CqmomQuadrature::abscissa` | [642](../../../src/public_api.cpp#L642) |
| `CqmomQuadrature::node` | [655](../../../src/public_api.cpp#L655) |
| `to_cqmom_status` | [668](../../../src/public_api.cpp#L668) |
| `to_cqmom_diagnostics` | [726](../../../src/public_api.cpp#L726), [684](../../../src/public_api.cpp#L684) |
| `terminal_paths_for_public_quadrature` | [738](../../../src/public_api.cpp#L738) |
| `failed_phase_space_dimension` | [760](../../../src/public_api.cpp#L760) |
| `to_cqmom_result` | [770](../../../src/public_api.cpp#L770) |
| `compute_cqmom` | [801](../../../src/public_api.cpp#L801), [796](../../../src/public_api.cpp#L796), [948](../../../src/public_api.cpp#L948) |
| `compute_cqmom_internal` | [910](../../../src/public_api.cpp#L910), [849](../../../src/public_api.cpp#L849) |
| `invalid_cqmom_hierarchy_result` | [786](../../../src/public_api.cpp#L786) |

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
| `IndexedMomentSet::indices` | [70](../../../src/recursive_quadrature.cpp#L70) |
| `IndexedMomentSet::values` | [75](../../../src/recursive_quadrature.cpp#L75) |
| `IndexedMomentSet::size` | [80](../../../src/recursive_quadrature.cpp#L80) |
| `IndexedMomentSet::contains` | [85](../../../src/recursive_quadrature.cpp#L85) |
| `IndexedMomentSet::offset` | [102](../../../src/recursive_quadrature.cpp#L102) |
| `IndexedMomentSet::moment` | [113](../../../src/recursive_quadrature.cpp#L113) |
| `CqmomMomentView::CqmomMomentView` | [118](../../../src/recursive_quadrature.cpp#L118), [120](../../../src/recursive_quadrature.cpp#L120) |
| `CqmomMomentView::moment` | [113](../../../src/recursive_quadrature.cpp#L113) |
| `create_indexed_cqmom_moment_set` | [131](../../../src/recursive_quadrature.cpp#L131) |
| `create_ordered_root_layout` | [139](../../../src/recursive_quadrature.cpp#L139) |
| `append_unique_index` | [213](../../../src/recursive_quadrature.cpp#L213) |
| `local_branch_offset` | [233](../../../src/recursive_quadrature.cpp#L233) |
| `create_child_layout_case` | [262](../../../src/recursive_quadrature.cpp#L262) |
| `create_branch_moment_layout` | [310](../../../src/recursive_quadrature.cpp#L310) |
| `PreparedCqmomMetadata::PreparedCqmomMetadata` | [373](../../../src/recursive_quadrature.cpp#L373) |
| `create_prepared_cqmom_metadata` | [385](../../../src/recursive_quadrature.cpp#L385) |
| `prepared_cqmom_metadata_is_compatible` | [400](../../../src/recursive_quadrature.cpp#L400) |
| `extract_wheeler_moments` | [415](../../../src/recursive_quadrature.cpp#L415) |
| `child_layout_for_parent_node_count` | [434](../../../src/recursive_quadrature.cpp#L434) |
| `make_branch_rhs_matrix` | [448](../../../src/recursive_quadrature.cpp#L448) |
| `validate_positive_extents` | [477](../../../src/recursive_quadrature.cpp#L477) |
| `checked_multiply` | [489](../../../src/recursive_quadrature.cpp#L489) |
| `default_conditioning_order` | [499](../../../src/recursive_quadrature.cpp#L499) |
| `validate_conditioning_order` | [506](../../../src/recursive_quadrature.cpp#L506) |
| `ordered_extents` | [528](../../../src/recursive_quadrature.cpp#L528) |
| `from_ordered_index` | [539](../../../src/recursive_quadrature.cpp#L539) |
| `validate_abscissas` | [549](../../../src/recursive_quadrature.cpp#L549) |
| `valid_conditional_options` | [590](../../../src/recursive_quadrature.cpp#L590) |
| `valid_diagnostic_level` | [597](../../../src/recursive_quadrature.cpp#L597) |
| `records_success_branch_diagnostics` | [609](../../../src/recursive_quadrature.cpp#L609) |
| `create_ordered_root_slice` | [614](../../../src/recursive_quadrature.cpp#L614) |
| `one_dimensional_moments` | [629](../../../src/recursive_quadrature.cpp#L629) |
| `reconstruct_child_slices` | [634](../../../src/recursive_quadrature.cpp#L634) |
| `create_recursive_result` | [726](../../../src/recursive_quadrature.cpp#L726) |
| `record_failure` | [731](../../../src/recursive_quadrature.cpp#L731) |
| `validate_recursive_request_preflight` | [740](../../../src/recursive_quadrature.cpp#L740) |
| `resolve_cqmom_wheeler_options` | [764](../../../src/recursive_quadrature.cpp#L764) |
| `operator<` | [780](../../../src/recursive_quadrature.cpp#L780) |
| `path_has_prefix` | [805](../../../src/recursive_quadrature.cpp#L805) |
| `branch_limit_invalidated` | [821](../../../src/recursive_quadrature.cpp#L821) |
| `validate_active_branch_topology` | [833](../../../src/recursive_quadrature.cpp#L833) |
| `recursive_branch_node_count_limit` | [905](../../../src/recursive_quadrature.cpp#L905) |
| `compute_branch_quadrature` | [923](../../../src/recursive_quadrature.cpp#L923) |
| `validate_coordinate_support` | [928](../../../src/recursive_quadrature.cpp#L928) |
| `support_interval_for_dimension` | [955](../../../src/recursive_quadrature.cpp#L955) |
| `record_bounded_candidate_diagnostics` | [965](../../../src/recursive_quadrature.cpp#L965) |
| `branch_growth_trial_allowed` | [983](../../../src/recursive_quadrature.cpp#L983) |
| `product_extents` | [997](../../../src/recursive_quadrature.cpp#L997) |
| `flatten_multi_index` | [1009](../../../src/recursive_quadrature.cpp#L1009) |
| `unflatten_multi_index` | [1029](../../../src/recursive_quadrature.cpp#L1029) |
| `enumerate_multi_indices` | [1179](../../../src/recursive_quadrature.cpp#L1179) |
| `PhaseSpaceHierarchy::PhaseSpaceHierarchy` | [1191](../../../src/recursive_quadrature.cpp#L1191), [1198](../../../src/recursive_quadrature.cpp#L1198) |
| `PhaseSpaceHierarchy::number_of_dimensions` | [1448](../../../src/recursive_quadrature.cpp#L1448) |
| `PhaseSpaceHierarchy::node_count` | [1204](../../../src/recursive_quadrature.cpp#L1204) |
| `PhaseSpaceHierarchy::total_nodes` | [1453](../../../src/recursive_quadrature.cpp#L1453) |
| `PhaseSpaceHierarchy::moment_extents` | [1213](../../../src/recursive_quadrature.cpp#L1213) |
| `PhaseSpaceHierarchy::flatten_node_index` | [1224](../../../src/recursive_quadrature.cpp#L1224) |
| `PhaseSpaceHierarchy::unflatten_node_index` | [1233](../../../src/recursive_quadrature.cpp#L1233) |
| `create_cqmom_moment_layout` | [1238](../../../src/recursive_quadrature.cpp#L1238) |
| `CqmomMomentData::CqmomMomentData` | [1300](../../../src/recursive_quadrature.cpp#L1300) |
| `CqmomMomentData::contains_required_index` | [1305](../../../src/recursive_quadrature.cpp#L1305) |
| `CqmomMomentData::moment` | [1310](../../../src/recursive_quadrature.cpp#L1310) |
| `extract_cqmom_moments` | [1326](../../../src/recursive_quadrature.cpp#L1326) |
| `create_cqmom_moment_data` | [1353](../../../src/recursive_quadrature.cpp#L1353) |
| `enumerate_quadrature_node_multi_indices` | [1375](../../../src/recursive_quadrature.cpp#L1375) |
| `rectangular_terminal_paths` | [1390](../../../src/recursive_quadrature.cpp#L1390) |
| `RecursiveQuadrature::RecursiveQuadrature` | [1425](../../../src/recursive_quadrature.cpp#L1425), [1432](../../../src/recursive_quadrature.cpp#L1432) |
| `RecursiveQuadrature::hierarchy` | [1443](../../../src/recursive_quadrature.cpp#L1443) |
| `RecursiveQuadrature::number_of_dimensions` | [1448](../../../src/recursive_quadrature.cpp#L1448) |
| `RecursiveQuadrature::total_nodes` | [1453](../../../src/recursive_quadrature.cpp#L1453) |
| `RecursiveQuadrature::weights` | [1458](../../../src/recursive_quadrature.cpp#L1458) |
| `RecursiveQuadrature::abscissas` | [1463](../../../src/recursive_quadrature.cpp#L1463) |
| `RecursiveQuadrature::adaptive_branch_counts_enabled` | [1468](../../../src/recursive_quadrature.cpp#L1468) |
| `RecursiveQuadrature::has_rectangular_storage` | [1473](../../../src/recursive_quadrature.cpp#L1473) |
| `RecursiveQuadrature::tree_nodes` | [1478](../../../src/recursive_quadrature.cpp#L1478) |
| `RecursiveQuadrature::terminal_node_ids` | [1483](../../../src/recursive_quadrature.cpp#L1483) |
| `RecursiveQuadrature::terminal_paths` | [1390](../../../src/recursive_quadrature.cpp#L1390) |
| `RecursiveQuadrature::weight` | [1493](../../../src/recursive_quadrature.cpp#L1493) |
| `RecursiveQuadrature::abscissa` | [1502](../../../src/recursive_quadrature.cpp#L1502) |
| `RecursiveQuadrature::node` | [1515](../../../src/recursive_quadrature.cpp#L1515) |
| `RecursiveQuadrature::set_weight` | [1525](../../../src/recursive_quadrature.cpp#L1525) |
| `RecursiveQuadrature::set_abscissa` | [1534](../../../src/recursive_quadrature.cpp#L1534) |
| `RecursiveQuadrature::set_adaptive_branch_counts_enabled` | [1547](../../../src/recursive_quadrature.cpp#L1547) |
| `RecursiveQuadrature::add_tree_node` | [1560](../../../src/recursive_quadrature.cpp#L1560) |
| `RecursiveQuadrature::set_terminal_quadrature` | [1572](../../../src/recursive_quadrature.cpp#L1572) |
| `compute_recursive_quadrature_from_moment_view` | [1598](../../../src/recursive_quadrature.cpp#L1598) |
| `compute_recursive_quadrature_from_validated_moment_view` | [1938](../../../src/recursive_quadrature.cpp#L1938) |
| `compute_recursive_quadrature_from_moment_data` | [1947](../../../src/recursive_quadrature.cpp#L1947) |
| `compute_recursive_quadrature_from_trusted_moment_data` | [1981](../../../src/recursive_quadrature.cpp#L1981) |
| `compute_recursive_quadrature` | [2013](../../../src/recursive_quadrature.cpp#L2013) |
| `operator==` | [1048](../../../src/recursive_quadrature.cpp#L1048), [1058](../../../src/recursive_quadrature.cpp#L1058), [1063](../../../src/recursive_quadrature.cpp#L1063), [1169](../../../src/recursive_quadrature.cpp#L1169) |
| `operator!=` | [1053](../../../src/recursive_quadrature.cpp#L1053), [1068](../../../src/recursive_quadrature.cpp#L1068), [1073](../../../src/recursive_quadrature.cpp#L1073), [1174](../../../src/recursive_quadrature.cpp#L1174) |
| `MultiIndexTable::MultiIndexTable` | [1078](../../../src/recursive_quadrature.cpp#L1078), [1094](../../../src/recursive_quadrature.cpp#L1094) |
| `MultiIndexTable::at` | [1112](../../../src/recursive_quadrature.cpp#L1112) |
| `MultiIndexTable::push_back` | [1121](../../../src/recursive_quadrature.cpp#L1121), [1126](../../../src/recursive_quadrature.cpp#L1126) |
| `MultiIndexTable::pop_back` | [1139](../../../src/recursive_quadrature.cpp#L1139) |
| `MultiIndexTable::reserve` | [1148](../../../src/recursive_quadrature.cpp#L1148) |
| `MultiIndexTable::clear` | [1153](../../../src/recursive_quadrature.cpp#L1153) |
| `MultiIndexTable::to_multi_indices` | [1158](../../../src/recursive_quadrature.cpp#L1158) |
| `validate_cqmom_moment_data_matches_hierarchy` | [1962](../../../src/recursive_quadrature.cpp#L1962) |

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
