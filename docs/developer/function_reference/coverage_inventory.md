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

The inventory includes **991 definition locations**, grouped by named function
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
| `ConditioningPathView::ConditioningPathView` | [70](../../../include/qmomkit/conditioning_path.hpp#L70) |
| `ConditioningPathView::empty` | [75](../../../include/qmomkit/conditioning_path.hpp#L75) |
| `ConditioningPathView::size` | [79](../../../include/qmomkit/conditioning_path.hpp#L79) |
| `ConditioningPathView::begin` | [86](../../../include/qmomkit/conditioning_path.hpp#L86) |
| `ConditioningPathView::end` | [90](../../../include/qmomkit/conditioning_path.hpp#L90) |
| `ConditioningPathView::data` | [94](../../../include/qmomkit/conditioning_path.hpp#L94) |
| `ConditioningPathView::to_conditioning_path` | [100](../../../include/qmomkit/conditioning_path.hpp#L100) |
| `ConditioningPathTable::const_iterator::operator*` | [143](../../../include/qmomkit/conditioning_path.hpp#L143) |
| `ConditioningPathTable::const_iterator::operator++` | [147](../../../include/qmomkit/conditioning_path.hpp#L147), [152](../../../include/qmomkit/conditioning_path.hpp#L152) |
| `ConditioningPathTable::const_iterator::operator==` | [158](../../../include/qmomkit/conditioning_path.hpp#L158) |
| `ConditioningPathTable::const_iterator::operator!=` | [162](../../../include/qmomkit/conditioning_path.hpp#L162) |
| `ConditioningPathTable::size` | [183](../../../include/qmomkit/conditioning_path.hpp#L183) |
| `ConditioningPathTable::levels` | [188](../../../include/qmomkit/conditioning_path.hpp#L188) |
| `ConditioningPathTable::empty` | [192](../../../include/qmomkit/conditioning_path.hpp#L192) |
| `ConditioningPathTable::operator[]` | [196](../../../include/qmomkit/conditioning_path.hpp#L196) |
| `ConditioningPathTable::begin` | [202](../../../include/qmomkit/conditioning_path.hpp#L202) |
| `ConditioningPathTable::end` | [206](../../../include/qmomkit/conditioning_path.hpp#L206) |
| `ConditioningPathTable::node_indices` | [212](../../../include/qmomkit/conditioning_path.hpp#L212) |

### `include/qmomkit/configuration.hpp`

[Source](../../../include/qmomkit/configuration.hpp) · [Responsibilities,
contracts, and evidence](initialization_and_configuration.md)

| Function or overload family | Definition lines |
|---|---|
| `ConfigurationResult` | [27](../../../include/qmomkit/configuration.hpp#L27) |
| `failure` | [29](../../../include/qmomkit/configuration.hpp#L29) |
| `ok` | [36](../../../include/qmomkit/configuration.hpp#L36) |
| `operator bool` | [41](../../../include/qmomkit/configuration.hpp#L41) |
| `value` | [46](../../../include/qmomkit/configuration.hpp#L46), [56](../../../include/qmomkit/configuration.hpp#L56) |
| `ConfigurationResult::value` | [51](../../../include/qmomkit/configuration.hpp#L51) |
| `error` | [61](../../../include/qmomkit/configuration.hpp#L61) |

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
| `CqmomResult::ok` | [390](../../../include/qmomkit/cqmom.hpp#L390) |
| `CqmomStateCandidate::ok` | [295](../../../include/qmomkit/cqmom.hpp#L295) |
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
| `QmomMomentSourceUpdateResult::ok` | [289](../../../include/qmomkit/moment_source.hpp#L289) |
| `CqmomMomentSourceUpdateResult::ok` | [302](../../../include/qmomkit/moment_source.hpp#L302) |

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
| `QmomNodeAdvectionResult::ok` | [104](../../../include/qmomkit/node_advection.hpp#L104) |
| `CqmomNodeAdvectionResult::ok` | [104](../../../include/qmomkit/node_advection.hpp#L104) |

### `include/qmomkit/phase_space_evolution.hpp`

[Source](../../../include/qmomkit/phase_space_evolution.hpp) ·
[Responsibilities, contracts, and evidence](phase_space_evolution.md)

| Function or overload family | Definition lines |
|---|---|
| `QmomPhaseSpaceEvolutionResult::ok` | [144](../../../include/qmomkit/phase_space_evolution.hpp#L144) |
| `CqmomPhaseSpaceEvolutionResult::ok` | [166](../../../include/qmomkit/phase_space_evolution.hpp#L166) |

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
| `QmomResult::ok` | [43](../../../include/qmomkit/qmom.hpp#L43) |
| `QmomState::moment_values` | [74](../../../include/qmomkit/qmom.hpp#L74) |
| `QmomState::quadrature` | [79](../../../include/qmomkit/qmom.hpp#L79) |
| `QmomState::active_node_count` | [84](../../../include/qmomkit/qmom.hpp#L84) |
| `QmomState::maximum_node_count` | [89](../../../include/qmomkit/qmom.hpp#L89) |
| `QmomState::qmom_options` | [94](../../../include/qmomkit/qmom.hpp#L94) |
| `QmomState::abscissa_bounds` | [99](../../../include/qmomkit/qmom.hpp#L99) |
| `QmomStateCandidate::ok` | [155](../../../include/qmomkit/qmom.hpp#L155) |

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
| `SspRk3Result::ok` | [88](../../../include/qmomkit/runge_kutta.hpp#L88) |

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
| `accepts_unnumbered_realizability` | [82](../../../src/configuration.cpp#L82) |
| `accepts_coordinate_realizability` | [87](../../../src/configuration.cpp#L87) |
| `trim` | [105](../../../src/configuration.cpp#L105) |
| `add_error` | [116](../../../src/configuration.cpp#L116) |
| `parse_size` | [122](../../../src/configuration.cpp#L122) |
| `parse_index` | [138](../../../src/configuration.cpp#L138) |
| `parse_coordinate_realizability_section` | [143](../../../src/configuration.cpp#L143) |
| `is_coordinate_realizability_key` | [152](../../../src/configuration.cpp#L152) |
| `is_configured_coordinate_realizability_key` | [167](../../../src/configuration.cpp#L167) |
| `is_decimal_number` | [180](../../../src/configuration.cpp#L180) |
| `parse_double` | [231](../../../src/configuration.cpp#L231) |
| `parse_size_list` | [247](../../../src/configuration.cpp#L247) |
| `parse_double_list` | [263](../../../src/configuration.cpp#L263) |
| `parse_optional_double_list` | [279](../../../src/configuration.cpp#L279) |
| `validate_names` | [301](../../../src/configuration.cpp#L301) |
| `parse_moment_key` | [322](../../../src/configuration.cpp#L322) |
| `mixture_component_prefix` | [343](../../../src/configuration.cpp#L343) |
| `is_mixture_component_section` | [348](../../../src/configuration.cpp#L348) |
| `is_independent_product_dimension_section` | [359](../../../src/configuration.cpp#L359) |
| `parse_file` | [369](../../../src/configuration.cpp#L369) |
| `require_entry` | [452](../../../src/configuration.cpp#L452) |
| `find_entry` | [463](../../../src/configuration.cpp#L463) |
| `reject_unknown_keys` | [473](../../../src/configuration.cpp#L473) |
| `reject_unknown_sections` | [492](../../../src/configuration.cpp#L492) |
| `validate_format_version` | [512](../../../src/configuration.cpp#L512) |
| `parse_phase_space_node_counts` | [521](../../../src/configuration.cpp#L521) |
| `validate_phase_space_names` | [533](../../../src/configuration.cpp#L533) |
| `validate_qmom_phase_space` | [543](../../../src/configuration.cpp#L543) |
| `validate_cqmom_phase_space` | [556](../../../src/configuration.cpp#L556) |
| `validate_optional_qmom_phase_space` | [576](../../../src/configuration.cpp#L576) |
| `validate_optional_cqmom_phase_space` | [592](../../../src/configuration.cpp#L592) |
| `parse_required_double` | [608](../../../src/configuration.cpp#L608) |
| `parse_optional_double` | [623](../../../src/configuration.cpp#L623) |
| `parse_optional_bool` | [638](../../../src/configuration.cpp#L638) |
| `parse_optional_size` | [659](../../../src/configuration.cpp#L659) |
| `parse_optional_positive_double` | [674](../../../src/configuration.cpp#L674) |
| `infer_qmom_initial_maximum_order` | [691](../../../src/configuration.cpp#L691) |
| `parse_optional_positive_size` | [718](../../../src/configuration.cpp#L718) |
| `distribution_parameters` | [736](../../../src/configuration.cpp#L736) |
| `is_distribution_parameter` | [743](../../../src/configuration.cpp#L743) |
| `reject_other_family_parameters` | [751](../../../src/configuration.cpp#L751) |
| `parse_distribution` | [766](../../../src/configuration.cpp#L766), [860](../../../src/configuration.cpp#L860) |
| `component_key_allowed` | [865](../../../src/configuration.cpp#L865) |
| `marginal_key_allowed` | [875](../../../src/configuration.cpp#L875) |
| `has_explicit_qmom_moments` | [889](../../../src/configuration.cpp#L889) |
| `has_explicit_cqmom_moments` | [894](../../../src/configuration.cpp#L894) |
| `qmom_allowed_keys` | [899](../../../src/configuration.cpp#L899) |
| `append_time_integration_allowed_keys` | [911](../../../src/configuration.cpp#L911) |
| `append_output_allowed_keys` | [933](../../../src/configuration.cpp#L933) |
| `with_runtime_allowed_keys` | [942](../../../src/configuration.cpp#L942) |
| `time_integration_sections` | [949](../../../src/configuration.cpp#L949) |
| `reject_unknown_qmom_sections` | [959](../../../src/configuration.cpp#L959) |
| `reject_unknown_qmom_keys` | [972](../../../src/configuration.cpp#L972) |
| `reject_unconsumed_indexed_entries` | [1006](../../../src/configuration.cpp#L1006) |
| `parse_mixture` | [1029](../../../src/configuration.cpp#L1029) |
| `has_mixture_component_sections` | [1111](../../../src/configuration.cpp#L1111) |
| `cqmom_explicit_allowed_keys` | [1123](../../../src/configuration.cpp#L1123) |
| `cqmom_initial_allowed_keys` | [1131](../../../src/configuration.cpp#L1131) |
| `reject_unknown_cqmom_sections` | [1140](../../../src/configuration.cpp#L1140) |
| `reject_unknown_cqmom_keys` | [1153](../../../src/configuration.cpp#L1153) |
| `reject_duplicate_indexed_section` | [991](../../../src/configuration.cpp#L991) |
| `parse_independent_product` | [1172](../../../src/configuration.cpp#L1172) |
| `checked_product` | [1266](../../../src/configuration.cpp#L1266) |
| `parse_explicit_tensor` | [1277](../../../src/configuration.cpp#L1277) |
| `format_configuration_error` | [1346](../../../src/configuration.cpp#L1346) |
| `format_configuration_errors` | [1362](../../../src/configuration.cpp#L1362) |
| `failure_result` | [1371](../../../src/configuration.cpp#L1371), [1376](../../../src/configuration.cpp#L1376) |
| `file_failure` | [1381](../../../src/configuration.cpp#L1381) |
| `validate_replacement_tensor_structure` | [1386](../../../src/configuration.cpp#L1386) |
| `format_moment_index` | [1409](../../../src/configuration.cpp#L1409) |
| `read_configuration_file` | [1425](../../../src/configuration.cpp#L1425) |
| `has_realizability_section` | [1447](../../../src/configuration.cpp#L1447) |
| `has_unnumbered_realizability_section` | [1457](../../../src/configuration.cpp#L1457) |
| `coordinate_realizability_sections` | [1467](../../../src/configuration.cpp#L1467) |
| `is_realizability_key` | [1481](../../../src/configuration.cpp#L1481) |
| `realizability_policy_fields` | [1486](../../../src/configuration.cpp#L1486) |
| `is_node_advection_key` | [1492](../../../src/configuration.cpp#L1492) |
| `node_advection_allowed_keys` | [1497](../../../src/configuration.cpp#L1497) |
| `parse_qmom_options_from_parsed` | [1506](../../../src/configuration.cpp#L1506) |
| `reject_unknown_qmom_option_keys` | [1534](../../../src/configuration.cpp#L1534) |
| `parse_realizability_check_mode_value` | [1552](../../../src/configuration.cpp#L1552) |
| `parse_realizability_check_mode` | [1571](../../../src/configuration.cpp#L1571) |
| `parse_moment_support` | [1582](../../../src/configuration.cpp#L1582) |
| `parse_phase_space_evolution_scheme_value` | [1605](../../../src/configuration.cpp#L1605) |
| `parse_node_advection_scheme_value` | [1624](../../../src/configuration.cpp#L1624) |
| `parse_node_advection_options_from_parsed` | [1643](../../../src/configuration.cpp#L1643) |
| `parse_realizability_policy_from_parsed` | [1691](../../../src/configuration.cpp#L1691) |
| `parse_coordinate_realizability_policies_from_parsed` | [1783](../../../src/configuration.cpp#L1783) |
| `parse_cqmom_diagnostic_level_value` | [1843](../../../src/configuration.cpp#L1843) |
| `parse_cqmom_diagnostic_level` | [1861](../../../src/configuration.cpp#L1861) |
| `parse_vandermonde_solver_value` | [1871](../../../src/configuration.cpp#L1871) |
| `parse_vandermonde_solver` | [1920](../../../src/configuration.cpp#L1920) |
| `has_output_section` | [1930](../../../src/configuration.cpp#L1930) |
| `is_output_key` | [1940](../../../src/configuration.cpp#L1940) |
| `validate_optional_phase_space_evolution_scheme` | [1945](../../../src/configuration.cpp#L1945) |
| `validate_optional_node_advection_options` | [1954](../../../src/configuration.cpp#L1954) |
| `validate_optional_realizability_policy` | [1959](../../../src/configuration.cpp#L1959) |
| `validate_cqmom_coordinate_realizability_ranks` | [1975](../../../src/configuration.cpp#L1975) |
| `reject_unknown_phase_space_evolution_sections` | [1989](../../../src/configuration.cpp#L1989) |
| `reject_unknown_phase_space_evolution_keys` | [2005](../../../src/configuration.cpp#L2005) |
| `parse_phase_space_evolution_scheme_from_parsed` | [2026](../../../src/configuration.cpp#L2026) |
| `validate_optional_cqmom_diagnostic_level` | [2046](../../../src/configuration.cpp#L2046) |
| `validate_optional_vandermonde_solver` | [2051](../../../src/configuration.cpp#L2051) |
| `parse_coordinate_support` | [2057](../../../src/configuration.cpp#L2057) |
| `parse_cqmom_case_settings` | [2120](../../../src/configuration.cpp#L2120) |
| `parse_output_options_from_parsed` | [2210](../../../src/configuration.cpp#L2210) |
| `parse_qmom_initial_condition_from_parsed` | [2308](../../../src/configuration.cpp#L2308) |
| `parse_cqmom_initial_condition_from_parsed` | [2406](../../../src/configuration.cpp#L2406) |
| `initialize_qmom_moments` | [2529](../../../src/configuration.cpp#L2529) |
| `initialize_cqmom_moment_tensor` | [2573](../../../src/configuration.cpp#L2573) |
| `initialize_cqmom_moments` | [2612](../../../src/configuration.cpp#L2612) |
| `parse_qmom_initial_condition` | [2659](../../../src/configuration.cpp#L2659) |
| `parse_cqmom_initial_condition` | [2665](../../../src/configuration.cpp#L2665) |
| `create_qmom_case` | [2671](../../../src/configuration.cpp#L2671) |
| `create_cqmom_case` | [2698](../../../src/configuration.cpp#L2698) |
| `with_initial_moments` | [2710](../../../src/configuration.cpp#L2710), [2739](../../../src/configuration.cpp#L2739), [2781](../../../src/configuration.cpp#L2781), [2788](../../../src/configuration.cpp#L2788) |
| `with_prepared_metadata_of` | [2761](../../../src/configuration.cpp#L2761) |
| `parse_qmom` | [2795](../../../src/configuration.cpp#L2795) |
| `parse_cqmom` | [2845](../../../src/configuration.cpp#L2845) |
| `parse_cqmom_prepared` | [2870](../../../src/configuration.cpp#L2870) |
| `parse_realizability_policy` | [2898](../../../src/configuration.cpp#L2898) |
| `parse_coordinate_realizability_policies` | [2936](../../../src/configuration.cpp#L2936) |
| `parse_phase_space_evolution_scheme` | [2956](../../../src/configuration.cpp#L2956) |
| `parse_node_advection_options` | [2975](../../../src/configuration.cpp#L2975) |
| `parse_qmom_phase_space_evolution_options` | [2994](../../../src/configuration.cpp#L2994) |
| `parse_cqmom_phase_space_evolution_options` | [3023](../../../src/configuration.cpp#L3023) |
| `parse_output_options` | [3049](../../../src/configuration.cpp#L3049) |
| `load_configuration` | [3068](../../../src/configuration.cpp#L3068) |
| `load_qmom` | [3080](../../../src/configuration.cpp#L3080) |
| `load_cqmom` | [3085](../../../src/configuration.cpp#L3085) |
| `load_cqmom_prepared` | [3090](../../../src/configuration.cpp#L3090) |
| `load_output_options` | [3095](../../../src/configuration.cpp#L3095) |
| `load_node_advection_options` | [3100](../../../src/configuration.cpp#L3100) |
| `load_realizability_policy` | [3105](../../../src/configuration.cpp#L3105) |
| `load_qmom_phase_space_evolution_options` | [3110](../../../src/configuration.cpp#L3110) |
| `load_cqmom_phase_space_evolution_options` | [3115](../../../src/configuration.cpp#L3115) |
| `parse_vandermonde_rcond` | [1902](../../../src/configuration.cpp#L1902) |

### `src/coordinate_scales.cpp`

[Source](../../../src/coordinate_scales.cpp) · [Responsibilities,
contracts, and evidence](qmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `finite_positive` | [11](../../../src/coordinate_scales.cpp#L11) |
| `convert_to_double` | [16](../../../src/coordinate_scales.cpp#L16) |
| `coordinate_scale` | [27](../../../src/coordinate_scales.cpp#L27) |
| `cqmom_coordinate_scales` | [52](../../../src/coordinate_scales.cpp#L52) |

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
| `ScaledProduct::ScaledProduct` | [214](../../../src/distributions.cpp#L214) |
| `ScaledProduct::multiply` | [221](../../../src/distributions.cpp#L221) |
| `ScaledProduct::value` | [232](../../../src/distributions.cpp#L232) |
| `multiply_beta_step` | [258](../../../src/distributions.cpp#L258) |
| `multiply_gamma_step` | [264](../../../src/distributions.cpp#L264) |
| `multiply_chi_squared_step` | [270](../../../src/distributions.cpp#L270) |
| `evaluate_product_moment` | [278](../../../src/distributions.cpp#L278) |
| `half_step_gamma_ratio` | [295](../../../src/distributions.cpp#L295) |
| `chi_first_moment` | [317](../../../src/distributions.cpp#L317) |
| `multiply_chi_step` | [323](../../../src/distributions.cpp#L323) |
| `evaluate_chi_moment` | [328](../../../src/distributions.cpp#L328) |
| `evaluate_weibull_moment` | [343](../../../src/distributions.cpp#L343) |
| `compute_from_products` | [353](../../../src/distributions.cpp#L353) |
| `compute_from_log_moments` | [377](../../../src/distributions.cpp#L377) |
| `check_selected_common_inputs` | [398](../../../src/distributions.cpp#L398) |
| `finish_selected_orders` | [430](../../../src/distributions.cpp#L430) |
| `zero_selected_moments` | [435](../../../src/distributions.cpp#L435) |
| `selected_failure_from_moment` | [440](../../../src/distributions.cpp#L440) |
| `NormalExpectationsResult::ok` | [452](../../../src/distributions.cpp#L452) |
| `normal_expectations` | [458](../../../src/distributions.cpp#L458) |
| `compute_normal_selected` | [502](../../../src/distributions.cpp#L502) |
| `compute_log_normal_selected` | [540](../../../src/distributions.cpp#L540) |
| `compute_beta_selected` | [569](../../../src/distributions.cpp#L569) |
| `compute_gamma_selected` | [596](../../../src/distributions.cpp#L596) |
| `compute_chi_selected` | [623](../../../src/distributions.cpp#L623) |
| `compute_chi_squared_selected` | [650](../../../src/distributions.cpp#L650) |
| `compute_weibull_selected` | [677](../../../src/distributions.cpp#L677) |
| `compute_normal` | [704](../../../src/distributions.cpp#L704) |
| `compute_log_normal` | [737](../../../src/distributions.cpp#L737) |
| `compute_beta` | [750](../../../src/distributions.cpp#L750) |
| `compute_gamma` | [761](../../../src/distributions.cpp#L761) |
| `compute_chi` | [772](../../../src/distributions.cpp#L772) |
| `compute_chi_squared` | [808](../../../src/distributions.cpp#L808) |
| `compute_weibull` | [820](../../../src/distributions.cpp#L820) |
| `RawMomentVisitor::operator()` | [837](../../../src/distributions.cpp#L837), [842](../../../src/distributions.cpp#L842), [847](../../../src/distributions.cpp#L847), [852](../../../src/distributions.cpp#L852), [857](../../../src/distributions.cpp#L857), [862](../../../src/distributions.cpp#L862), [867](../../../src/distributions.cpp#L867) |
| `SelectedRawMomentVisitor::operator()` | [880](../../../src/distributions.cpp#L880), [885](../../../src/distributions.cpp#L885), [890](../../../src/distributions.cpp#L890), [895](../../../src/distributions.cpp#L895), [900](../../../src/distributions.cpp#L900), [905](../../../src/distributions.cpp#L905), [910](../../../src/distributions.cpp#L910) |
| `add_compensated` | [916](../../../src/distributions.cpp#L916) |
| `add_selected_compensated` | [947](../../../src/distributions.cpp#L947) |
| `contains_order` | [981](../../../src/distributions.cpp#L981) |
| `order_position` | [993](../../../src/distributions.cpp#L993) |
| `compute_independent_product_marginal_moments` | [1005](../../../src/distributions.cpp#L1005), [1015](../../../src/distributions.cpp#L1015) |
| `first_requested_index_with_order` | [1025](../../../src/distributions.cpp#L1025) |
| `make_selected_raw_moment_result` | [1039](../../../src/distributions.cpp#L1039) |
| `compute_raw_moments` | [1083](../../../src/distributions.cpp#L1083), [1099](../../../src/distributions.cpp#L1099), [1109](../../../src/distributions.cpp#L1109), [1179](../../../src/distributions.cpp#L1179), [1260](../../../src/distributions.cpp#L1260), [1363](../../../src/distributions.cpp#L1363) |

### `src/internal/bjorck_pereyra.cpp`

[Source](../../../src/internal/bjorck_pereyra.cpp) · [Responsibilities,
contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_all_finite` | [16](../../../src/internal/bjorck_pereyra.cpp#L16) |
| `failure_result` | [31](../../../src/internal/bjorck_pereyra.cpp#L31) |
| `is_singular` | [46](../../../src/internal/bjorck_pereyra.cpp#L46) |
| `solve_bjorck_pereyra` | [80](../../../src/internal/bjorck_pereyra.cpp#L80) |

### `src/internal/cqmom_case.hpp`

[Source](../../../src/internal/cqmom_case.hpp) · [Responsibilities,
contracts, and evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `with_moments` | [32](../../../src/internal/cqmom_case.hpp#L32) |
| `extract_coordinate_marginal` | [40](../../../src/internal/cqmom_case.hpp#L40) |

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
| `validate_cqmom_quadrature_storage` | [131](../../../src/internal/cqmom_validation.hpp#L131) |
| `validate_cqmom_moment_rate` | [161](../../../src/internal/cqmom_validation.hpp#L161) |
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
| `make_identity` | [126](../../../src/internal/linear_algebra.cpp#L126) |
| `residual_norm` | [136](../../../src/internal/linear_algebra.cpp#L136) |
| `condition_estimate` | [156](../../../src/internal/linear_algebra.cpp#L156) |
| `power_of_two_shift` | [172](../../../src/internal/linear_algebra.cpp#L172) |
| `rotate` | [196](../../../src/internal/linear_algebra.cpp#L196) |
| `symmetric_jacobi_eigen` | [241](../../../src/internal/linear_algebra.cpp#L241) |
| `matrix_all_finite` | [16](../../../src/internal/linear_algebra.cpp#L16) |

### `src/internal/linear_algebra.hpp`

[Source](../../../src/internal/linear_algebra.hpp) · [Responsibilities,
contracts, and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `DenseMatrix::DenseMatrix` | [19](../../../src/internal/linear_algebra.hpp#L19) |
| `DenseMatrix::rows` | [21](../../../src/internal/linear_algebra.hpp#L21) |
| `DenseMatrix::cols` | [26](../../../src/internal/linear_algebra.hpp#L26) |
| `DenseMatrix::operator()` | [31](../../../src/internal/linear_algebra.hpp#L31), [40](../../../src/internal/linear_algebra.hpp#L40) |
| `DenseMatrix::storage_size` | [50](../../../src/internal/linear_algebra.hpp#L50) |

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
| `RecursiveQuadratureResult::ok` | [237](../../../src/internal/recursive_quadrature.hpp#L237) |

### `src/internal/vandermonde.cpp`

[Source](../../../src/internal/vandermonde.cpp) · [Responsibilities, contracts,
and evidence](conditional_systems.md)

| Function or overload family | Definition lines |
|---|---|
| `matrix_all_finite` | [20](../../../src/internal/vandermonde.cpp#L20) |
| `failure_result` | [35](../../../src/internal/vandermonde.cpp#L35) |
| `select_automatic_solver` | [40](../../../src/internal/vandermonde.cpp#L40) |
| `resolve_solver` | [45](../../../src/internal/vandermonde.cpp#L45) |
| `vandermonde_result` | [61](../../../src/internal/vandermonde.cpp#L61) |
| `qr_linear_solve_result` | [68](../../../src/internal/vandermonde.cpp#L68) |
| `svd_linear_solve_result` | [85](../../../src/internal/vandermonde.cpp#L85) |
| `residual_is_accepted` | [102](../../../src/internal/vandermonde.cpp#L102) |
| `weighted_vandermonde_input_failure` | [117](../../../src/internal/vandermonde.cpp#L117) |
| `make_scaled_vandermonde_system` | [168](../../../src/internal/vandermonde.cpp#L168) |
| `divide_by_weights` | [204](../../../src/internal/vandermonde.cpp#L204) |
| `make_conditional_vandermonde_matrix` | [226](../../../src/internal/vandermonde.cpp#L226) |
| `solve_conditional_vandermonde_system` | [240](../../../src/internal/vandermonde.cpp#L240) |

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
| `status_from_ssp_rk3` | [25](../../../src/moment_source.cpp#L25) |
| `message_from_status` | [39](../../../src/moment_source.cpp#L39) |
| `diagnostics_from_ssp_rk3` | [57](../../../src/moment_source.cpp#L57) |
| `invalid_result` | [68](../../../src/moment_source.cpp#L68) |
| `invalid_moment_data_result` | [75](../../../src/moment_source.cpp#L75) |
| `failed_moment_data_result` | [82](../../../src/moment_source.cpp#L82) |
| `multi_index_precedes` | [89](../../../src/moment_source.cpp#L89) |
| `validate_source_context` | [94](../../../src/moment_source.cpp#L94) |
| `validate_source_generalized_order` | [109](../../../src/moment_source.cpp#L109) |
| `validate_fixed_cqmom_moment_data` | [121](../../../src/moment_source.cpp#L121) |
| `validate_rate` | [128](../../../src/moment_source.cpp#L128) |
| `validate_source_rate` | [135](../../../src/moment_source.cpp#L135) |
| `validate_boundary_rate` | [147](../../../src/moment_source.cpp#L147) |
| `validate_cqmom_boundary_rate` | [159](../../../src/moment_source.cpp#L159) |
| `validate_cqmom_stage_state` | [166](../../../src/moment_source.cpp#L166) |
| `has_coordinate_realizability_policy` | [178](../../../src/moment_source.cpp#L178) |
| `accepts_coordinate_realizability_policies` | [185](../../../src/moment_source.cpp#L185) |
| `boundary_callback_error` | [206](../../../src/moment_source.cpp#L206) |
| `cqmom_boundary_callback_error` | [211](../../../src/moment_source.cpp#L211) |
| `validate_request` | [216](../../../src/moment_source.cpp#L216), [238](../../../src/moment_source.cpp#L238) |
| `validate_cqmom_source_update_controls` | [248](../../../src/moment_source.cpp#L248) |
| `validate_moment_data_request` | [264](../../../src/moment_source.cpp#L264), [281](../../../src/moment_source.cpp#L281) |
| `validate_trusted_moment_data_request` | [272](../../../src/moment_source.cpp#L272), [291](../../../src/moment_source.cpp#L291) |
| `compose_source_and_boundary` | [301](../../../src/moment_source.cpp#L301) |
| `advance_qmom_moments_with_source` | [347](../../../src/moment_source.cpp#L347), [858](../../../src/moment_source.cpp#L858), [878](../../../src/moment_source.cpp#L878) |
| `compose_cqmom_source_and_boundary` | [400](../../../src/moment_source.cpp#L400) |
| `advance_validated_cqmom_moments_with_source` | [427](../../../src/moment_source.cpp#L427) |
| `validate_qmom_moment_source_options` | [491](../../../src/moment_source.cpp#L491) |
| `validate_cqmom_moment_source_options` | [506](../../../src/moment_source.cpp#L506) |
| `QmomSourceContextAccess::create` | [530](../../../src/moment_source.cpp#L530) |
| `CqmomSourceContextAccess::create` | [536](../../../src/moment_source.cpp#L536) |
| `QmomSourceContext::QmomSourceContext` | [546](../../../src/moment_source.cpp#L546), [554](../../../src/moment_source.cpp#L554) |
| `QmomSourceContext::time` | [567](../../../src/moment_source.cpp#L567) |
| `QmomSourceContext::moment_count` | [572](../../../src/moment_source.cpp#L572) |
| `QmomSourceContext::size` | [577](../../../src/moment_source.cpp#L577) |
| `QmomSourceContext::values` | [582](../../../src/moment_source.cpp#L582) |
| `QmomSourceContext::value` | [587](../../../src/moment_source.cpp#L587) |
| `QmomSourceContext::moment` | [592](../../../src/moment_source.cpp#L592) |
| `QmomSourceContext::generalized_moment` | [609](../../../src/moment_source.cpp#L609) |
| `QmomSourceContext::for_each_node` | [615](../../../src/moment_source.cpp#L615) |
| `QmomSourceContext::reconstructed_stage_quadrature` | [629](../../../src/moment_source.cpp#L629) |
| `CqmomSourceContext::CqmomSourceContext` | [647](../../../src/moment_source.cpp#L647), [658](../../../src/moment_source.cpp#L658) |
| `CqmomSourceContext::time` | [672](../../../src/moment_source.cpp#L672) |
| `CqmomSourceContext::moment_count` | [677](../../../src/moment_source.cpp#L677) |
| `CqmomSourceContext::dimension_count` | [682](../../../src/moment_source.cpp#L682) |
| `CqmomSourceContext::size` | [687](../../../src/moment_source.cpp#L687) |
| `CqmomSourceContext::indices` | [692](../../../src/moment_source.cpp#L692) |
| `CqmomSourceContext::values` | [697](../../../src/moment_source.cpp#L697) |
| `CqmomSourceContext::index` | [702](../../../src/moment_source.cpp#L702) |
| `CqmomSourceContext::value` | [707](../../../src/moment_source.cpp#L707), [762](../../../src/moment_source.cpp#L762) |
| `CqmomSourceContext::sorted_index_positions` | [712](../../../src/moment_source.cpp#L712) |
| `CqmomSourceContext::position` | [749](../../../src/moment_source.cpp#L749) |
| `CqmomSourceContext::moment` | [772](../../../src/moment_source.cpp#L772) |
| `CqmomSourceContext::generalized_moment` | [788](../../../src/moment_source.cpp#L788) |
| `CqmomSourceContext::for_each_terminal_node` | [801](../../../src/moment_source.cpp#L801) |
| `CqmomSourceContext::reconstructed_stage_quadrature` | [822](../../../src/moment_source.cpp#L822) |
| `CqmomSourceContext::hierarchy` | [853](../../../src/moment_source.cpp#L853) |
| `advance_cqmom_moments_with_source` | [894](../../../src/moment_source.cpp#L894), [909](../../../src/moment_source.cpp#L909) |
| `advance_trusted_qmom_moments_with_source` | [928](../../../src/moment_source.cpp#L928), [953](../../../src/moment_source.cpp#L953) |
| `advance_trusted_cqmom_moments_with_source` | [980](../../../src/moment_source.cpp#L980), [995](../../../src/moment_source.cpp#L995), [1014](../../../src/moment_source.cpp#L1014), [1032](../../../src/moment_source.cpp#L1032) |
| `assign_source_rates` | [342](../../../src/moment_source.cpp#L342) |

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
| `integrate_nodes_ssp_rk3` | [217](../../../src/node_advection.cpp#L217) |
| `integrate_nodes` | [288](../../../src/node_advection.cpp#L288) |
| `scheme_name` | [319](../../../src/node_advection.cpp#L319) |
| `PooledCoordinateBuffer::PooledCoordinateBuffer` | [349](../../../src/node_advection.cpp#L349) |
| `PooledCoordinateBuffer::~PooledCoordinateBuffer` | [370](../../../src/node_advection.cpp#L370) |
| `PooledCoordinateBuffer::free_list` | [382](../../../src/node_advection.cpp#L382) |
| `make_phase_space_velocity` | [476](../../../src/node_advection.cpp#L476) |
| `advance_nodes` | [393](../../../src/node_advection.cpp#L393), [421](../../../src/node_advection.cpp#L421), [436](../../../src/node_advection.cpp#L436) |
| `format_node_advection_options` | [516](../../../src/node_advection.cpp#L516) |

### `src/output.cpp`

[Source](../../../src/output.cpp) · [Responsibilities, contracts, and
evidence](utilities.md)

| Function or overload family | Definition lines |
|---|---|
| `open_output_file` | [23](../../../src/output.cpp#L23) |
| `flush_or_throw` | [37](../../../src/output.cpp#L37) |
| `qmom_label` | [47](../../../src/output.cpp#L47) |
| `column_width_for` | [64](../../../src/output.cpp#L64) |
| `write_column_header` | [69](../../../src/output.cpp#L69) |
| `cqmom_label` | [82](../../../src/output.cpp#L82) |
| `sequential_indices` | [137](../../../src/output.cpp#L137) |
| `qmom_labels` | [153](../../../src/output.cpp#L153) |
| `cqmom_labels` | [169](../../../src/output.cpp#L169) |
| `flatten_selected_indices` | [185](../../../src/output.cpp#L185) |
| `all_cqmom_labels` | [196](../../../src/output.cpp#L196) |
| `cqmom_abscissa_label` | [201](../../../src/output.cpp#L201) |
| `qmom_quadrature_labels` | [206](../../../src/output.cpp#L206) |
| `cqmom_quadrature_labels` | [216](../../../src/output.cpp#L216) |
| `MomentWriter::for_qmom` | [240](../../../src/output.cpp#L240), [246](../../../src/output.cpp#L246), [252](../../../src/output.cpp#L252), [262](../../../src/output.cpp#L262) |
| `MomentWriter::for_cqmom` | [273](../../../src/output.cpp#L273), [280](../../../src/output.cpp#L280), [286](../../../src/output.cpp#L286), [297](../../../src/output.cpp#L297), [302](../../../src/output.cpp#L302), [311](../../../src/output.cpp#L311), [321](../../../src/output.cpp#L321) |
| `MomentWriter::MomentWriter` | [326](../../../src/output.cpp#L326), [344](../../../src/output.cpp#L344), [362](../../../src/output.cpp#L362), [379](../../../src/output.cpp#L379) |
| `MomentWriter::layout_mismatch_message` | [395](../../../src/output.cpp#L395) |
| `MomentWriter::write_header` | [421](../../../src/output.cpp#L421) |
| `MomentWriter::write_row` | [426](../../../src/output.cpp#L426) |
| `MomentWriter::write_aligned_row` | [456](../../../src/output.cpp#L456) |
| `MomentWriter::write` | [484](../../../src/output.cpp#L484), [503](../../../src/output.cpp#L503), [526](../../../src/output.cpp#L526), [545](../../../src/output.cpp#L545), [564](../../../src/output.cpp#L564), [587](../../../src/output.cpp#L587), [615](../../../src/output.cpp#L615), [634](../../../src/output.cpp#L634), [648](../../../src/output.cpp#L648), [667](../../../src/output.cpp#L667) |
| `QuadratureWriter::for_qmom` | [681](../../../src/output.cpp#L681) |
| `QuadratureWriter::for_cqmom` | [686](../../../src/output.cpp#L686) |
| `QuadratureWriter::QuadratureWriter` | [691](../../../src/output.cpp#L691) |
| `QuadratureWriter::write_header` | [707](../../../src/output.cpp#L707) |
| `QuadratureWriter::write_row` | [712](../../../src/output.cpp#L712) |
| `QuadratureWriter::write_qmom_row` | [735](../../../src/output.cpp#L735) |
| `QuadratureWriter::write` | [774](../../../src/output.cpp#L774), [779](../../../src/output.cpp#L779), [784](../../../src/output.cpp#L784) |
| `centered_column_label` | [52](../../../src/output.cpp#L52) |
| `generalized_moment_label` | [99](../../../src/output.cpp#L99) |
| `append_generalized_labels` | [116](../../../src/output.cpp#L116) |
| `validate_generalized_moment_orders` | [125](../../../src/output.cpp#L125) |

### `src/phase_space_evolution.cpp`

[Source](../../../src/phase_space_evolution.cpp) · [Responsibilities,
contracts, and evidence](phase_space_evolution.md)

| Function or overload family | Definition lines |
|---|---|
| `apply_batch` | [49](../../../src/phase_space_evolution.cpp#L49) |
| `validate_time_input` | [62](../../../src/phase_space_evolution.cpp#L62) |
| `validate_phase_space_evolution_scheme` | [74](../../../src/phase_space_evolution.cpp#L74) |
| `execute_phase_space_evolution_sequence` | [87](../../../src/phase_space_evolution.cpp#L87) |
| `validate_qmom_input` | [122](../../../src/phase_space_evolution.cpp#L122) |
| `validate_node_advection_options` | [165](../../../src/phase_space_evolution.cpp#L165) |
| `validate_qmom_source_options` | [182](../../../src/phase_space_evolution.cpp#L182) |
| `validate_boundary` | [191](../../../src/phase_space_evolution.cpp#L191) |
| `validate_cqmom_boundary` | [199](../../../src/phase_space_evolution.cpp#L199) |
| `successful_reconstruction` | [207](../../../src/phase_space_evolution.cpp#L207) |
| `record_advection_diagnostics` | [212](../../../src/phase_space_evolution.cpp#L212) |
| `record_source_update_diagnostics` | [228](../../../src/phase_space_evolution.cpp#L228) |
| `qmom_from_quadrature` | [235](../../../src/phase_space_evolution.cpp#L235) |
| `add_advection_totals` | [253](../../../src/phase_space_evolution.cpp#L253) |
| `finalize_diagnostics` | [262](../../../src/phase_space_evolution.cpp#L262) |
| `record_failure_diagnostics` | [268](../../../src/phase_space_evolution.cpp#L268) |
| `qmom_failure` | [275](../../../src/phase_space_evolution.cpp#L275) |
| `select_post_source_qmom` | [283](../../../src/phase_space_evolution.cpp#L283) |
| `validate_qmom_advection_input` | [326](../../../src/phase_space_evolution.cpp#L326) |
| `validate_qmom_moment_source_input` | [339](../../../src/phase_space_evolution.cpp#L339) |
| `advance_qmom_advection_stage` | [375](../../../src/phase_space_evolution.cpp#L375) |
| `advance_qmom_moment_source_stage` | [425](../../../src/phase_space_evolution.cpp#L425) |
| `operator<` | [495](../../../src/phase_space_evolution.cpp#L495) |
| `same_active_branch_count` | [508](../../../src/phase_space_evolution.cpp#L508) |
| `active_branch_key` | [515](../../../src/phase_space_evolution.cpp#L515) |
| `path_has_prefix` | [520](../../../src/phase_space_evolution.cpp#L520) |
| `has_changed_strict_ancestor` | [536](../../../src/phase_space_evolution.cpp#L536) |
| `record_cqmom_node_count_transitions` | [554](../../../src/phase_space_evolution.cpp#L554) |
| `advance_qmom_phase_space_impl` | [652](../../../src/phase_space_evolution.cpp#L652) |
| `advance_qmom_advection` | [699](../../../src/phase_space_evolution.cpp#L699) |
| `advance_qmom_moment_source_impl` | [725](../../../src/phase_space_evolution.cpp#L725) |
| `advance_qmom_moment_source` | [753](../../../src/phase_space_evolution.cpp#L753), [759](../../../src/phase_space_evolution.cpp#L759) |
| `advance_qmom_phase_space` | [766](../../../src/phase_space_evolution.cpp#L766), [773](../../../src/phase_space_evolution.cpp#L773) |
| `initialize_cqmom_diagnostics` | [783](../../../src/phase_space_evolution.cpp#L783), [789](../../../src/phase_space_evolution.cpp#L789) |
| `validate_cqmom_moments_for_phase_space_evolution` | [794](../../../src/phase_space_evolution.cpp#L794) |
| `validate_cqmom_state_for_phase_space_evolution` | [801](../../../src/phase_space_evolution.cpp#L801) |
| `validate_cqmom_advection_input` | [806](../../../src/phase_space_evolution.cpp#L806) |
| `validate_cqmom_moment_source_input` | [816](../../../src/phase_space_evolution.cpp#L816) |
| `validate_cqmom_phase_space_evolution_input` | [828](../../../src/phase_space_evolution.cpp#L828) |
| `CqmomStageOutcome::ok` | [855](../../../src/phase_space_evolution.cpp#L855) |
| `take_cqmom_advection_active_counts` | [869](../../../src/phase_space_evolution.cpp#L869) |
| `cqmom_state_failure` | [878](../../../src/phase_space_evolution.cpp#L878) |
| `cqmom_success` | [886](../../../src/phase_space_evolution.cpp#L886) |
| `materialize_cqmom_phase_space_evolution_result` | [891](../../../src/phase_space_evolution.cpp#L891) |
| `cqmom_state_from_accepted_quadrature` | [911](../../../src/phase_space_evolution.cpp#L911) |
| `prepare_cqmom_phase_space_evolution_working_state` | [920](../../../src/phase_space_evolution.cpp#L920) |
| `advance_cqmom_advection_stage` | [956](../../../src/phase_space_evolution.cpp#L956) |
| `advance_cqmom_moment_source_stage` | [1015](../../../src/phase_space_evolution.cpp#L1015) |
| `advance_cqmom_advection_impl` | [1078](../../../src/phase_space_evolution.cpp#L1078) |
| `advance_cqmom_moment_source_impl` | [1085](../../../src/phase_space_evolution.cpp#L1085) |
| `advance_cqmom_phase_space_state_impl` | [1092](../../../src/phase_space_evolution.cpp#L1092) |
| `advance_cqmom_advection` | [1157](../../../src/phase_space_evolution.cpp#L1157) |
| `advance_cqmom_moment_source_public_impl` | [1180](../../../src/phase_space_evolution.cpp#L1180) |
| `advance_cqmom_moment_source` | [1204](../../../src/phase_space_evolution.cpp#L1204), [1210](../../../src/phase_space_evolution.cpp#L1210) |
| `advance_cqmom_phase_space` | [1217](../../../src/phase_space_evolution.cpp#L1217), [1224](../../../src/phase_space_evolution.cpp#L1224) |
| `advance_qmom_advection_batch` | [1232](../../../src/phase_space_evolution.cpp#L1232) |
| `advance_cqmom_advection_batch` | [1240](../../../src/phase_space_evolution.cpp#L1240) |
| `advance_qmom_moment_source_batch` | [1248](../../../src/phase_space_evolution.cpp#L1248), [1264](../../../src/phase_space_evolution.cpp#L1264) |
| `advance_cqmom_moment_source_batch` | [1256](../../../src/phase_space_evolution.cpp#L1256), [1273](../../../src/phase_space_evolution.cpp#L1273) |
| `advance_qmom_phase_space_batch` | [1282](../../../src/phase_space_evolution.cpp#L1282), [1300](../../../src/phase_space_evolution.cpp#L1300) |
| `advance_cqmom_phase_space_batch` | [1291](../../../src/phase_space_evolution.cpp#L1291), [1310](../../../src/phase_space_evolution.cpp#L1310) |
| `advance_qmom_advection_batch_with_index` | [1320](../../../src/phase_space_evolution.cpp#L1320) |
| `advance_cqmom_advection_batch_with_index` | [1337](../../../src/phase_space_evolution.cpp#L1337) |
| `advance_qmom_moment_source_batch_with_index` | [1354](../../../src/phase_space_evolution.cpp#L1354), [1384](../../../src/phase_space_evolution.cpp#L1384) |
| `advance_cqmom_moment_source_batch_with_index` | [1369](../../../src/phase_space_evolution.cpp#L1369), [1405](../../../src/phase_space_evolution.cpp#L1405) |
| `advance_qmom_phase_space_batch_with_index` | [1426](../../../src/phase_space_evolution.cpp#L1426), [1471](../../../src/phase_space_evolution.cpp#L1471) |
| `advance_cqmom_phase_space_batch_with_index` | [1448](../../../src/phase_space_evolution.cpp#L1448), [1499](../../../src/phase_space_evolution.cpp#L1499) |

### `src/public_api.cpp`

[Source](../../../src/public_api.cpp) · [Responsibilities, contracts, and
evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `to_wheeler_options` | [29](../../../src/public_api.cpp#L29) |
| `compute_qmom_from_moments` | [36](../../../src/public_api.cpp#L36) |
| `to_qmom_status` | [53](../../../src/public_api.cpp#L53) |
| `hierarchy` | [70](../../../src/public_api.cpp#L70) |
| `qmom_options` | [75](../../../src/public_api.cpp#L75) |
| `conditional_options` | [80](../../../src/public_api.cpp#L80) |
| `diagnostic_level` | [85](../../../src/public_api.cpp#L85) |
| `vandermonde_solver` | [90](../../../src/public_api.cpp#L90) |
| `coordinate_support` | [95](../../../src/public_api.cpp#L95) |
| `diagnostic_level_name` | [100](../../../src/public_api.cpp#L100) |
| `vandermonde_solver_name` | [114](../../../src/public_api.cpp#L114) |
| `compute_qmom` | [135](../../../src/public_api.cpp#L135), [140](../../../src/public_api.cpp#L140), [146](../../../src/public_api.cpp#L146) |
| `QmomState::QmomState` | [151](../../../src/public_api.cpp#L151) |
| `below_number_density_cutoff` | [162](../../../src/public_api.cpp#L162) |
| `policy_verdict` | [167](../../../src/public_api.cpp#L167) |
| `policy_accepts` | [177](../../../src/public_api.cpp#L177) |
| `qmom_state_candidate` | [187](../../../src/public_api.cpp#L187) |
| `checked_qmom_state_candidate` | [241](../../../src/public_api.cpp#L241) |
| `QmomCase::state_candidate` | [254](../../../src/public_api.cpp#L254), [259](../../../src/public_api.cpp#L259) |
| `QmomCase::initial_state` | [264](../../../src/public_api.cpp#L264) |
| `format_cqmom_settings` | [288](../../../src/public_api.cpp#L288) |
| `CqmomPreparedCase::CqmomPreparedCase` | [298](../../../src/public_api.cpp#L298), [300](../../../src/public_api.cpp#L300), [305](../../../src/public_api.cpp#L305) |
| `CqmomPreparedCase::operator=` | [309](../../../src/public_api.cpp#L309) |
| `CqmomPreparedCase::hierarchy` | [324](../../../src/public_api.cpp#L324) |
| `CqmomPreparedCase::qmom_options` | [329](../../../src/public_api.cpp#L329) |
| `CqmomPreparedCase::conditional_options` | [334](../../../src/public_api.cpp#L334) |
| `CqmomPreparedCase::diagnostic_level` | [339](../../../src/public_api.cpp#L339) |
| `CqmomPreparedCase::vandermonde_solver` | [344](../../../src/public_api.cpp#L344) |
| `CqmomPreparedCase::abscissa_bounds` | [349](../../../src/public_api.cpp#L349) |
| `cqmom_state_candidate` | [362](../../../src/public_api.cpp#L362) |
| `initial_cqmom_state` | [441](../../../src/public_api.cpp#L441) |
| `prepared_cqmom_state_candidate` | [454](../../../src/public_api.cpp#L454) |
| `CqmomCase::initial_state` | [477](../../../src/public_api.cpp#L477) |
| `CqmomPreparedCase::state_candidate` | [484](../../../src/public_api.cpp#L484), [489](../../../src/public_api.cpp#L489) |
| `CqmomPreparedCase::initial_state` | [495](../../../src/public_api.cpp#L495) |
| `CqmomState::CqmomState` | [500](../../../src/public_api.cpp#L500), [502](../../../src/public_api.cpp#L502) |
| `CqmomState::operator=` | [504](../../../src/public_api.cpp#L504) |
| `CqmomState::hierarchy` | [519](../../../src/public_api.cpp#L519) |
| `CqmomState::quadrature` | [524](../../../src/public_api.cpp#L524) |
| `CqmomState::qmom_options` | [533](../../../src/public_api.cpp#L533) |
| `CqmomState::conditional_options` | [538](../../../src/public_api.cpp#L538) |
| `CqmomState::diagnostic_level` | [543](../../../src/public_api.cpp#L543) |
| `CqmomState::vandermonde_solver` | [548](../../../src/public_api.cpp#L548) |
| `CqmomState::abscissa_bounds` | [553](../../../src/public_api.cpp#L553) |
| `CqmomState::active_branch_counts` | [558](../../../src/public_api.cpp#L558) |
| `CqmomState::moment_indices` | [563](../../../src/public_api.cpp#L563) |
| `CqmomState::moment_values` | [568](../../../src/public_api.cpp#L568) |
| `CqmomState::moment` | [573](../../../src/public_api.cpp#L573) |
| `validate_prepared_metadata` | [586](../../../src/public_api.cpp#L586) |
| `make_prepared_metadata` | [595](../../../src/public_api.cpp#L595) |
| `CqmomConstructionAccess::create_prepared_case` | [606](../../../src/public_api.cpp#L606) |
| `CqmomConstructionAccess::create_state` | [612](../../../src/public_api.cpp#L612), [617](../../../src/public_api.cpp#L617), [623](../../../src/public_api.cpp#L623), [631](../../../src/public_api.cpp#L631) |
| `CqmomConstructionAccess::case_definition` | [638](../../../src/public_api.cpp#L638) |
| `CqmomConstructionAccess::state_definition` | [643](../../../src/public_api.cpp#L643) |
| `CqmomConstructionAccess::state_active_branch_counts` | [648](../../../src/public_api.cpp#L648) |
| `CqmomConstructionAccess::state_accepted_quadrature` | [653](../../../src/public_api.cpp#L653) |
| `CqmomConstructionAccess::case_prepared_metadata` | [658](../../../src/public_api.cpp#L658) |
| `CqmomConstructionAccess::state_prepared_metadata` | [663](../../../src/public_api.cpp#L663) |
| `CqmomConstructionAccess::create_quadrature` | [668](../../../src/public_api.cpp#L668) |
| `CqmomQuadrature::CqmomQuadrature` | [678](../../../src/public_api.cpp#L678) |
| `CqmomQuadrature::weight` | [685](../../../src/public_api.cpp#L685) |
| `CqmomQuadrature::abscissa` | [694](../../../src/public_api.cpp#L694) |
| `CqmomQuadrature::node` | [707](../../../src/public_api.cpp#L707) |
| `to_cqmom_status` | [720](../../../src/public_api.cpp#L720) |
| `to_cqmom_diagnostics` | [736](../../../src/public_api.cpp#L736), [765](../../../src/public_api.cpp#L765) |
| `terminal_paths_for_public_quadrature` | [776](../../../src/public_api.cpp#L776) |
| `failed_phase_space_dimension` | [798](../../../src/public_api.cpp#L798) |
| `to_cqmom_result` | [808](../../../src/public_api.cpp#L808) |
| `compute_cqmom` | [828](../../../src/public_api.cpp#L828), [891](../../../src/public_api.cpp#L891) |
| `compute_cqmom_internal` | [852](../../../src/public_api.cpp#L852) |
| `invalid_cqmom_hierarchy_result` | [818](../../../src/public_api.cpp#L818) |

### `src/realizability.cpp`

[Source](../../../src/realizability.cpp) · [Responsibilities, contracts, and
evidence](realizability.md)

| Function or overload family | Definition lines |
|---|---|
| `RealizabilityResult::realizable` | [21](../../../src/realizability.cpp#L21) |
| `RealizabilityResult::on_boundary` | [26](../../../src/realizability.cpp#L26) |
| `support_name` | [83](../../../src/realizability.cpp#L83) |
| `make_result` | [97](../../../src/realizability.cpp#L97) |
| `vector_scale` | [103](../../../src/realizability.cpp#L103) |
| `classify_eigenvalue` | [131](../../../src/realizability.cpp#L131) |
| `hankel_matrix` | [148](../../../src/realizability.cpp#L148) |
| `hankel_uncertainty` | [166](../../../src/realizability.cpp#L166) |
| `analyze_matrix` | [179](../../../src/realizability.cpp#L179) |
| `update_limiting_matrix` | [265](../../../src/realizability.cpp#L265) |
| `incorporate_matrix` | [279](../../../src/realizability.cpp#L279) |
| `lower_localizing_sequence` | [308](../../../src/realizability.cpp#L308) |
| `upper_localizing_sequence` | [321](../../../src/realizability.cpp#L321) |
| `interval_localizing_sequence` | [336](../../../src/realizability.cpp#L336) |
| `extension_column` | [350](../../../src/realizability.cpp#L350) |
| `euclidean_norm` | [361](../../../src/realizability.cpp#L361) |
| `projection_on_eigenvector` | [371](../../../src/realizability.cpp#L371) |
| `kernel_projection_tolerance` | [387](../../../src/realizability.cpp#L387) |
| `check_range_consistency` | [405](../../../src/realizability.cpp#L405) |
| `incorporate_range_consistency` | [431](../../../src/realizability.cpp#L431) |
| `numerical_rank` | [456](../../../src/realizability.cpp#L456) |
| `singular_as_given` | [466](../../../src/realizability.cpp#L466) |
| `finish_classification` | [472](../../../src/realizability.cpp#L472) |
| `tolerance_error` | [488](../../../src/realizability.cpp#L488) |
| `interval_error` | [493](../../../src/realizability.cpp#L493) |
| `validate_common_inputs` | [500](../../../src/realizability.cpp#L500) |
| `check_hamburger` | [534](../../../src/realizability.cpp#L534) |
| `check_stieltjes` | [588](../../../src/realizability.cpp#L588) |
| `check_hausdorff` | [644](../../../src/realizability.cpp#L644) |
| `TestedMoments::bound` | [714](../../../src/realizability.cpp#L714) |
| `standardized_moments` | [720](../../../src/realizability.cpp#L720) |
| `with_support_context` | [794](../../../src/realizability.cpp#L794) |
| `check_realizability` | [805](../../../src/realizability.cpp#L805), [834](../../../src/realizability.cpp#L834) |
| `realizability_policy_error` | [852](../../../src/realizability.cpp#L852) |

### `src/recursive_quadrature.cpp`

[Source](../../../src/recursive_quadrature.cpp) · [Responsibilities, contracts,
and evidence](cqmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `create_ordered_root_layout` | [36](../../../src/recursive_quadrature.cpp#L36) |
| `append_unique_index` | [110](../../../src/recursive_quadrature.cpp#L110) |
| `local_branch_offset` | [130](../../../src/recursive_quadrature.cpp#L130) |
| `create_child_layout_case` | [159](../../../src/recursive_quadrature.cpp#L159) |
| `create_branch_moment_layout` | [207](../../../src/recursive_quadrature.cpp#L207) |
| `PreparedCqmomMetadata::PreparedCqmomMetadata` | [270](../../../src/recursive_quadrature.cpp#L270) |
| `create_prepared_cqmom_metadata` | [282](../../../src/recursive_quadrature.cpp#L282) |
| `prepared_cqmom_metadata_is_compatible` | [297](../../../src/recursive_quadrature.cpp#L297) |
| `extract_wheeler_moments` | [312](../../../src/recursive_quadrature.cpp#L312) |
| `child_layout_for_parent_node_count` | [331](../../../src/recursive_quadrature.cpp#L331) |
| `make_branch_rhs_matrix` | [345](../../../src/recursive_quadrature.cpp#L345) |
| `validate_positive_extents` | [374](../../../src/recursive_quadrature.cpp#L374) |
| `checked_multiply` | [386](../../../src/recursive_quadrature.cpp#L386) |
| `default_conditioning_order` | [396](../../../src/recursive_quadrature.cpp#L396) |
| `validate_conditioning_order` | [403](../../../src/recursive_quadrature.cpp#L403) |
| `ordered_extents` | [425](../../../src/recursive_quadrature.cpp#L425) |
| `from_ordered_index` | [436](../../../src/recursive_quadrature.cpp#L436) |
| `validate_abscissas` | [446](../../../src/recursive_quadrature.cpp#L446) |
| `valid_conditional_options` | [487](../../../src/recursive_quadrature.cpp#L487) |
| `valid_diagnostic_level` | [494](../../../src/recursive_quadrature.cpp#L494) |
| `records_success_branch_diagnostics` | [506](../../../src/recursive_quadrature.cpp#L506) |
| `create_ordered_root_slice` | [514](../../../src/recursive_quadrature.cpp#L514) |
| `one_dimensional_moments` | [522](../../../src/recursive_quadrature.cpp#L522) |
| `reconstruct_child_slices` | [527](../../../src/recursive_quadrature.cpp#L527) |
| `create_recursive_result` | [619](../../../src/recursive_quadrature.cpp#L619) |
| `record_failure` | [625](../../../src/recursive_quadrature.cpp#L625) |
| `validate_recursive_request_preflight` | [634](../../../src/recursive_quadrature.cpp#L634) |
| `path_has_prefix` | [656](../../../src/recursive_quadrature.cpp#L656) |
| `branch_limit_invalidated` | [672](../../../src/recursive_quadrature.cpp#L672) |
| `parent_path_less` | [684](../../../src/recursive_quadrature.cpp#L684) |
| `validate_active_branch_counts` | [694](../../../src/recursive_quadrature.cpp#L694) |
| `recursive_branch_node_count_limit` | [773](../../../src/recursive_quadrature.cpp#L773) |
| `compute_branch_quadrature` | [790](../../../src/recursive_quadrature.cpp#L790) |
| `validate_coordinate_support` | [796](../../../src/recursive_quadrature.cpp#L796) |
| `support_interval_for_dimension` | [823](../../../src/recursive_quadrature.cpp#L823) |
| `record_bounded_candidate_diagnostics` | [833](../../../src/recursive_quadrature.cpp#L833) |
| `branch_growth_trial_allowed` | [854](../../../src/recursive_quadrature.cpp#L854) |
| `product_extents` | [864](../../../src/recursive_quadrature.cpp#L864) |
| `flatten_multi_index` | [876](../../../src/recursive_quadrature.cpp#L876) |
| `unflatten_multi_index` | [896](../../../src/recursive_quadrature.cpp#L896) |
| `enumerate_multi_indices` | [1046](../../../src/recursive_quadrature.cpp#L1046) |
| `PhaseSpaceHierarchy::PhaseSpaceHierarchy` | [1058](../../../src/recursive_quadrature.cpp#L1058), [1065](../../../src/recursive_quadrature.cpp#L1065) |
| `PhaseSpaceHierarchy::node_count` | [1071](../../../src/recursive_quadrature.cpp#L1071) |
| `PhaseSpaceHierarchy::moment_extents` | [1080](../../../src/recursive_quadrature.cpp#L1080) |
| `PhaseSpaceHierarchy::flatten_node_index` | [1091](../../../src/recursive_quadrature.cpp#L1091) |
| `PhaseSpaceHierarchy::unflatten_node_index` | [1100](../../../src/recursive_quadrature.cpp#L1100) |
| `create_cqmom_moment_layout` | [1105](../../../src/recursive_quadrature.cpp#L1105) |
| `CqmomMomentData::CqmomMomentData` | [1167](../../../src/recursive_quadrature.cpp#L1167) |
| `CqmomMomentData::contains_required_index` | [1172](../../../src/recursive_quadrature.cpp#L1172) |
| `CqmomMomentData::moment` | [1177](../../../src/recursive_quadrature.cpp#L1177) |
| `extract_cqmom_moments` | [1193](../../../src/recursive_quadrature.cpp#L1193) |
| `create_cqmom_moment_data` | [1220](../../../src/recursive_quadrature.cpp#L1220) |
| `enumerate_quadrature_node_multi_indices` | [1242](../../../src/recursive_quadrature.cpp#L1242) |
| `rectangular_terminal_paths` | [1257](../../../src/recursive_quadrature.cpp#L1257) |
| `RecursiveQuadrature::RecursiveQuadrature` | [1292](../../../src/recursive_quadrature.cpp#L1292), [1299](../../../src/recursive_quadrature.cpp#L1299) |
| `RecursiveQuadrature::hierarchy` | [1310](../../../src/recursive_quadrature.cpp#L1310) |
| `RecursiveQuadrature::number_of_dimensions` | [1315](../../../src/recursive_quadrature.cpp#L1315) |
| `RecursiveQuadrature::total_nodes` | [1320](../../../src/recursive_quadrature.cpp#L1320) |
| `RecursiveQuadrature::weights` | [1325](../../../src/recursive_quadrature.cpp#L1325) |
| `RecursiveQuadrature::abscissas` | [1330](../../../src/recursive_quadrature.cpp#L1330) |
| `RecursiveQuadrature::adaptive_branch_counts_enabled` | [1335](../../../src/recursive_quadrature.cpp#L1335) |
| `RecursiveQuadrature::has_rectangular_storage` | [1340](../../../src/recursive_quadrature.cpp#L1340) |
| `RecursiveQuadrature::tree_nodes` | [1345](../../../src/recursive_quadrature.cpp#L1345) |
| `RecursiveQuadrature::terminal_node_ids` | [1350](../../../src/recursive_quadrature.cpp#L1350) |
| `RecursiveQuadrature::terminal_paths` | [1355](../../../src/recursive_quadrature.cpp#L1355) |
| `RecursiveQuadrature::weight` | [1360](../../../src/recursive_quadrature.cpp#L1360) |
| `RecursiveQuadrature::abscissa` | [1369](../../../src/recursive_quadrature.cpp#L1369) |
| `RecursiveQuadrature::node` | [1382](../../../src/recursive_quadrature.cpp#L1382) |
| `RecursiveQuadrature::set_weight` | [1392](../../../src/recursive_quadrature.cpp#L1392) |
| `RecursiveQuadrature::set_abscissa` | [1401](../../../src/recursive_quadrature.cpp#L1401) |
| `RecursiveQuadrature::set_adaptive_branch_counts_enabled` | [1414](../../../src/recursive_quadrature.cpp#L1414) |
| `RecursiveQuadrature::add_tree_node` | [1427](../../../src/recursive_quadrature.cpp#L1427) |
| `RecursiveQuadrature::set_terminal_quadrature` | [1439](../../../src/recursive_quadrature.cpp#L1439) |
| `compute_recursive_quadrature_from_validated_moment_data` | [1472](../../../src/recursive_quadrature.cpp#L1472) |
| `compute_recursive_quadrature_from_moment_data` | [1826](../../../src/recursive_quadrature.cpp#L1826) |
| `compute_recursive_quadrature_from_trusted_moment_data` | [1852](../../../src/recursive_quadrature.cpp#L1852) |
| `operator==` | [915](../../../src/recursive_quadrature.cpp#L915), [925](../../../src/recursive_quadrature.cpp#L925), [930](../../../src/recursive_quadrature.cpp#L930), [1036](../../../src/recursive_quadrature.cpp#L1036) |
| `operator!=` | [920](../../../src/recursive_quadrature.cpp#L920), [935](../../../src/recursive_quadrature.cpp#L935), [940](../../../src/recursive_quadrature.cpp#L940), [1041](../../../src/recursive_quadrature.cpp#L1041) |
| `MultiIndexTable::MultiIndexTable` | [945](../../../src/recursive_quadrature.cpp#L945), [961](../../../src/recursive_quadrature.cpp#L961) |
| `MultiIndexTable::at` | [979](../../../src/recursive_quadrature.cpp#L979) |
| `MultiIndexTable::push_back` | [988](../../../src/recursive_quadrature.cpp#L988), [993](../../../src/recursive_quadrature.cpp#L993) |
| `MultiIndexTable::pop_back` | [1006](../../../src/recursive_quadrature.cpp#L1006) |
| `MultiIndexTable::reserve` | [1015](../../../src/recursive_quadrature.cpp#L1015) |
| `MultiIndexTable::clear` | [1020](../../../src/recursive_quadrature.cpp#L1020) |
| `MultiIndexTable::to_multi_indices` | [1025](../../../src/recursive_quadrature.cpp#L1025) |

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
| `minimum_step` | [68](../../../src/runge_kutta.cpp#L68) |
| `error_norm` | [80](../../../src/runge_kutta.cpp#L80) |
| `invalid_embedded_rk_result` | [100](../../../src/runge_kutta.cpp#L100) |
| `stage_accepted` | [105](../../../src/runge_kutta.cpp#L105) |
| `try_ssp_rk3_step` | [110](../../../src/runge_kutta.cpp#L110) |
| `ssp_rk3_step` | [172](../../../src/runge_kutta.cpp#L172) |
| `adaptive_rk_options_error` | [219](../../../src/runge_kutta.cpp#L219) |
| `integrate_embedded_rk` | [246](../../../src/runge_kutta.cpp#L246) |

### `src/wheeler.cpp`

[Source](../../../src/wheeler.cpp) · [Responsibilities, contracts, and
evidence](qmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `one_node_result` | [42](../../../src/wheeler.cpp#L42) |
| `negative_mass_result` | [55](../../../src/wheeler.cpp#L55) |
| `compute_recurrence` | [63](../../../src/wheeler.cpp#L63) |
| `make_jacobi_matrix` | [93](../../../src/wheeler.cpp#L93) |
| `recurrence_is_finite` | [114](../../../src/wheeler.cpp#L114) |
| `has_negative_recurrence_b` | [119](../../../src/wheeler.cpp#L119) |
| `sigma_tolerances` | [133](../../../src/wheeler.cpp#L133) |
| `compute_weights` | [161](../../../src/wheeler.cpp#L161) |
| `adaptive_accepts` | [173](../../../src/wheeler.cpp#L173) |
| `required_moments_are_finite` | [226](../../../src/wheeler.cpp#L226) |
| `attempt_quadrature_count` | [238](../../../src/wheeler.cpp#L238) |
| `compute_quadrature_at_scale` | [350](../../../src/wheeler.cpp#L350) |
| `compute_quadrature` | [423](../../../src/wheeler.cpp#L423) |
| `compute_conditional_quadrature` | [431](../../../src/wheeler.cpp#L431) |
| `select_bounded_quadrature` | [437](../../../src/wheeler.cpp#L437) |

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
| [src/internal/coordinate_scales.hpp](../../../src/internal/coordinate_scales.hpp) | [qmom reconstruction](qmom_reconstruction.md) |
| [src/internal/cqmom_case.hpp](../../../src/internal/cqmom_case.hpp) | [state and representations](state_and_representations.md) |
| [src/internal/embedded_runge_kutta.hpp](../../../src/internal/embedded_runge_kutta.hpp) | [time integration](time_integration.md) |
| [src/internal/householder_qr.hpp](../../../src/internal/householder_qr.hpp) | [conditional systems](conditional_systems.md) |
| [src/internal/jacobi_svd.hpp](../../../src/internal/jacobi_svd.hpp) | [conditional systems](conditional_systems.md) |
| [src/internal/moment_reconstruction.hpp](../../../src/internal/moment_reconstruction.hpp) | [moment evaluation](moment_evaluation.md) |
| [src/internal/moment_source.hpp](../../../src/internal/moment_source.hpp) | [moment sources](moment_sources.md) |
| [src/internal/vandermonde.hpp](../../../src/internal/vandermonde.hpp) | [conditional systems](conditional_systems.md) |
