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
| `CqmomResult::ok` | [388](../../../include/qmomkit/cqmom.hpp#L388) |
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
| `QmomMomentSourceUpdateResult::ok` | [287](../../../include/qmomkit/moment_source.hpp#L287) |
| `CqmomMomentSourceUpdateResult::ok` | [300](../../../include/qmomkit/moment_source.hpp#L300) |

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
| `accepts_unnumbered_realizability` | [61](../../../src/configuration.cpp#L61) |
| `accepts_coordinate_realizability` | [66](../../../src/configuration.cpp#L66) |
| `trim` | [84](../../../src/configuration.cpp#L84) |
| `add_error` | [95](../../../src/configuration.cpp#L95) |
| `parse_size` | [101](../../../src/configuration.cpp#L101) |
| `parse_index` | [117](../../../src/configuration.cpp#L117) |
| `parse_coordinate_realizability_section` | [122](../../../src/configuration.cpp#L122) |
| `is_coordinate_realizability_key` | [131](../../../src/configuration.cpp#L131) |
| `is_configured_coordinate_realizability_key` | [146](../../../src/configuration.cpp#L146) |
| `is_decimal_number` | [159](../../../src/configuration.cpp#L159) |
| `parse_double` | [210](../../../src/configuration.cpp#L210) |
| `parse_size_list` | [226](../../../src/configuration.cpp#L226) |
| `parse_double_list` | [242](../../../src/configuration.cpp#L242) |
| `parse_optional_double_list` | [258](../../../src/configuration.cpp#L258) |
| `validate_names` | [280](../../../src/configuration.cpp#L280) |
| `parse_moment_key` | [301](../../../src/configuration.cpp#L301) |
| `mixture_component_prefix` | [322](../../../src/configuration.cpp#L322) |
| `is_mixture_component_section` | [327](../../../src/configuration.cpp#L327) |
| `is_independent_product_dimension_section` | [338](../../../src/configuration.cpp#L338) |
| `parse_file` | [348](../../../src/configuration.cpp#L348) |
| `require_entry` | [431](../../../src/configuration.cpp#L431) |
| `find_entry` | [442](../../../src/configuration.cpp#L442) |
| `reject_unknown_keys` | [452](../../../src/configuration.cpp#L452) |
| `reject_unknown_sections` | [471](../../../src/configuration.cpp#L471) |
| `validate_format_version` | [491](../../../src/configuration.cpp#L491) |
| `parse_phase_space_node_counts` | [500](../../../src/configuration.cpp#L500) |
| `validate_phase_space_names` | [512](../../../src/configuration.cpp#L512) |
| `validate_qmom_phase_space` | [522](../../../src/configuration.cpp#L522) |
| `validate_cqmom_phase_space` | [535](../../../src/configuration.cpp#L535) |
| `validate_optional_qmom_phase_space` | [555](../../../src/configuration.cpp#L555) |
| `validate_optional_cqmom_phase_space` | [571](../../../src/configuration.cpp#L571) |
| `parse_required_double` | [587](../../../src/configuration.cpp#L587) |
| `parse_optional_double` | [602](../../../src/configuration.cpp#L602) |
| `parse_optional_bool` | [617](../../../src/configuration.cpp#L617) |
| `parse_optional_size` | [638](../../../src/configuration.cpp#L638) |
| `parse_optional_positive_double` | [653](../../../src/configuration.cpp#L653) |
| `infer_qmom_initial_maximum_order` | [670](../../../src/configuration.cpp#L670) |
| `parse_optional_positive_size` | [697](../../../src/configuration.cpp#L697) |
| `parse_distribution` | [714](../../../src/configuration.cpp#L714), [801](../../../src/configuration.cpp#L801) |
| `component_key_allowed` | [806](../../../src/configuration.cpp#L806) |
| `marginal_key_allowed` | [818](../../../src/configuration.cpp#L818) |
| `has_explicit_qmom_moments` | [834](../../../src/configuration.cpp#L834) |
| `has_explicit_cqmom_moments` | [839](../../../src/configuration.cpp#L839) |
| `qmom_allowed_keys` | [844](../../../src/configuration.cpp#L844) |
| `append_time_integration_allowed_keys` | [868](../../../src/configuration.cpp#L868) |
| `append_output_allowed_keys` | [890](../../../src/configuration.cpp#L890) |
| `with_runtime_allowed_keys` | [899](../../../src/configuration.cpp#L899) |
| `time_integration_sections` | [906](../../../src/configuration.cpp#L906) |
| `reject_unknown_qmom_sections` | [916](../../../src/configuration.cpp#L916) |
| `reject_unknown_qmom_keys` | [929](../../../src/configuration.cpp#L929) |
| `reject_unconsumed_indexed_entries` | [963](../../../src/configuration.cpp#L963) |
| `parse_mixture` | [986](../../../src/configuration.cpp#L986) |
| `has_mixture_component_sections` | [1068](../../../src/configuration.cpp#L1068) |
| `cqmom_explicit_allowed_keys` | [1080](../../../src/configuration.cpp#L1080) |
| `cqmom_initial_allowed_keys` | [1088](../../../src/configuration.cpp#L1088) |
| `reject_unknown_cqmom_sections` | [1097](../../../src/configuration.cpp#L1097) |
| `reject_unknown_cqmom_keys` | [1110](../../../src/configuration.cpp#L1110) |
| `reject_duplicate_indexed_section` | [948](../../../src/configuration.cpp#L948) |
| `parse_independent_product` | [1129](../../../src/configuration.cpp#L1129) |
| `checked_product` | [1223](../../../src/configuration.cpp#L1223) |
| `parse_explicit_tensor` | [1234](../../../src/configuration.cpp#L1234) |
| `format_configuration_error` | [1303](../../../src/configuration.cpp#L1303) |
| `format_configuration_errors` | [1319](../../../src/configuration.cpp#L1319) |
| `failure_result` | [1328](../../../src/configuration.cpp#L1328), [1333](../../../src/configuration.cpp#L1333) |
| `file_failure` | [1338](../../../src/configuration.cpp#L1338) |
| `validate_replacement_tensor_structure` | [1343](../../../src/configuration.cpp#L1343) |
| `format_moment_index` | [1366](../../../src/configuration.cpp#L1366) |
| `read_configuration_file` | [1382](../../../src/configuration.cpp#L1382) |
| `has_realizability_section` | [1404](../../../src/configuration.cpp#L1404) |
| `has_unnumbered_realizability_section` | [1414](../../../src/configuration.cpp#L1414) |
| `coordinate_realizability_sections` | [1424](../../../src/configuration.cpp#L1424) |
| `is_realizability_key` | [1438](../../../src/configuration.cpp#L1438) |
| `realizability_policy_fields` | [1443](../../../src/configuration.cpp#L1443) |
| `is_node_advection_key` | [1449](../../../src/configuration.cpp#L1449) |
| `node_advection_allowed_keys` | [1454](../../../src/configuration.cpp#L1454) |
| `parse_qmom_options_from_parsed` | [1463](../../../src/configuration.cpp#L1463) |
| `reject_unknown_qmom_option_keys` | [1491](../../../src/configuration.cpp#L1491) |
| `parse_realizability_check_mode_value` | [1509](../../../src/configuration.cpp#L1509) |
| `parse_realizability_check_mode` | [1528](../../../src/configuration.cpp#L1528) |
| `parse_moment_support` | [1539](../../../src/configuration.cpp#L1539) |
| `parse_phase_space_evolution_scheme_value` | [1562](../../../src/configuration.cpp#L1562) |
| `parse_node_advection_scheme_value` | [1581](../../../src/configuration.cpp#L1581) |
| `parse_node_advection_options_from_parsed` | [1600](../../../src/configuration.cpp#L1600) |
| `parse_realizability_policy_from_parsed` | [1648](../../../src/configuration.cpp#L1648) |
| `parse_coordinate_realizability_policies_from_parsed` | [1729](../../../src/configuration.cpp#L1729) |
| `parse_cqmom_diagnostic_level_value` | [1789](../../../src/configuration.cpp#L1789) |
| `parse_cqmom_diagnostic_level` | [1807](../../../src/configuration.cpp#L1807) |
| `parse_vandermonde_solver_value` | [1817](../../../src/configuration.cpp#L1817) |
| `parse_vandermonde_solver` | [1866](../../../src/configuration.cpp#L1866) |
| `has_output_section` | [1876](../../../src/configuration.cpp#L1876) |
| `is_output_key` | [1886](../../../src/configuration.cpp#L1886) |
| `validate_optional_phase_space_evolution_scheme` | [1891](../../../src/configuration.cpp#L1891) |
| `validate_optional_node_advection_options` | [1900](../../../src/configuration.cpp#L1900) |
| `validate_optional_realizability_policy` | [1905](../../../src/configuration.cpp#L1905) |
| `validate_cqmom_coordinate_realizability_ranks` | [1921](../../../src/configuration.cpp#L1921) |
| `reject_unknown_phase_space_evolution_sections` | [1935](../../../src/configuration.cpp#L1935) |
| `reject_unknown_phase_space_evolution_keys` | [1951](../../../src/configuration.cpp#L1951) |
| `parse_phase_space_evolution_scheme_from_parsed` | [1972](../../../src/configuration.cpp#L1972) |
| `validate_optional_cqmom_diagnostic_level` | [1992](../../../src/configuration.cpp#L1992) |
| `validate_optional_vandermonde_solver` | [1997](../../../src/configuration.cpp#L1997) |
| `parse_coordinate_support` | [2003](../../../src/configuration.cpp#L2003) |
| `parse_cqmom_case_settings` | [2066](../../../src/configuration.cpp#L2066) |
| `parse_output_options_from_parsed` | [2156](../../../src/configuration.cpp#L2156) |
| `parse_qmom_initial_condition_from_parsed` | [2254](../../../src/configuration.cpp#L2254) |
| `parse_cqmom_initial_condition_from_parsed` | [2352](../../../src/configuration.cpp#L2352) |
| `initialize_qmom_moments` | [2475](../../../src/configuration.cpp#L2475) |
| `initialize_cqmom_moment_tensor` | [2519](../../../src/configuration.cpp#L2519) |
| `initialize_cqmom_moments` | [2558](../../../src/configuration.cpp#L2558) |
| `parse_qmom_initial_condition` | [2605](../../../src/configuration.cpp#L2605) |
| `parse_cqmom_initial_condition` | [2611](../../../src/configuration.cpp#L2611) |
| `create_qmom_case` | [2617](../../../src/configuration.cpp#L2617) |
| `create_cqmom_case` | [2644](../../../src/configuration.cpp#L2644) |
| `with_initial_moments` | [2656](../../../src/configuration.cpp#L2656), [2685](../../../src/configuration.cpp#L2685), [2727](../../../src/configuration.cpp#L2727), [2734](../../../src/configuration.cpp#L2734) |
| `with_prepared_metadata_of` | [2707](../../../src/configuration.cpp#L2707) |
| `parse_qmom` | [2741](../../../src/configuration.cpp#L2741) |
| `parse_cqmom` | [2791](../../../src/configuration.cpp#L2791) |
| `parse_cqmom_prepared` | [2816](../../../src/configuration.cpp#L2816) |
| `parse_realizability_policy` | [2844](../../../src/configuration.cpp#L2844) |
| `parse_coordinate_realizability_policies` | [2882](../../../src/configuration.cpp#L2882) |
| `parse_phase_space_evolution_scheme` | [2902](../../../src/configuration.cpp#L2902) |
| `parse_node_advection_options` | [2921](../../../src/configuration.cpp#L2921) |
| `parse_qmom_phase_space_evolution_options` | [2940](../../../src/configuration.cpp#L2940) |
| `parse_cqmom_phase_space_evolution_options` | [2969](../../../src/configuration.cpp#L2969) |
| `parse_output_options` | [2995](../../../src/configuration.cpp#L2995) |
| `load_configuration` | [3014](../../../src/configuration.cpp#L3014) |
| `load_qmom` | [3026](../../../src/configuration.cpp#L3026) |
| `load_cqmom` | [3031](../../../src/configuration.cpp#L3031) |
| `load_cqmom_prepared` | [3036](../../../src/configuration.cpp#L3036) |
| `load_output_options` | [3041](../../../src/configuration.cpp#L3041) |
| `load_node_advection_options` | [3046](../../../src/configuration.cpp#L3046) |
| `load_realizability_policy` | [3051](../../../src/configuration.cpp#L3051) |
| `load_qmom_phase_space_evolution_options` | [3056](../../../src/configuration.cpp#L3056) |
| `load_cqmom_phase_space_evolution_options` | [3061](../../../src/configuration.cpp#L3061) |
| `parse_vandermonde_rcond` | [1848](../../../src/configuration.cpp#L1848) |

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
| `status_from_ssp_rk3` | [24](../../../src/moment_source.cpp#L24) |
| `message_from_status` | [38](../../../src/moment_source.cpp#L38) |
| `diagnostics_from_ssp_rk3` | [56](../../../src/moment_source.cpp#L56) |
| `invalid_result` | [67](../../../src/moment_source.cpp#L67) |
| `invalid_moment_data_result` | [74](../../../src/moment_source.cpp#L74) |
| `failed_moment_data_result` | [81](../../../src/moment_source.cpp#L81) |
| `multi_index_precedes` | [88](../../../src/moment_source.cpp#L88) |
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
| `accepts_coordinate_realizability_policies` | [196](../../../src/moment_source.cpp#L196) |
| `boundary_callback_error` | [217](../../../src/moment_source.cpp#L217) |
| `cqmom_boundary_callback_error` | [222](../../../src/moment_source.cpp#L222) |
| `validate_request` | [227](../../../src/moment_source.cpp#L227), [257](../../../src/moment_source.cpp#L257) |
| `validate_cqmom_source_update_controls` | [267](../../../src/moment_source.cpp#L267) |
| `validate_moment_data_request` | [288](../../../src/moment_source.cpp#L288), [305](../../../src/moment_source.cpp#L305) |
| `validate_trusted_moment_data_request` | [296](../../../src/moment_source.cpp#L296), [315](../../../src/moment_source.cpp#L315) |
| `compose_source_and_boundary` | [325](../../../src/moment_source.cpp#L325) |
| `advance_qmom_moments_with_source` | [371](../../../src/moment_source.cpp#L371), [838](../../../src/moment_source.cpp#L838), [858](../../../src/moment_source.cpp#L858) |
| `compose_cqmom_source_and_boundary` | [424](../../../src/moment_source.cpp#L424) |
| `advance_validated_cqmom_moments_with_source` | [451](../../../src/moment_source.cpp#L451) |
| `QmomSourceContextAccess::create` | [515](../../../src/moment_source.cpp#L515) |
| `CqmomSourceContextAccess::create` | [521](../../../src/moment_source.cpp#L521) |
| `QmomSourceContext::QmomSourceContext` | [531](../../../src/moment_source.cpp#L531), [539](../../../src/moment_source.cpp#L539) |
| `QmomSourceContext::time` | [552](../../../src/moment_source.cpp#L552) |
| `QmomSourceContext::moment_count` | [557](../../../src/moment_source.cpp#L557) |
| `QmomSourceContext::size` | [562](../../../src/moment_source.cpp#L562) |
| `QmomSourceContext::values` | [567](../../../src/moment_source.cpp#L567) |
| `QmomSourceContext::value` | [572](../../../src/moment_source.cpp#L572) |
| `QmomSourceContext::moment` | [577](../../../src/moment_source.cpp#L577) |
| `QmomSourceContext::generalized_moment` | [594](../../../src/moment_source.cpp#L594) |
| `QmomSourceContext::for_each_node` | [600](../../../src/moment_source.cpp#L600) |
| `QmomSourceContext::reconstructed_stage_quadrature` | [614](../../../src/moment_source.cpp#L614) |
| `CqmomSourceContext::CqmomSourceContext` | [632](../../../src/moment_source.cpp#L632), [643](../../../src/moment_source.cpp#L643) |
| `CqmomSourceContext::time` | [657](../../../src/moment_source.cpp#L657) |
| `CqmomSourceContext::moment_count` | [662](../../../src/moment_source.cpp#L662) |
| `CqmomSourceContext::dimension_count` | [667](../../../src/moment_source.cpp#L667) |
| `CqmomSourceContext::size` | [672](../../../src/moment_source.cpp#L672) |
| `CqmomSourceContext::indices` | [677](../../../src/moment_source.cpp#L677) |
| `CqmomSourceContext::values` | [682](../../../src/moment_source.cpp#L682) |
| `CqmomSourceContext::index` | [687](../../../src/moment_source.cpp#L687) |
| `CqmomSourceContext::value` | [692](../../../src/moment_source.cpp#L692), [747](../../../src/moment_source.cpp#L747) |
| `CqmomSourceContext::sorted_index_positions` | [697](../../../src/moment_source.cpp#L697) |
| `CqmomSourceContext::position` | [734](../../../src/moment_source.cpp#L734) |
| `CqmomSourceContext::moment` | [757](../../../src/moment_source.cpp#L757) |
| `CqmomSourceContext::generalized_moment` | [773](../../../src/moment_source.cpp#L773) |
| `CqmomSourceContext::for_each_terminal_node` | [786](../../../src/moment_source.cpp#L786) |
| `CqmomSourceContext::reconstructed_stage_quadrature` | [807](../../../src/moment_source.cpp#L807) |
| `CqmomSourceContext::hierarchy` | [833](../../../src/moment_source.cpp#L833) |
| `advance_cqmom_moments_with_source` | [874](../../../src/moment_source.cpp#L874), [889](../../../src/moment_source.cpp#L889) |
| `advance_trusted_qmom_moments_with_source` | [908](../../../src/moment_source.cpp#L908), [933](../../../src/moment_source.cpp#L933) |
| `advance_trusted_cqmom_moments_with_source` | [960](../../../src/moment_source.cpp#L960), [975](../../../src/moment_source.cpp#L975), [994](../../../src/moment_source.cpp#L994), [1012](../../../src/moment_source.cpp#L1012) |
| `assign_source_rates` | [366](../../../src/moment_source.cpp#L366) |

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
| `apply_batch` | [47](../../../src/phase_space_evolution.cpp#L47) |
| `validate_time_input` | [60](../../../src/phase_space_evolution.cpp#L60) |
| `validate_phase_space_evolution_scheme` | [72](../../../src/phase_space_evolution.cpp#L72) |
| `execute_phase_space_evolution_sequence` | [85](../../../src/phase_space_evolution.cpp#L85) |
| `validate_qmom_input` | [120](../../../src/phase_space_evolution.cpp#L120) |
| `validate_boundary` | [160](../../../src/phase_space_evolution.cpp#L160) |
| `validate_cqmom_boundary` | [168](../../../src/phase_space_evolution.cpp#L168) |
| `successful_reconstruction` | [176](../../../src/phase_space_evolution.cpp#L176) |
| `record_advection_diagnostics` | [181](../../../src/phase_space_evolution.cpp#L181) |
| `record_source_update_diagnostics` | [197](../../../src/phase_space_evolution.cpp#L197) |
| `qmom_from_quadrature` | [204](../../../src/phase_space_evolution.cpp#L204) |
| `add_advection_totals` | [222](../../../src/phase_space_evolution.cpp#L222) |
| `finalize_diagnostics` | [231](../../../src/phase_space_evolution.cpp#L231) |
| `record_failure_diagnostics` | [237](../../../src/phase_space_evolution.cpp#L237) |
| `qmom_failure` | [244](../../../src/phase_space_evolution.cpp#L244) |
| `select_post_source_qmom` | [252](../../../src/phase_space_evolution.cpp#L252) |
| `validate_qmom_advection_input` | [295](../../../src/phase_space_evolution.cpp#L295) |
| `validate_qmom_moment_source_input` | [308](../../../src/phase_space_evolution.cpp#L308) |
| `advance_qmom_advection_stage` | [344](../../../src/phase_space_evolution.cpp#L344) |
| `advance_qmom_moment_source_stage` | [394](../../../src/phase_space_evolution.cpp#L394) |
| `operator<` | [464](../../../src/phase_space_evolution.cpp#L464) |
| `same_active_branch_count` | [477](../../../src/phase_space_evolution.cpp#L477) |
| `active_branch_key` | [484](../../../src/phase_space_evolution.cpp#L484) |
| `path_has_prefix` | [489](../../../src/phase_space_evolution.cpp#L489) |
| `has_changed_strict_ancestor` | [505](../../../src/phase_space_evolution.cpp#L505) |
| `record_cqmom_node_count_transitions` | [523](../../../src/phase_space_evolution.cpp#L523) |
| `advance_qmom_phase_space_impl` | [621](../../../src/phase_space_evolution.cpp#L621) |
| `advance_qmom_advection` | [666](../../../src/phase_space_evolution.cpp#L666) |
| `advance_qmom_moment_source_impl` | [691](../../../src/phase_space_evolution.cpp#L691) |
| `advance_qmom_moment_source` | [718](../../../src/phase_space_evolution.cpp#L718), [724](../../../src/phase_space_evolution.cpp#L724) |
| `advance_qmom_phase_space` | [731](../../../src/phase_space_evolution.cpp#L731), [738](../../../src/phase_space_evolution.cpp#L738) |
| `initialize_cqmom_diagnostics` | [748](../../../src/phase_space_evolution.cpp#L748), [754](../../../src/phase_space_evolution.cpp#L754) |
| `validate_cqmom_moments_for_phase_space_evolution` | [759](../../../src/phase_space_evolution.cpp#L759) |
| `validate_cqmom_state_for_phase_space_evolution` | [766](../../../src/phase_space_evolution.cpp#L766) |
| `validate_cqmom_advection_input` | [771](../../../src/phase_space_evolution.cpp#L771) |
| `validate_cqmom_moment_source_input` | [781](../../../src/phase_space_evolution.cpp#L781) |
| `validate_cqmom_phase_space_evolution_input` | [793](../../../src/phase_space_evolution.cpp#L793) |
| `CqmomStageOutcome::ok` | [820](../../../src/phase_space_evolution.cpp#L820) |
| `take_cqmom_advection_active_counts` | [834](../../../src/phase_space_evolution.cpp#L834) |
| `cqmom_state_failure` | [843](../../../src/phase_space_evolution.cpp#L843) |
| `cqmom_success` | [851](../../../src/phase_space_evolution.cpp#L851) |
| `materialize_cqmom_phase_space_evolution_result` | [856](../../../src/phase_space_evolution.cpp#L856) |
| `cqmom_state_from_accepted_quadrature` | [876](../../../src/phase_space_evolution.cpp#L876) |
| `prepare_cqmom_phase_space_evolution_working_state` | [885](../../../src/phase_space_evolution.cpp#L885) |
| `advance_cqmom_advection_stage` | [921](../../../src/phase_space_evolution.cpp#L921) |
| `advance_cqmom_moment_source_stage` | [980](../../../src/phase_space_evolution.cpp#L980) |
| `advance_cqmom_advection_impl` | [1043](../../../src/phase_space_evolution.cpp#L1043) |
| `advance_cqmom_moment_source_impl` | [1050](../../../src/phase_space_evolution.cpp#L1050) |
| `advance_cqmom_phase_space_state_impl` | [1057](../../../src/phase_space_evolution.cpp#L1057) |
| `advance_cqmom_advection` | [1120](../../../src/phase_space_evolution.cpp#L1120) |
| `advance_cqmom_moment_source_public_impl` | [1142](../../../src/phase_space_evolution.cpp#L1142) |
| `advance_cqmom_moment_source` | [1165](../../../src/phase_space_evolution.cpp#L1165), [1171](../../../src/phase_space_evolution.cpp#L1171) |
| `advance_cqmom_phase_space` | [1178](../../../src/phase_space_evolution.cpp#L1178), [1185](../../../src/phase_space_evolution.cpp#L1185) |
| `advance_qmom_advection_batch` | [1193](../../../src/phase_space_evolution.cpp#L1193) |
| `advance_cqmom_advection_batch` | [1201](../../../src/phase_space_evolution.cpp#L1201) |
| `advance_qmom_moment_source_batch` | [1209](../../../src/phase_space_evolution.cpp#L1209), [1225](../../../src/phase_space_evolution.cpp#L1225) |
| `advance_cqmom_moment_source_batch` | [1217](../../../src/phase_space_evolution.cpp#L1217), [1234](../../../src/phase_space_evolution.cpp#L1234) |
| `advance_qmom_phase_space_batch` | [1243](../../../src/phase_space_evolution.cpp#L1243), [1261](../../../src/phase_space_evolution.cpp#L1261) |
| `advance_cqmom_phase_space_batch` | [1252](../../../src/phase_space_evolution.cpp#L1252), [1271](../../../src/phase_space_evolution.cpp#L1271) |
| `advance_qmom_advection_batch_with_index` | [1281](../../../src/phase_space_evolution.cpp#L1281) |
| `advance_cqmom_advection_batch_with_index` | [1298](../../../src/phase_space_evolution.cpp#L1298) |
| `advance_qmom_moment_source_batch_with_index` | [1315](../../../src/phase_space_evolution.cpp#L1315), [1345](../../../src/phase_space_evolution.cpp#L1345) |
| `advance_cqmom_moment_source_batch_with_index` | [1330](../../../src/phase_space_evolution.cpp#L1330), [1366](../../../src/phase_space_evolution.cpp#L1366) |
| `advance_qmom_phase_space_batch_with_index` | [1387](../../../src/phase_space_evolution.cpp#L1387), [1432](../../../src/phase_space_evolution.cpp#L1432) |
| `advance_cqmom_phase_space_batch_with_index` | [1409](../../../src/phase_space_evolution.cpp#L1409), [1460](../../../src/phase_space_evolution.cpp#L1460) |

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
| `qmom_state_candidate` | [171](../../../src/public_api.cpp#L171) |
| `checked_qmom_state_candidate` | [229](../../../src/public_api.cpp#L229) |
| `QmomCase::state_candidate` | [242](../../../src/public_api.cpp#L242), [247](../../../src/public_api.cpp#L247) |
| `QmomCase::initial_state` | [252](../../../src/public_api.cpp#L252) |
| `format_cqmom_settings` | [276](../../../src/public_api.cpp#L276) |
| `CqmomPreparedCase::CqmomPreparedCase` | [286](../../../src/public_api.cpp#L286), [288](../../../src/public_api.cpp#L288), [293](../../../src/public_api.cpp#L293) |
| `CqmomPreparedCase::operator=` | [297](../../../src/public_api.cpp#L297) |
| `CqmomPreparedCase::hierarchy` | [312](../../../src/public_api.cpp#L312) |
| `CqmomPreparedCase::qmom_options` | [317](../../../src/public_api.cpp#L317) |
| `CqmomPreparedCase::conditional_options` | [322](../../../src/public_api.cpp#L322) |
| `CqmomPreparedCase::diagnostic_level` | [327](../../../src/public_api.cpp#L327) |
| `CqmomPreparedCase::vandermonde_solver` | [332](../../../src/public_api.cpp#L332) |
| `CqmomPreparedCase::abscissa_bounds` | [337](../../../src/public_api.cpp#L337) |
| `cqmom_state_candidate` | [350](../../../src/public_api.cpp#L350) |
| `initial_cqmom_state` | [433](../../../src/public_api.cpp#L433) |
| `prepared_cqmom_state_candidate` | [446](../../../src/public_api.cpp#L446) |
| `CqmomCase::initial_state` | [469](../../../src/public_api.cpp#L469) |
| `CqmomPreparedCase::state_candidate` | [476](../../../src/public_api.cpp#L476), [481](../../../src/public_api.cpp#L481) |
| `CqmomPreparedCase::initial_state` | [487](../../../src/public_api.cpp#L487) |
| `CqmomState::CqmomState` | [492](../../../src/public_api.cpp#L492), [494](../../../src/public_api.cpp#L494) |
| `CqmomState::operator=` | [496](../../../src/public_api.cpp#L496) |
| `CqmomState::hierarchy` | [511](../../../src/public_api.cpp#L511) |
| `CqmomState::quadrature` | [516](../../../src/public_api.cpp#L516) |
| `CqmomState::qmom_options` | [525](../../../src/public_api.cpp#L525) |
| `CqmomState::conditional_options` | [530](../../../src/public_api.cpp#L530) |
| `CqmomState::diagnostic_level` | [535](../../../src/public_api.cpp#L535) |
| `CqmomState::vandermonde_solver` | [540](../../../src/public_api.cpp#L540) |
| `CqmomState::abscissa_bounds` | [545](../../../src/public_api.cpp#L545) |
| `CqmomState::active_branch_counts` | [550](../../../src/public_api.cpp#L550) |
| `CqmomState::moment_indices` | [555](../../../src/public_api.cpp#L555) |
| `CqmomState::moment_values` | [560](../../../src/public_api.cpp#L560) |
| `CqmomState::moment` | [565](../../../src/public_api.cpp#L565) |
| `validate_prepared_metadata` | [578](../../../src/public_api.cpp#L578) |
| `make_prepared_metadata` | [587](../../../src/public_api.cpp#L587) |
| `CqmomConstructionAccess::create_prepared_case` | [598](../../../src/public_api.cpp#L598) |
| `CqmomConstructionAccess::create_state` | [604](../../../src/public_api.cpp#L604), [609](../../../src/public_api.cpp#L609), [615](../../../src/public_api.cpp#L615), [623](../../../src/public_api.cpp#L623) |
| `CqmomConstructionAccess::case_definition` | [630](../../../src/public_api.cpp#L630) |
| `CqmomConstructionAccess::state_definition` | [635](../../../src/public_api.cpp#L635) |
| `CqmomConstructionAccess::state_active_branch_counts` | [640](../../../src/public_api.cpp#L640) |
| `CqmomConstructionAccess::state_accepted_quadrature` | [645](../../../src/public_api.cpp#L645) |
| `CqmomConstructionAccess::case_prepared_metadata` | [650](../../../src/public_api.cpp#L650) |
| `CqmomConstructionAccess::state_prepared_metadata` | [655](../../../src/public_api.cpp#L655) |
| `CqmomConstructionAccess::create_quadrature` | [660](../../../src/public_api.cpp#L660) |
| `CqmomQuadrature::CqmomQuadrature` | [670](../../../src/public_api.cpp#L670) |
| `CqmomQuadrature::weight` | [677](../../../src/public_api.cpp#L677) |
| `CqmomQuadrature::abscissa` | [686](../../../src/public_api.cpp#L686) |
| `CqmomQuadrature::node` | [699](../../../src/public_api.cpp#L699) |
| `to_cqmom_status` | [712](../../../src/public_api.cpp#L712) |
| `to_cqmom_diagnostics` | [728](../../../src/public_api.cpp#L728), [757](../../../src/public_api.cpp#L757) |
| `terminal_paths_for_public_quadrature` | [768](../../../src/public_api.cpp#L768) |
| `failed_phase_space_dimension` | [790](../../../src/public_api.cpp#L790) |
| `to_cqmom_result` | [800](../../../src/public_api.cpp#L800) |
| `compute_cqmom` | [820](../../../src/public_api.cpp#L820), [883](../../../src/public_api.cpp#L883) |
| `compute_cqmom_internal` | [844](../../../src/public_api.cpp#L844) |
| `invalid_cqmom_hierarchy_result` | [810](../../../src/public_api.cpp#L810) |

### `src/realizability.cpp`

[Source](../../../src/realizability.cpp) · [Responsibilities, contracts, and
evidence](realizability.md)

| Function or overload family | Definition lines |
|---|---|
| `RealizabilityResult::realizable` | [19](../../../src/realizability.cpp#L19) |
| `RealizabilityResult::on_boundary` | [24](../../../src/realizability.cpp#L24) |
| `support_name` | [81](../../../src/realizability.cpp#L81) |
| `make_result` | [95](../../../src/realizability.cpp#L95) |
| `vector_scale` | [101](../../../src/realizability.cpp#L101) |
| `classify_eigenvalue` | [129](../../../src/realizability.cpp#L129) |
| `hankel_matrix` | [146](../../../src/realizability.cpp#L146) |
| `hankel_uncertainty` | [164](../../../src/realizability.cpp#L164) |
| `analyze_matrix` | [177](../../../src/realizability.cpp#L177) |
| `update_limiting_matrix` | [263](../../../src/realizability.cpp#L263) |
| `incorporate_matrix` | [277](../../../src/realizability.cpp#L277) |
| `lower_localizing_sequence` | [306](../../../src/realizability.cpp#L306) |
| `upper_localizing_sequence` | [319](../../../src/realizability.cpp#L319) |
| `interval_localizing_sequence` | [334](../../../src/realizability.cpp#L334) |
| `extension_column` | [348](../../../src/realizability.cpp#L348) |
| `euclidean_norm` | [359](../../../src/realizability.cpp#L359) |
| `projection_on_eigenvector` | [369](../../../src/realizability.cpp#L369) |
| `kernel_projection_tolerance` | [385](../../../src/realizability.cpp#L385) |
| `check_range_consistency` | [403](../../../src/realizability.cpp#L403) |
| `incorporate_range_consistency` | [429](../../../src/realizability.cpp#L429) |
| `numerical_rank` | [454](../../../src/realizability.cpp#L454) |
| `singular_as_given` | [464](../../../src/realizability.cpp#L464) |
| `finish_classification` | [470](../../../src/realizability.cpp#L470) |
| `validate_common_inputs` | [486](../../../src/realizability.cpp#L486) |
| `check_hamburger` | [521](../../../src/realizability.cpp#L521) |
| `check_stieltjes` | [575](../../../src/realizability.cpp#L575) |
| `check_hausdorff` | [631](../../../src/realizability.cpp#L631) |
| `TestedMoments::bound` | [701](../../../src/realizability.cpp#L701) |
| `standardized_moments` | [707](../../../src/realizability.cpp#L707) |
| `with_support_context` | [781](../../../src/realizability.cpp#L781) |
| `check_realizability` | [792](../../../src/realizability.cpp#L792), [821](../../../src/realizability.cpp#L821) |

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
