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
| `CqmomResult::ok` | [331](../../../include/qmomkit/cqmom.hpp#L331) |
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
| `QmomMomentSourceUpdateResult::ok` | [277](../../../include/qmomkit/moment_source.hpp#L277) |
| `CqmomMomentSourceUpdateResult::ok` | [290](../../../include/qmomkit/moment_source.hpp#L290) |

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
| `QmomResult::ok` | [40](../../../include/qmomkit/qmom.hpp#L40) |
| `QmomState::moment_values` | [70](../../../include/qmomkit/qmom.hpp#L70) |
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
| `SspRk3Result::ok` | [87](../../../include/qmomkit/runge_kutta.hpp#L87) |

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
| `append_time_integration_allowed_keys` | [869](../../../src/configuration.cpp#L869) |
| `append_output_allowed_keys` | [892](../../../src/configuration.cpp#L892) |
| `with_runtime_allowed_keys` | [901](../../../src/configuration.cpp#L901) |
| `time_integration_sections` | [908](../../../src/configuration.cpp#L908) |
| `reject_unknown_qmom_sections` | [918](../../../src/configuration.cpp#L918) |
| `reject_unknown_qmom_keys` | [931](../../../src/configuration.cpp#L931) |
| `reject_unconsumed_indexed_entries` | [966](../../../src/configuration.cpp#L966) |
| `parse_mixture` | [989](../../../src/configuration.cpp#L989) |
| `has_mixture_component_sections` | [1059](../../../src/configuration.cpp#L1059) |
| `cqmom_explicit_allowed_keys` | [1071](../../../src/configuration.cpp#L1071) |
| `cqmom_initial_allowed_keys` | [1080](../../../src/configuration.cpp#L1080) |
| `reject_unknown_cqmom_sections` | [1089](../../../src/configuration.cpp#L1089) |
| `reject_unknown_cqmom_keys` | [1102](../../../src/configuration.cpp#L1102) |
| `reject_duplicate_indexed_section` | [951](../../../src/configuration.cpp#L951) |
| `parse_independent_product` | [1113](../../../src/configuration.cpp#L1113) |
| `checked_product` | [1207](../../../src/configuration.cpp#L1207) |
| `parse_explicit_tensor` | [1218](../../../src/configuration.cpp#L1218) |
| `format_configuration_error` | [1281](../../../src/configuration.cpp#L1281) |
| `format_configuration_errors` | [1297](../../../src/configuration.cpp#L1297) |
| `failure_result` | [1306](../../../src/configuration.cpp#L1306), [1311](../../../src/configuration.cpp#L1311) |
| `file_failure` | [1316](../../../src/configuration.cpp#L1316) |
| `validate_replacement_tensor_structure` | [1321](../../../src/configuration.cpp#L1321) |
| `format_moment_index` | [1344](../../../src/configuration.cpp#L1344) |
| `read_configuration_file` | [1360](../../../src/configuration.cpp#L1360) |
| `has_realizability_section` | [1383](../../../src/configuration.cpp#L1383) |
| `has_unnumbered_realizability_section` | [1393](../../../src/configuration.cpp#L1393) |
| `coordinate_realizability_sections` | [1403](../../../src/configuration.cpp#L1403) |
| `is_realizability_key` | [1417](../../../src/configuration.cpp#L1417) |
| `realizability_policy_fields` | [1422](../../../src/configuration.cpp#L1422) |
| `is_node_advection_key` | [1428](../../../src/configuration.cpp#L1428) |
| `node_advection_allowed_keys` | [1433](../../../src/configuration.cpp#L1433) |
| `parse_qmom_options_from_parsed` | [1442](../../../src/configuration.cpp#L1442) |
| `parse_automatic_scaling_from_parsed` | [1470](../../../src/configuration.cpp#L1470) |
| `reject_unknown_qmom_option_keys` | [1489](../../../src/configuration.cpp#L1489) |
| `parse_realizability_check_mode_value` | [1507](../../../src/configuration.cpp#L1507) |
| `parse_realizability_check_mode` | [1526](../../../src/configuration.cpp#L1526) |
| `parse_moment_support` | [1537](../../../src/configuration.cpp#L1537) |
| `parse_phase_space_evolution_scheme_value` | [1560](../../../src/configuration.cpp#L1560) |
| `parse_node_advection_scheme_value` | [1579](../../../src/configuration.cpp#L1579) |
| `parse_node_advection_options_from_parsed` | [1598](../../../src/configuration.cpp#L1598) |
| `parse_realizability_policy_from_parsed` | [1646](../../../src/configuration.cpp#L1646) |
| `parse_coordinate_realizability_policies_from_parsed` | [1727](../../../src/configuration.cpp#L1727) |
| `parse_cqmom_diagnostic_level_value` | [1767](../../../src/configuration.cpp#L1767) |
| `parse_cqmom_diagnostic_level` | [1785](../../../src/configuration.cpp#L1785) |
| `parse_vandermonde_solver_value` | [1795](../../../src/configuration.cpp#L1795) |
| `parse_vandermonde_solver` | [1844](../../../src/configuration.cpp#L1844) |
| `has_output_section` | [1854](../../../src/configuration.cpp#L1854) |
| `is_output_key` | [1864](../../../src/configuration.cpp#L1864) |
| `validate_optional_phase_space_evolution_scheme` | [1869](../../../src/configuration.cpp#L1869) |
| `validate_optional_node_advection_options` | [1878](../../../src/configuration.cpp#L1878) |
| `validate_optional_realizability_policy` | [1883](../../../src/configuration.cpp#L1883) |
| `validate_cqmom_coordinate_realizability_ranks` | [1899](../../../src/configuration.cpp#L1899) |
| `reject_unknown_phase_space_evolution_sections` | [1913](../../../src/configuration.cpp#L1913) |
| `reject_unknown_phase_space_evolution_keys` | [1929](../../../src/configuration.cpp#L1929) |
| `parse_phase_space_evolution_scheme_from_parsed` | [1950](../../../src/configuration.cpp#L1950) |
| `validate_optional_cqmom_diagnostic_level` | [1970](../../../src/configuration.cpp#L1970) |
| `validate_optional_vandermonde_solver` | [1975](../../../src/configuration.cpp#L1975) |
| `parse_coordinate_support` | [1981](../../../src/configuration.cpp#L1981) |
| `parse_cqmom_case_settings` | [2044](../../../src/configuration.cpp#L2044) |
| `parse_output_options_from_parsed` | [2161](../../../src/configuration.cpp#L2161) |
| `parse_qmom_initial_condition_from_parsed` | [2259](../../../src/configuration.cpp#L2259) |
| `parse_cqmom_initial_condition_from_parsed` | [2357](../../../src/configuration.cpp#L2357) |
| `initialize_qmom_moments` | [2480](../../../src/configuration.cpp#L2480) |
| `initialize_cqmom_moment_tensor` | [2524](../../../src/configuration.cpp#L2524) |
| `initialize_cqmom_moments` | [2563](../../../src/configuration.cpp#L2563) |
| `parse_qmom_initial_condition` | [2610](../../../src/configuration.cpp#L2610) |
| `parse_cqmom_initial_condition` | [2616](../../../src/configuration.cpp#L2616) |
| `create_qmom_case` | [2622](../../../src/configuration.cpp#L2622) |
| `create_cqmom_case` | [2649](../../../src/configuration.cpp#L2649) |
| `with_initial_moments` | [2667](../../../src/configuration.cpp#L2667), [2696](../../../src/configuration.cpp#L2696), [2738](../../../src/configuration.cpp#L2738), [2745](../../../src/configuration.cpp#L2745) |
| `with_prepared_metadata_of` | [2718](../../../src/configuration.cpp#L2718) |
| `parse_qmom` | [2752](../../../src/configuration.cpp#L2752) |
| `parse_cqmom` | [2807](../../../src/configuration.cpp#L2807) |
| `parse_cqmom_prepared` | [2832](../../../src/configuration.cpp#L2832) |
| `parse_realizability_policy` | [2861](../../../src/configuration.cpp#L2861) |
| `parse_coordinate_realizability_policies` | [2899](../../../src/configuration.cpp#L2899) |
| `parse_phase_space_evolution_scheme` | [2919](../../../src/configuration.cpp#L2919) |
| `parse_node_advection_options` | [2938](../../../src/configuration.cpp#L2938) |
| `parse_qmom_phase_space_evolution_options` | [2957](../../../src/configuration.cpp#L2957) |
| `parse_cqmom_phase_space_evolution_options` | [2987](../../../src/configuration.cpp#L2987) |
| `parse_output_options` | [3013](../../../src/configuration.cpp#L3013) |
| `load_configuration` | [3032](../../../src/configuration.cpp#L3032) |
| `load_qmom` | [3044](../../../src/configuration.cpp#L3044) |
| `load_cqmom` | [3049](../../../src/configuration.cpp#L3049) |
| `load_cqmom_prepared` | [3054](../../../src/configuration.cpp#L3054) |
| `load_output_options` | [3059](../../../src/configuration.cpp#L3059) |
| `load_node_advection_options` | [3064](../../../src/configuration.cpp#L3064) |
| `load_realizability_policy` | [3069](../../../src/configuration.cpp#L3069) |
| `load_qmom_phase_space_evolution_options` | [3074](../../../src/configuration.cpp#L3074) |
| `load_cqmom_phase_space_evolution_options` | [3079](../../../src/configuration.cpp#L3079) |
| `parse_vandermonde_rcond` | [1826](../../../src/configuration.cpp#L1826) |

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
| `with_moments` | [30](../../../src/internal/cqmom_case.hpp#L30) |

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
| `make_identity` | [124](../../../src/internal/linear_algebra.cpp#L124) |
| `residual_norm` | [134](../../../src/internal/linear_algebra.cpp#L134) |
| `condition_estimate` | [154](../../../src/internal/linear_algebra.cpp#L154) |
| `power_of_two_shift` | [170](../../../src/internal/linear_algebra.cpp#L170) |
| `rotate` | [194](../../../src/internal/linear_algebra.cpp#L194) |
| `symmetric_jacobi_eigen` | [239](../../../src/internal/linear_algebra.cpp#L239) |
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
| `RecursiveQuadratureResult::ok` | [227](../../../src/internal/recursive_quadrature.hpp#L227) |

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
| `advance_qmom_moments_with_source` | [405](../../../src/moment_source.cpp#L405), [873](../../../src/moment_source.cpp#L873), [893](../../../src/moment_source.cpp#L893) |
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
| `QmomSourceContext::for_each_node` | [634](../../../src/moment_source.cpp#L634) |
| `CqmomSourceContext::CqmomSourceContext` | [666](../../../src/moment_source.cpp#L666), [677](../../../src/moment_source.cpp#L677) |
| `CqmomSourceContext::time` | [691](../../../src/moment_source.cpp#L691) |
| `CqmomSourceContext::moment_count` | [696](../../../src/moment_source.cpp#L696) |
| `CqmomSourceContext::dimension_count` | [701](../../../src/moment_source.cpp#L701) |
| `CqmomSourceContext::size` | [706](../../../src/moment_source.cpp#L706) |
| `CqmomSourceContext::indices` | [711](../../../src/moment_source.cpp#L711) |
| `CqmomSourceContext::values` | [716](../../../src/moment_source.cpp#L716) |
| `CqmomSourceContext::index` | [721](../../../src/moment_source.cpp#L721) |
| `CqmomSourceContext::value` | [726](../../../src/moment_source.cpp#L726), [781](../../../src/moment_source.cpp#L781) |
| `CqmomSourceContext::position` | [768](../../../src/moment_source.cpp#L768) |
| `CqmomSourceContext::moment` | [791](../../../src/moment_source.cpp#L791) |
| `CqmomSourceContext::generalized_moment` | [807](../../../src/moment_source.cpp#L807) |
| `CqmomSourceContext::for_each_terminal_node` | [820](../../../src/moment_source.cpp#L820) |
| `CqmomSourceContext::reconstructed_stage_quadrature` | [841](../../../src/moment_source.cpp#L841) |
| `CqmomSourceContext::hierarchy` | [868](../../../src/moment_source.cpp#L868) |
| `advance_cqmom_moments_with_source` | [909](../../../src/moment_source.cpp#L909), [924](../../../src/moment_source.cpp#L924) |
| `advance_trusted_qmom_moments_with_source` | [943](../../../src/moment_source.cpp#L943), [968](../../../src/moment_source.cpp#L968) |
| `advance_trusted_cqmom_moments_with_source` | [995](../../../src/moment_source.cpp#L995), [1010](../../../src/moment_source.cpp#L1010), [1029](../../../src/moment_source.cpp#L1029), [1047](../../../src/moment_source.cpp#L1047) |
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
| `integrate_nodes_ssp_rk3` | [208](../../../src/node_advection.cpp#L208) |
| `integrate_nodes` | [279](../../../src/node_advection.cpp#L279) |
| `scheme_name` | [310](../../../src/node_advection.cpp#L310) |
| `PooledCoordinateBuffer::PooledCoordinateBuffer` | [340](../../../src/node_advection.cpp#L340) |
| `PooledCoordinateBuffer::~PooledCoordinateBuffer` | [361](../../../src/node_advection.cpp#L361) |
| `PooledCoordinateBuffer::free_list` | [373](../../../src/node_advection.cpp#L373) |
| `make_phase_space_velocity` | [467](../../../src/node_advection.cpp#L467) |
| `advance_nodes` | [384](../../../src/node_advection.cpp#L384), [412](../../../src/node_advection.cpp#L412), [423](../../../src/node_advection.cpp#L423) |
| `format_node_advection_options` | [507](../../../src/node_advection.cpp#L507) |

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
| `MomentWriter::write_aligned_row` | [428](../../../src/output.cpp#L428) |
| `MomentWriter::write` | [459](../../../src/output.cpp#L459), [478](../../../src/output.cpp#L478), [501](../../../src/output.cpp#L501), [520](../../../src/output.cpp#L520), [539](../../../src/output.cpp#L539), [562](../../../src/output.cpp#L562), [590](../../../src/output.cpp#L590), [609](../../../src/output.cpp#L609), [623](../../../src/output.cpp#L623), [642](../../../src/output.cpp#L642) |
| `QuadratureWriter::for_qmom` | [656](../../../src/output.cpp#L656) |
| `QuadratureWriter::for_cqmom` | [661](../../../src/output.cpp#L661) |
| `QuadratureWriter::QuadratureWriter` | [666](../../../src/output.cpp#L666) |
| `QuadratureWriter::write_header` | [682](../../../src/output.cpp#L682) |
| `QuadratureWriter::write_row` | [687](../../../src/output.cpp#L687) |
| `QuadratureWriter::write_qmom_row` | [714](../../../src/output.cpp#L714) |
| `QuadratureWriter::write` | [752](../../../src/output.cpp#L752), [757](../../../src/output.cpp#L757), [762](../../../src/output.cpp#L762) |
| `centered_column_label` | [24](../../../src/output.cpp#L24) |
| `generalized_moment_label` | [71](../../../src/output.cpp#L71) |
| `append_generalized_labels` | [87](../../../src/output.cpp#L87) |
| `validate_generalized_moment_orders` | [96](../../../src/output.cpp#L96) |

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
| `validate_boundary` | [162](../../../src/phase_space_evolution.cpp#L162) |
| `validate_cqmom_boundary` | [170](../../../src/phase_space_evolution.cpp#L170) |
| `successful_reconstruction` | [178](../../../src/phase_space_evolution.cpp#L178) |
| `record_advection_diagnostics` | [183](../../../src/phase_space_evolution.cpp#L183) |
| `record_source_update_diagnostics` | [199](../../../src/phase_space_evolution.cpp#L199) |
| `qmom_from_quadrature` | [206](../../../src/phase_space_evolution.cpp#L206) |
| `add_advection_totals` | [224](../../../src/phase_space_evolution.cpp#L224) |
| `finalize_diagnostics` | [233](../../../src/phase_space_evolution.cpp#L233) |
| `record_failure_diagnostics` | [239](../../../src/phase_space_evolution.cpp#L239) |
| `qmom_failure` | [246](../../../src/phase_space_evolution.cpp#L246) |
| `select_post_source_qmom` | [254](../../../src/phase_space_evolution.cpp#L254) |
| `validate_qmom_advection_input` | [318](../../../src/phase_space_evolution.cpp#L318) |
| `validate_qmom_moment_source_input` | [331](../../../src/phase_space_evolution.cpp#L331) |
| `advance_qmom_advection_stage` | [367](../../../src/phase_space_evolution.cpp#L367) |
| `advance_qmom_moment_source_stage` | [417](../../../src/phase_space_evolution.cpp#L417) |
| `operator<` | [487](../../../src/phase_space_evolution.cpp#L487) |
| `active_branch_key` | [500](../../../src/phase_space_evolution.cpp#L500) |
| `path_has_prefix` | [505](../../../src/phase_space_evolution.cpp#L505) |
| `has_changed_strict_ancestor` | [521](../../../src/phase_space_evolution.cpp#L521) |
| `record_cqmom_node_count_transitions` | [539](../../../src/phase_space_evolution.cpp#L539) |
| `advance_qmom_phase_space_impl` | [629](../../../src/phase_space_evolution.cpp#L629) |
| `advance_qmom_advection` | [674](../../../src/phase_space_evolution.cpp#L674) |
| `advance_qmom_moment_source_impl` | [729](../../../src/phase_space_evolution.cpp#L729) |
| `advance_qmom_moment_source` | [726](../../../src/phase_space_evolution.cpp#L726), [732](../../../src/phase_space_evolution.cpp#L732) |
| `advance_qmom_phase_space` | [739](../../../src/phase_space_evolution.cpp#L739), [746](../../../src/phase_space_evolution.cpp#L746) |
| `initialize_cqmom_diagnostics` | [756](../../../src/phase_space_evolution.cpp#L756), [762](../../../src/phase_space_evolution.cpp#L762) |
| `validate_cqmom_moments_for_phase_space_evolution` | [767](../../../src/phase_space_evolution.cpp#L767) |
| `validate_cqmom_state_for_phase_space_evolution` | [774](../../../src/phase_space_evolution.cpp#L774) |
| `validate_cqmom_advection_input` | [779](../../../src/phase_space_evolution.cpp#L779) |
| `validate_cqmom_moment_source_input` | [789](../../../src/phase_space_evolution.cpp#L789) |
| `validate_cqmom_phase_space_evolution_input` | [801](../../../src/phase_space_evolution.cpp#L801) |
| `CqmomStageOutcome::ok` | [828](../../../src/phase_space_evolution.cpp#L828) |
| `take_cqmom_advection_active_counts` | [842](../../../src/phase_space_evolution.cpp#L842) |
| `cqmom_state_failure` | [851](../../../src/phase_space_evolution.cpp#L851) |
| `cqmom_success` | [859](../../../src/phase_space_evolution.cpp#L859) |
| `materialize_cqmom_phase_space_evolution_result` | [864](../../../src/phase_space_evolution.cpp#L864) |
| `cqmom_state_from_accepted_quadrature` | [884](../../../src/phase_space_evolution.cpp#L884) |
| `prepare_cqmom_phase_space_evolution_working_state` | [893](../../../src/phase_space_evolution.cpp#L893) |
| `advance_cqmom_advection_stage` | [931](../../../src/phase_space_evolution.cpp#L931) |
| `advance_cqmom_moment_source_stage` | [990](../../../src/phase_space_evolution.cpp#L990) |
| `advance_cqmom_advection_impl` | [1054](../../../src/phase_space_evolution.cpp#L1054) |
| `advance_cqmom_moment_source_impl` | [1061](../../../src/phase_space_evolution.cpp#L1061) |
| `advance_cqmom_phase_space_state_impl` | [1068](../../../src/phase_space_evolution.cpp#L1068) |
| `advance_cqmom_advection` | [1131](../../../src/phase_space_evolution.cpp#L1131) |
| `advance_cqmom_moment_source_public_impl` | [1153](../../../src/phase_space_evolution.cpp#L1153) |
| `advance_cqmom_moment_source` | [1176](../../../src/phase_space_evolution.cpp#L1176), [1182](../../../src/phase_space_evolution.cpp#L1182) |
| `advance_cqmom_phase_space` | [1189](../../../src/phase_space_evolution.cpp#L1189), [1196](../../../src/phase_space_evolution.cpp#L1196) |
| `advance_qmom_advection_batch` | [1204](../../../src/phase_space_evolution.cpp#L1204) |
| `advance_cqmom_advection_batch` | [1212](../../../src/phase_space_evolution.cpp#L1212) |
| `advance_qmom_moment_source_batch` | [1220](../../../src/phase_space_evolution.cpp#L1220), [1236](../../../src/phase_space_evolution.cpp#L1236) |
| `advance_cqmom_moment_source_batch` | [1228](../../../src/phase_space_evolution.cpp#L1228), [1245](../../../src/phase_space_evolution.cpp#L1245) |
| `advance_qmom_phase_space_batch` | [1254](../../../src/phase_space_evolution.cpp#L1254), [1272](../../../src/phase_space_evolution.cpp#L1272) |
| `advance_cqmom_phase_space_batch` | [1263](../../../src/phase_space_evolution.cpp#L1263), [1282](../../../src/phase_space_evolution.cpp#L1282) |
| `advance_qmom_advection_batch_with_index` | [1292](../../../src/phase_space_evolution.cpp#L1292) |
| `advance_cqmom_advection_batch_with_index` | [1309](../../../src/phase_space_evolution.cpp#L1309) |
| `advance_qmom_moment_source_batch_with_index` | [1326](../../../src/phase_space_evolution.cpp#L1326), [1356](../../../src/phase_space_evolution.cpp#L1356) |
| `advance_cqmom_moment_source_batch_with_index` | [1341](../../../src/phase_space_evolution.cpp#L1341), [1377](../../../src/phase_space_evolution.cpp#L1377) |
| `advance_qmom_phase_space_batch_with_index` | [1398](../../../src/phase_space_evolution.cpp#L1398), [1443](../../../src/phase_space_evolution.cpp#L1443) |
| `advance_cqmom_phase_space_batch_with_index` | [1420](../../../src/phase_space_evolution.cpp#L1420), [1471](../../../src/phase_space_evolution.cpp#L1471) |

### `src/public_api.cpp`

[Source](../../../src/public_api.cpp) · [Responsibilities, contracts, and
evidence](state_and_representations.md)

| Function or overload family | Definition lines |
|---|---|
| `to_wheeler_options` | [27](../../../src/public_api.cpp#L27) |
| `compute_qmom_unscaled` | [34](../../../src/public_api.cpp#L34) |
| `compute_qmom_with_automatic_scaling` | [51](../../../src/public_api.cpp#L51) |
| `to_qmom_status` | [73](../../../src/public_api.cpp#L73) |
| `hierarchy` | [90](../../../src/public_api.cpp#L90) |
| `qmom_options` | [95](../../../src/public_api.cpp#L95) |
| `conditional_options` | [100](../../../src/public_api.cpp#L100) |
| `adaptive_branch_counts` | [105](../../../src/public_api.cpp#L105) |
| `diagnostic_level` | [110](../../../src/public_api.cpp#L110) |
| `vandermonde_solver` | [115](../../../src/public_api.cpp#L115) |
| `coordinate_support` | [120](../../../src/public_api.cpp#L120) |
| `diagnostic_level_name` | [125](../../../src/public_api.cpp#L125) |
| `vandermonde_solver_name` | [139](../../../src/public_api.cpp#L139) |
| `compute_qmom` | [160](../../../src/public_api.cpp#L160), [165](../../../src/public_api.cpp#L165), [171](../../../src/public_api.cpp#L171) |
| `QmomState::QmomState` | [176](../../../src/public_api.cpp#L176) |
| `QmomCase::initial_state` | [183](../../../src/public_api.cpp#L183) |
| `format_cqmom_settings` | [207](../../../src/public_api.cpp#L207) |
| `CqmomPreparedCase::CqmomPreparedCase` | [218](../../../src/public_api.cpp#L218), [220](../../../src/public_api.cpp#L220), [225](../../../src/public_api.cpp#L225), [238](../../../src/public_api.cpp#L238) |
| `CqmomPreparedCase::operator=` | [229](../../../src/public_api.cpp#L229), [240](../../../src/public_api.cpp#L240) |
| `CqmomPreparedCase::hierarchy` | [244](../../../src/public_api.cpp#L244) |
| `CqmomPreparedCase::qmom_options` | [249](../../../src/public_api.cpp#L249) |
| `CqmomPreparedCase::conditional_options` | [254](../../../src/public_api.cpp#L254) |
| `CqmomPreparedCase::adaptive_branch_counts` | [259](../../../src/public_api.cpp#L259) |
| `CqmomPreparedCase::diagnostic_level` | [264](../../../src/public_api.cpp#L264) |
| `CqmomPreparedCase::vandermonde_solver` | [269](../../../src/public_api.cpp#L269) |
| `CqmomPreparedCase::abscissa_bounds` | [274](../../../src/public_api.cpp#L274) |
| `initial_cqmom_state` | [284](../../../src/public_api.cpp#L284) |
| `CqmomCase::initial_state` | [301](../../../src/public_api.cpp#L301) |
| `CqmomPreparedCase::initial_state` | [308](../../../src/public_api.cpp#L308) |
| `CqmomState::CqmomState` | [313](../../../src/public_api.cpp#L313), [315](../../../src/public_api.cpp#L315), [326](../../../src/public_api.cpp#L326) |
| `CqmomState::operator=` | [317](../../../src/public_api.cpp#L317), [328](../../../src/public_api.cpp#L328) |
| `CqmomState::hierarchy` | [332](../../../src/public_api.cpp#L332) |
| `CqmomState::quadrature` | [337](../../../src/public_api.cpp#L337) |
| `CqmomState::qmom_options` | [346](../../../src/public_api.cpp#L346) |
| `CqmomState::conditional_options` | [351](../../../src/public_api.cpp#L351) |
| `CqmomState::adaptive_branch_counts` | [356](../../../src/public_api.cpp#L356) |
| `CqmomState::diagnostic_level` | [361](../../../src/public_api.cpp#L361) |
| `CqmomState::vandermonde_solver` | [366](../../../src/public_api.cpp#L366) |
| `CqmomState::abscissa_bounds` | [371](../../../src/public_api.cpp#L371) |
| `CqmomState::active_branch_counts` | [376](../../../src/public_api.cpp#L376) |
| `CqmomState::moment_indices` | [381](../../../src/public_api.cpp#L381) |
| `CqmomState::moment_values` | [386](../../../src/public_api.cpp#L386) |
| `CqmomState::moment` | [391](../../../src/public_api.cpp#L391) |
| `validate_prepared_metadata` | [404](../../../src/public_api.cpp#L404) |
| `make_prepared_metadata` | [413](../../../src/public_api.cpp#L413) |
| `CqmomConstructionAccess::create_prepared_case` | [424](../../../src/public_api.cpp#L424) |
| `CqmomConstructionAccess::create_state` | [430](../../../src/public_api.cpp#L430), [435](../../../src/public_api.cpp#L435), [441](../../../src/public_api.cpp#L441), [449](../../../src/public_api.cpp#L449) |
| `CqmomConstructionAccess::case_definition` | [456](../../../src/public_api.cpp#L456) |
| `CqmomConstructionAccess::state_definition` | [461](../../../src/public_api.cpp#L461) |
| `CqmomConstructionAccess::state_active_branch_counts` | [466](../../../src/public_api.cpp#L466) |
| `CqmomConstructionAccess::state_accepted_quadrature` | [471](../../../src/public_api.cpp#L471) |
| `CqmomConstructionAccess::case_prepared_metadata` | [476](../../../src/public_api.cpp#L476) |
| `CqmomConstructionAccess::state_prepared_metadata` | [481](../../../src/public_api.cpp#L481) |
| `CqmomConstructionAccess::create_quadrature` | [486](../../../src/public_api.cpp#L486) |
| `CqmomQuadrature::CqmomQuadrature` | [496](../../../src/public_api.cpp#L496) |
| `CqmomQuadrature::weight` | [503](../../../src/public_api.cpp#L503) |
| `CqmomQuadrature::abscissa` | [512](../../../src/public_api.cpp#L512) |
| `CqmomQuadrature::node` | [525](../../../src/public_api.cpp#L525) |
| `to_cqmom_status` | [538](../../../src/public_api.cpp#L538) |
| `to_cqmom_diagnostics` | [554](../../../src/public_api.cpp#L554), [596](../../../src/public_api.cpp#L596) |
| `terminal_paths_for_public_quadrature` | [608](../../../src/public_api.cpp#L608) |
| `failed_phase_space_dimension` | [630](../../../src/public_api.cpp#L630) |
| `to_cqmom_result` | [640](../../../src/public_api.cpp#L640) |
| `compute_cqmom` | [666](../../../src/public_api.cpp#L666), [752](../../../src/public_api.cpp#L752) |
| `compute_cqmom_internal` | [714](../../../src/public_api.cpp#L714) |
| `invalid_cqmom_hierarchy_result` | [656](../../../src/public_api.cpp#L656) |

### `src/realizability.cpp`

[Source](../../../src/realizability.cpp) · [Responsibilities, contracts, and
evidence](realizability.md)

| Function or overload family | Definition lines |
|---|---|
| `RealizabilityResult::realizable` | [19](../../../src/realizability.cpp#L19) |
| `RealizabilityResult::on_boundary` | [24](../../../src/realizability.cpp#L24) |
| `support_name` | [82](../../../src/realizability.cpp#L82) |
| `make_result` | [96](../../../src/realizability.cpp#L96) |
| `vector_scale` | [102](../../../src/realizability.cpp#L102) |
| `classify_eigenvalue` | [124](../../../src/realizability.cpp#L124) |
| `hankel_matrix` | [141](../../../src/realizability.cpp#L141) |
| `hankel_uncertainty` | [159](../../../src/realizability.cpp#L159) |
| `analyze_matrix` | [172](../../../src/realizability.cpp#L172) |
| `update_limiting_matrix` | [258](../../../src/realizability.cpp#L258) |
| `incorporate_matrix` | [272](../../../src/realizability.cpp#L272) |
| `lower_localizing_sequence` | [301](../../../src/realizability.cpp#L301) |
| `upper_localizing_sequence` | [314](../../../src/realizability.cpp#L314) |
| `interval_localizing_sequence` | [329](../../../src/realizability.cpp#L329) |
| `extension_column` | [343](../../../src/realizability.cpp#L343) |
| `euclidean_norm` | [354](../../../src/realizability.cpp#L354) |
| `projection_on_eigenvector` | [364](../../../src/realizability.cpp#L364) |
| `kernel_projection_tolerance` | [380](../../../src/realizability.cpp#L380) |
| `check_range_consistency` | [399](../../../src/realizability.cpp#L399) |
| `incorporate_range_consistency` | [425](../../../src/realizability.cpp#L425) |
| `numerical_rank` | [450](../../../src/realizability.cpp#L450) |
| `singular_as_given` | [460](../../../src/realizability.cpp#L460) |
| `finish_classification` | [466](../../../src/realizability.cpp#L466) |
| `validate_common_inputs` | [482](../../../src/realizability.cpp#L482) |
| `check_hamburger` | [517](../../../src/realizability.cpp#L517) |
| `check_stieltjes` | [571](../../../src/realizability.cpp#L571) |
| `check_hausdorff` | [627](../../../src/realizability.cpp#L627) |
| `TestedMoments::bound` | [698](../../../src/realizability.cpp#L698) |
| `raw_moments` | [704](../../../src/realizability.cpp#L704) |
| `standardized_moments` | [712](../../../src/realizability.cpp#L712) |
| `tested_moments` | [787](../../../src/realizability.cpp#L787) |
| `with_support_context` | [792](../../../src/realizability.cpp#L792) |
| `check_realizability` | [803](../../../src/realizability.cpp#L803), [832](../../../src/realizability.cpp#L832) |

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
| `cqmom_coordinate_scales` | [121](../../../src/reconstruction_scaling.cpp#L121) |
| `make_cqmom_scaling_context` | [141](../../../src/reconstruction_scaling.cpp#L141) |
| `restore_abscissas` | [178](../../../src/reconstruction_scaling.cpp#L178) |
| `restore_cqmom_abscissas` | [193](../../../src/reconstruction_scaling.cpp#L193) |
| `restore_bounds` | [209](../../../src/reconstruction_scaling.cpp#L209) |
| `can_restore_cqmom_result` | [224](../../../src/reconstruction_scaling.cpp#L224) |

### `src/recursive_quadrature.cpp`

[Source](../../../src/recursive_quadrature.cpp) · [Responsibilities, contracts,
and evidence](cqmom_reconstruction.md)

| Function or overload family | Definition lines |
|---|---|
| `create_ordered_root_layout` | [37](../../../src/recursive_quadrature.cpp#L37) |
| `append_unique_index` | [111](../../../src/recursive_quadrature.cpp#L111) |
| `local_branch_offset` | [131](../../../src/recursive_quadrature.cpp#L131) |
| `create_child_layout_case` | [160](../../../src/recursive_quadrature.cpp#L160) |
| `create_branch_moment_layout` | [208](../../../src/recursive_quadrature.cpp#L208) |
| `PreparedCqmomMetadata::PreparedCqmomMetadata` | [271](../../../src/recursive_quadrature.cpp#L271) |
| `create_prepared_cqmom_metadata` | [283](../../../src/recursive_quadrature.cpp#L283) |
| `prepared_cqmom_metadata_is_compatible` | [298](../../../src/recursive_quadrature.cpp#L298) |
| `extract_wheeler_moments` | [313](../../../src/recursive_quadrature.cpp#L313) |
| `child_layout_for_parent_node_count` | [332](../../../src/recursive_quadrature.cpp#L332) |
| `make_branch_rhs_matrix` | [346](../../../src/recursive_quadrature.cpp#L346) |
| `validate_positive_extents` | [375](../../../src/recursive_quadrature.cpp#L375) |
| `checked_multiply` | [387](../../../src/recursive_quadrature.cpp#L387) |
| `default_conditioning_order` | [397](../../../src/recursive_quadrature.cpp#L397) |
| `validate_conditioning_order` | [404](../../../src/recursive_quadrature.cpp#L404) |
| `ordered_extents` | [426](../../../src/recursive_quadrature.cpp#L426) |
| `from_ordered_index` | [437](../../../src/recursive_quadrature.cpp#L437) |
| `validate_abscissas` | [447](../../../src/recursive_quadrature.cpp#L447) |
| `valid_conditional_options` | [488](../../../src/recursive_quadrature.cpp#L488) |
| `valid_diagnostic_level` | [495](../../../src/recursive_quadrature.cpp#L495) |
| `records_success_branch_diagnostics` | [507](../../../src/recursive_quadrature.cpp#L507) |
| `create_ordered_root_slice` | [515](../../../src/recursive_quadrature.cpp#L515) |
| `one_dimensional_moments` | [523](../../../src/recursive_quadrature.cpp#L523) |
| `reconstruct_child_slices` | [528](../../../src/recursive_quadrature.cpp#L528) |
| `create_recursive_result` | [620](../../../src/recursive_quadrature.cpp#L620) |
| `record_failure` | [625](../../../src/recursive_quadrature.cpp#L625) |
| `validate_recursive_request_preflight` | [634](../../../src/recursive_quadrature.cpp#L634) |
| `resolve_cqmom_wheeler_options` | [658](../../../src/recursive_quadrature.cpp#L658) |
| `operator<` | [674](../../../src/recursive_quadrature.cpp#L674) |
| `path_has_prefix` | [699](../../../src/recursive_quadrature.cpp#L699) |
| `branch_limit_invalidated` | [715](../../../src/recursive_quadrature.cpp#L715) |
| `validate_active_branch_topology` | [727](../../../src/recursive_quadrature.cpp#L727) |
| `recursive_branch_node_count_limit` | [799](../../../src/recursive_quadrature.cpp#L799) |
| `compute_branch_quadrature` | [817](../../../src/recursive_quadrature.cpp#L817) |
| `validate_coordinate_support` | [823](../../../src/recursive_quadrature.cpp#L823) |
| `support_interval_for_dimension` | [850](../../../src/recursive_quadrature.cpp#L850) |
| `record_bounded_candidate_diagnostics` | [860](../../../src/recursive_quadrature.cpp#L860) |
| `branch_growth_trial_allowed` | [882](../../../src/recursive_quadrature.cpp#L882) |
| `product_extents` | [892](../../../src/recursive_quadrature.cpp#L892) |
| `flatten_multi_index` | [904](../../../src/recursive_quadrature.cpp#L904) |
| `unflatten_multi_index` | [924](../../../src/recursive_quadrature.cpp#L924) |
| `enumerate_multi_indices` | [1074](../../../src/recursive_quadrature.cpp#L1074) |
| `PhaseSpaceHierarchy::PhaseSpaceHierarchy` | [1086](../../../src/recursive_quadrature.cpp#L1086), [1093](../../../src/recursive_quadrature.cpp#L1093) |
| `PhaseSpaceHierarchy::node_count` | [1099](../../../src/recursive_quadrature.cpp#L1099) |
| `PhaseSpaceHierarchy::moment_extents` | [1108](../../../src/recursive_quadrature.cpp#L1108) |
| `PhaseSpaceHierarchy::flatten_node_index` | [1119](../../../src/recursive_quadrature.cpp#L1119) |
| `PhaseSpaceHierarchy::unflatten_node_index` | [1128](../../../src/recursive_quadrature.cpp#L1128) |
| `create_cqmom_moment_layout` | [1133](../../../src/recursive_quadrature.cpp#L1133) |
| `CqmomMomentData::CqmomMomentData` | [1195](../../../src/recursive_quadrature.cpp#L1195) |
| `CqmomMomentData::contains_required_index` | [1200](../../../src/recursive_quadrature.cpp#L1200) |
| `CqmomMomentData::moment` | [1205](../../../src/recursive_quadrature.cpp#L1205) |
| `extract_cqmom_moments` | [1221](../../../src/recursive_quadrature.cpp#L1221) |
| `create_cqmom_moment_data` | [1248](../../../src/recursive_quadrature.cpp#L1248) |
| `enumerate_quadrature_node_multi_indices` | [1270](../../../src/recursive_quadrature.cpp#L1270) |
| `rectangular_terminal_paths` | [1285](../../../src/recursive_quadrature.cpp#L1285) |
| `RecursiveQuadrature::RecursiveQuadrature` | [1320](../../../src/recursive_quadrature.cpp#L1320), [1327](../../../src/recursive_quadrature.cpp#L1327) |
| `RecursiveQuadrature::hierarchy` | [1338](../../../src/recursive_quadrature.cpp#L1338) |
| `RecursiveQuadrature::number_of_dimensions` | [1343](../../../src/recursive_quadrature.cpp#L1343) |
| `RecursiveQuadrature::total_nodes` | [1348](../../../src/recursive_quadrature.cpp#L1348) |
| `RecursiveQuadrature::weights` | [1353](../../../src/recursive_quadrature.cpp#L1353) |
| `RecursiveQuadrature::abscissas` | [1358](../../../src/recursive_quadrature.cpp#L1358) |
| `RecursiveQuadrature::adaptive_branch_counts_enabled` | [1363](../../../src/recursive_quadrature.cpp#L1363) |
| `RecursiveQuadrature::has_rectangular_storage` | [1368](../../../src/recursive_quadrature.cpp#L1368) |
| `RecursiveQuadrature::tree_nodes` | [1373](../../../src/recursive_quadrature.cpp#L1373) |
| `RecursiveQuadrature::terminal_node_ids` | [1378](../../../src/recursive_quadrature.cpp#L1378) |
| `RecursiveQuadrature::terminal_paths` | [1383](../../../src/recursive_quadrature.cpp#L1383) |
| `RecursiveQuadrature::weight` | [1388](../../../src/recursive_quadrature.cpp#L1388) |
| `RecursiveQuadrature::abscissa` | [1397](../../../src/recursive_quadrature.cpp#L1397) |
| `RecursiveQuadrature::node` | [1410](../../../src/recursive_quadrature.cpp#L1410) |
| `RecursiveQuadrature::set_weight` | [1420](../../../src/recursive_quadrature.cpp#L1420) |
| `RecursiveQuadrature::set_abscissa` | [1429](../../../src/recursive_quadrature.cpp#L1429) |
| `RecursiveQuadrature::set_adaptive_branch_counts_enabled` | [1442](../../../src/recursive_quadrature.cpp#L1442) |
| `RecursiveQuadrature::add_tree_node` | [1455](../../../src/recursive_quadrature.cpp#L1455) |
| `RecursiveQuadrature::set_terminal_quadrature` | [1467](../../../src/recursive_quadrature.cpp#L1467) |
| `compute_recursive_quadrature_from_validated_moment_data` | [1501](../../../src/recursive_quadrature.cpp#L1501) |
| `compute_recursive_quadrature_from_moment_data` | [1832](../../../src/recursive_quadrature.cpp#L1832) |
| `compute_recursive_quadrature_from_trusted_moment_data` | [1860](../../../src/recursive_quadrature.cpp#L1860) |
| `operator==` | [943](../../../src/recursive_quadrature.cpp#L943), [953](../../../src/recursive_quadrature.cpp#L953), [958](../../../src/recursive_quadrature.cpp#L958), [1064](../../../src/recursive_quadrature.cpp#L1064) |
| `operator!=` | [948](../../../src/recursive_quadrature.cpp#L948), [963](../../../src/recursive_quadrature.cpp#L963), [968](../../../src/recursive_quadrature.cpp#L968), [1069](../../../src/recursive_quadrature.cpp#L1069) |
| `MultiIndexTable::MultiIndexTable` | [973](../../../src/recursive_quadrature.cpp#L973), [989](../../../src/recursive_quadrature.cpp#L989) |
| `MultiIndexTable::at` | [1007](../../../src/recursive_quadrature.cpp#L1007) |
| `MultiIndexTable::push_back` | [1016](../../../src/recursive_quadrature.cpp#L1016), [1021](../../../src/recursive_quadrature.cpp#L1021) |
| `MultiIndexTable::pop_back` | [1034](../../../src/recursive_quadrature.cpp#L1034) |
| `MultiIndexTable::reserve` | [1043](../../../src/recursive_quadrature.cpp#L1043) |
| `MultiIndexTable::clear` | [1048](../../../src/recursive_quadrature.cpp#L1048) |
| `MultiIndexTable::to_multi_indices` | [1053](../../../src/recursive_quadrature.cpp#L1053) |

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
| `recurrence_is_finite` | [114](../../../src/wheeler.cpp#L114) |
| `has_negative_recurrence_b` | [119](../../../src/wheeler.cpp#L119) |
| `adaptive_sigma_tolerances` | [133](../../../src/wheeler.cpp#L133) |
| `compute_weights` | [161](../../../src/wheeler.cpp#L161) |
| `adaptive_accepts` | [173](../../../src/wheeler.cpp#L173) |
| `required_moments_are_finite` | [226](../../../src/wheeler.cpp#L226) |
| `attempt_quadrature_count` | [238](../../../src/wheeler.cpp#L238) |
| `compute_quadrature_at_scale` | [330](../../../src/wheeler.cpp#L330) |
| `compute_quadrature` | [403](../../../src/wheeler.cpp#L403) |
| `compute_conditional_quadrature` | [411](../../../src/wheeler.cpp#L411) |
| `select_bounded_quadrature` | [417](../../../src/wheeler.cpp#L417) |

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
