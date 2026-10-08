# Output and supporting contracts

These utilities support numerical workflows without selecting physical models.
Most receive Level B/C treatment; shape/topology validators with meaningful
numerical boundaries are explained explicitly. The
[inventory](coverage_inventory.md) provides every definition/overload location.

## Moment and quadrature output — Level B

[output.cpp](../../../src/output.cpp) implements public `qmomkit::MomentWriter`
and `qmomkit::QuadratureWriter`, declared in
[output.hpp](../../../include/qmomkit/output.hpp). See
[output design](../output_design.md).

`MomentWriter::for_qmom` selects all orders or specified orders, or captures
the moments represented by a `QmomState` together with configured generalized
columns. `MomentWriter::for_cqmom` selects tensor extents, selected indices, a
required layout, a state with configured generalized columns, or a raw-moment
request evaluated from quadrature. Labels retain selected order and duplicates.
Constructors open/truncate the output stream in the classic locale
(`open_output_file`), validate writer layout, and write the header. They do not reconstruct a configured case.

`write_row` checks time and value shape/finiteness and emits
whitespace-separated fields at the writer's precision, 17 by default.
`write_aligned_row` expects values already aligned with the writer's selected
required/requested layout.
Public `write` overloads either select stored state or moment entries, evaluate
configured generalized values from a supplied quadrature, or reconstruct a
requested CQMOM moment from quadrature before writing. Each overload first
checks that the writer's layout matches its data argument; the map from
overload to accepting factory is one-to-one, and `layout_mismatch_message`
reports a mismatch by naming both the required and the actual factory. See the
[factory pairing table](../output_design.md#moment-writer-factory-pairing).
Every check of a moment row, or of the complete time block of a quadrature
write, runs before its first character is written, and values computed for a
row are computed before it as well, so a write rejected for its data leaves the
file unchanged. A write that is not stored throws after the flush, and the file
may then end in a partial row.

`QuadratureWriter::for_qmom` stores maximum node capacity; output may contain
fewer nodes. `for_cqmom` stores hierarchy. QMOM `write_qmom_row` validates
capacity, weight/coordinate sizes, and finite values. Overloads receiving a
QMOM result do not themselves reject its failure status before forwarding
numerical arrays. Therefore a diagnostic finite payload can be written; writing
is not certification of reconstruction success.

CQMOM `write` visits terminal nodes, emitting conditioning-path indices,
weight, and physical-coordinate values for each row. It validates shapes/path
rank and capacity, but does not independently infer or reconstruct the tree.
Fixed paths are aligned with physical row-major values under non-identity
conditioning order.

Evidence: [output unit tests](../../../tests/unit/output.cpp) compare headers
and rows, selected/raw requests, state-based generalized output, append
sequences, adaptive nodes, and malformed data. The inspected fixed-output case
uses identity conditioning. Tests whose expected rows copy the object's own
paths cannot independently establish path/value alignment.
`rejected_rows_leave_the_file_unchanged` checks that rejected writes leave the
file byte-identical, and `rows_that_cannot_be_stored_are_reported` that every
write that returns is in the file.

| Local helpers | Responsibility and semantics | Evidence |
|---|---|---|
| `qmom_label`, `cqmom_label`, `cqmom_abscissa_label` | Format order/coordinate metadata; labels do not impose physical units. | Output header tests |
| `sequential_indices`, `flatten_selected_indices` | Build all-order selection or map checked multi-indices to flat offsets in request order. | Selected-output tests |
| `qmom_labels`, `cqmom_labels`, `all_cqmom_labels`, `qmom_quadrature_labels`, `cqmom_quadrature_labels` | Construct header fields matching the selected writer layout. | Output header/row tests |
| `open_output_file`, `flush_or_throw` | Open the stream in the classic locale, so the global locale changes neither decimal point nor index digits; flush each row or quadrature block and throw if it cannot be stored, so no accepted row is lost in the buffer. | `writers_ignore_the_global_locale`, `rows_that_cannot_be_stored_are_reported` in output unit tests |
| Writer constructors, `write_header`, `layout_mismatch_message`, forwarding `write` overloads | Own stream/layout and adapt supported input representations. Opening truncates the file before the constructor's later checks. | Output failure tests |

## Shared finite and shape validation — Levels A/C

Definitions in
[cqmom_validation.hpp](../../../src/internal/cqmom_validation.hpp) are in
`qmomkit::internal` and are implementation-only.

`validate_fixed_cqmom_moment_layout` constructs the canonical layout and
delegates to `validate_cqmom_moment_layout_matches_expected`, which requires
exact hierarchy counts, conditioning order, and required-index order. Equal
sets in a different order are insufficient because the value vector is
positional. `validate_fixed_cqmom_moment_data` additionally checks matching
value count and finiteness via `validate_cqmom_moment_values_match_layout`.
`validate_cqmom_moment_data_matches_hierarchy` compares against the caller's
hierarchy, not only the layout's self-description. Trusted callers deliberately
reuse established layout invariants and can call the narrower value check.

`raw_moment_tensor_value_count` returns no count for zero
extents or multiplication overflow. `multi_index_within_extents` checks rank
and each bound; `same_cqmom_layout` compares counts/order/index enumeration.
These checks do not test realizability.

`validate_cqmom_quadrature_storage` checks nonempty weights, dimension rows,
per-row terminal count, and all finite weights/coordinates. It does **not**
require positive weights, path/value correspondence, or moment recovery.
Callers add operation-specific path/hierarchy checks where necessary.
`validate_cqmom_moment_rate` checks rate count/finiteness; it does not enforce
a model-specific conservation law.

Evidence: malformed layout, tensor, quadrature, and rate assertions in [public
API](../../../tests/unit/public_api.cpp), [moment
reconstruction](../../../tests/unit/moment_reconstruction.cpp),
[source](../../../tests/unit/moment_source.cpp). These test the relevant public
boundaries; they do not mean every trusted internal call revalidates
everything.

[numeric_utils.hpp](../../../src/internal/numeric_utils.hpp) supplies
`all_finite` (linear scan).
[phase_space_velocity_validation.hpp](../../../src/internal/phase_space_velocity_validation.hpp)
supplies `format_diagnostic_vector` (stable readable vector text) and
`validate_phase_space_velocity_result` plus message helpers (exact dimension
and finite velocity). They throw/contextualize invalid callback output; they
neither clip velocity nor infer coordinate semantics.

## Performance counters — Level B/C

[performance.hpp](../../../src/internal/performance.hpp) defines
`qmomkit::internal::performance` workload counters. With
`QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION`, `CounterSession::snapshot` copies
counters and `reset` clears them. `CounterScope` installs a thread-local active
session and its destructor restores the previous session, including nested
scopes. `active_counters` returns a pointer or null. This is thread-local
observation, not shared mutable numerical algorithm state.

Each `record_*` helper increments only the active session. Most record one
event; `record_vandermonde_solve` also classifies the selected solve kind, and
`record_embedded_rk_step` separates accepted/rejected steps. These count
instrumented call sites, not a mathematically complete operation count.
Disabled definitions are no-ops; `enabled` is false and snapshots are zero.
Both preprocessor branches were inspected; timing is not used as correctness
evidence.

Evidence: [instrumentation tests](../../../tests/unit/performance_instrumentation.cpp)
check counters, scopes, diagnostic separation, and numerical equivalence;
[instrumentation design](../performance_instrumentation.md) explains intended use.
No benchmark measurements or broad performance conclusions were produced here.
