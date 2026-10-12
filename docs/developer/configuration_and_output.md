# Configuration and output

[`src/configuration.cpp`](../../src/configuration.cpp) reads `.qmom` files
into the same public cases and options a C++ caller builds, in namespace
`qmomkit::config`. [`src/output.cpp`](../../src/output.cpp) writes moment and
quadrature histories as text tables. Neither owns physics or a schedule:
callbacks, the time loop and the choice of when to write are the host's. Every
key and file format is listed in the user
[Reference](../user/7_reference.md#configuration-file); this page explains
the rules behind them.

## Parser structure

Every public parser takes text in three steps:

1. `parse_file` splits the text into a `ParsedFile`: scalar entries keyed by
   the qualified name `section.key`, `M[...]` entries of `[moments]` in a
   separate map, and the sections with their lines.
2. The parser's validation rejects unknown and ineffective keys and checks
   values, appending to `ParsedFile::errors`.
3. Typed parsing builds the public object, which is returned only if no error
   was recorded.

Errors are `ConfigurationError` records with a `ConfigurationErrorCode`, a
source name, a line and a message; both types are internal to
`configuration.cpp`. The public `ConfigurationResult<T>` carries either the
value or one error string, the first recorded error formatted as
`<source>:<line>: <message>`, without the line when it is unknown. Each
`load_*` function reads exactly the named file and passes its text to the
`parse_*` counterpart, with the path as source name.

## Syntax

The rules a user sees, such as comments, plain decimals and canonical indices,
are in the [Reference](../user/7_reference.md#configuration-file). The
parser also fixes these, which are equally part of the file format:

- `[name]` opens a section. A mixture-component or product-dimension section,
  such as `[initial.mixture.component.0]`, may appear once; other sections may
  be reopened. A malformed header closes the open section, so the keys after it
  are not attributed to the section above it.
- A line is `key = value`; neither may be empty, and a key appears once.
  Surrounding whitespace is trimmed; quotes are kept as part of the value.
- Numbers must match the plain decimal grammar of `is_decimal_number` and are
  then read in the classic locale, whatever the host's locale. Non-finite,
  overflowing and subnormal values are rejected.
- Integers are unsigned and booleans `true` or `false`. Lists are
  comma-separated, with no empty entry and no trailing comma;
  `[output] generalized_moments` separates its exponent vectors with `|`.
- Indexed sections are looked up by their canonical index, so `00` never
  matches `0`.

## Every key is checked

The case parsers, `parse_qmom`, `parse_cqmom` and `parse_cqmom_prepared`, and
the initial-condition parsers `parse_qmom_initial_condition` and
`parse_cqmom_initial_condition` validate the whole file, including
`[time_integration]` and `[output]`, though they return only part of it. The
option parsers check their own sections, the QMOM evolution options also the
`[qmom]` keys, and ignore the rest, so one file can serve all of them.

An unknown key or section is an error, and so is a key the selected settings
would not read:

- a parameter of another distribution family (`reject_other_family_parameters`);
- a distribution key or `number_density` at the level of a `type = mixture`,
  where each component sets its own (`parse_mixture`);
- an indexed section or key that the selected initialization does not read,
  such as a mixture component under `type = distribution`
  (`reject_unconsumed_indexed_entries`);
- `lower` or `upper` in a realizability section without
  `support = bounded_interval`;
- the per-coordinate `[time_integration.realizability.N]` in a QMOM file, the
  unnumbered `[time_integration.realizability]` in a CQMOM file, and a
  coordinate index beyond the phase-space rank;
- `M[...]` entries in a QMOM file.

In a realizability section, `support` is required whenever `check_mode` is not
`none`.

One exception is deliberate: explicit moments take precedence over `[initial]`.
When a QMOM file has `[moments] values`, or a CQMOM file has `[moments] extents`
or `M[...]` entries, every key and section under `initial` is accepted without
being read or checked.

Two keys are accepted although they can be without effect.
`[vandermonde] rcond` is stored for every solver but read only by
`moore_penrose`, since a loaded `CqmomCase` can still change its solver.
`[initial] maximum_orders` is read only by
`initialize_cqmom_moment_tensor(...)`; the case parsers generate just the
moments the layout requires.

## Mapping to objects

The number of entries in `[phase_space] node_counts` selects the formulation:
one is QMOM, more than one CQMOM. Each case parser rejects a file of the
other formulation.

| `parse_<name>` and `load_<name>` | Builds |
| --- | --- |
| `qmom` | `QmomCase` |
| `cqmom` | `CqmomCase` |
| `cqmom_prepared` | `CqmomPreparedCase` |
| `qmom_phase_space_evolution_options` | `QmomPhaseSpaceEvolutionOptions`: scheme, node options, realizability policy |
| `cqmom_phase_space_evolution_options` | `CqmomPhaseSpaceEvolutionOptions`: scheme, node options, `source_options.coordinate_realizability_policies` |
| `node_advection_options` | `NodeAdvectionOptions` |
| `realizability_policy` | `RealizabilityPolicy` |
| `output_options` | `OutputOptions` |

`parse_coordinate_realizability_policies`,
`parse_phase_space_evolution_scheme`, `parse_qmom_initial_condition` and
`parse_cqmom_initial_condition` have no `load_*` counterpart.

The case parsers go through the same steps a caller can take separately:
`initialize_qmom_moments(...)` or `initialize_cqmom_moments(...)` turn an
`InitialConditionSpecification` into moments, and `create_qmom_case(...)` or
`create_cqmom_case(...)` build the case. A QMOM case needs at least $2N$
moments at this point. A CQMOM explicit tensor must be complete within its
`extents`, and the case keeps only the moments its layout requires
([CQMOM](cqmom.md)). Options without a key, such as `QmomOptions::cutoff` or
the source-step retries, keep their C++ defaults.

`with_initial_moments(...)` replaces the moments of a `CqmomCase` or
`CqmomPreparedCase`, given as a `RawMomentTensor` or as `CqmomMomentData`. It
keeps every setting of the base case, and the prepared overloads share its
branch structure. It validates structure only, as a configuration error;
whether the moments can be reconstructed is left to `compute_cqmom(...)` and
`initial_state()`.

## Adding a key

A setting gets a key when it is numerical configuration a user chooses per
case, such as a node count, a tolerance or a scheme. Physics, callbacks,
schedules and implementation machinery get none, and a setting rarely changed
can stay C++-only.

1. Add the qualified key to the lists that accept it:
   - a case key to `qmom_explicit_allowed_keys`, which
     `qmom_initial_allowed_keys` extends, for QMOM, and to
     `cqmom_explicit_allowed_keys`, which `cqmom_initial_allowed_keys`
     extends, for CQMOM; a key both accept goes into both, and a `[qmom]` key
     also into the hard-coded check of `reject_unknown_qmom_option_keys`,
     which the QMOM evolution-option parser applies;
   - an `[initial]` key to `qmom_initial_allowed_keys` or
     `cqmom_initial_allowed_keys`, and a distribution parameter to
     `distribution_parameters` ([Distributions](distributions.md#adding-a-family));
   - a `[time_integration]` key to `node_advection_allowed_keys` or
     `realizability_policy_fields`; a top-level one, as `scheme` is, to
     `append_time_integration_allowed_keys` and to the hard-coded check
     of `reject_unknown_phase_space_evolution_keys`, which the option
     parsers apply;
   - an `[output]` key to `output_allowed_keys`.

   `with_runtime_allowed_keys` adds the time-integration and output keys to the
   case parsers, and the option parsers check against the same lists, apart
   from those two hard-coded checks. A new section also goes into
   `reject_unknown_qmom_sections`, `reject_unknown_cqmom_sections` or
   `time_integration_sections`, and, for a `[time_integration]` section, into
   the section checks of
   `reject_unknown_phase_space_evolution_sections` and
   `parse_realizability_policy`.
2. Parse it in the function that builds its object, record an error with the
   line for an invalid value, and reject it where the other settings make it
   ineffective.
3. Write it to the C++ field, leaving the field's default when the key is
   absent, so the file and the API share one default.
4. Add it to the [Reference](../user/7_reference.md).
5. Test the value, the default, an invalid value and an ineffective use in
   [`tests/unit/configuration.cpp`](../../tests/unit/configuration.cpp); a
   failure seen through a file belongs in
   [`tests/verification/configuration/`](../../tests/verification/configuration/).

## Output writers

`MomentWriter` and `QuadratureWriter`, declared in
[`include/qmomkit/output.hpp`](../../include/qmomkit/output.hpp), are built by
static factories. Each opens exactly the given path and writes the header; each
`write(...)` appends one time row, or one row per active node.

| Factory | Rows |
| --- | --- |
| `MomentWriter::for_qmom(path, moment_count)`, `(path, moment_count, output_options)`, `(path, selected_orders)` | `MomentSet`, optionally with generalized values or a `Quadrature1D` |
| `MomentWriter::for_qmom(path, state, output_options)` | `QmomState` with its `QmomResult` |
| `MomentWriter::for_cqmom(path, extents)`, `(path, extents, selected)`, `(path, extents, selected, orders)` | `RawMomentTensor`, optionally with generalized values |
| `MomentWriter::for_cqmom(path, layout)`, `(path, layout, orders)` | `CqmomMomentData`, optionally with generalized values |
| `MomentWriter::for_cqmom(path, state, output_options)` | `CqmomState` with its `CqmomQuadrature` |
| `MomentWriter::for_cqmom(path, request)` | `CqmomQuadrature`, whose requested moments it reconstructs |
| `QuadratureWriter::for_qmom(path, maximum_node_count)` | `Quadrature1D` or `QmomResult` |
| `QuadratureWriter::for_cqmom(path, dimension_count)` | `CqmomQuadrature` |

The factory fixes the columns, so a row of a type that belongs to another
factory throws `std::invalid_argument`; for a `MomentWriter`, the message
names the factory the row needs (`layout_mismatch_message`).

The format contract, which readers of the files rely on:

- The header is one line: `# time` and one label per column, each centered in
  its column.
- Values are written in scientific notation with `precision` digits after the
  point, a trailing argument of every factory with default `17`. Node and path
  indices are written as integers.
- Numbers and labels use the classic locale, whatever the host's locale.
- Every check of a row, or of the whole time block of a quadrature, runs
  before its first character is written, so a rejected write leaves the file
  unchanged. Non-finite values and mismatched sizes are rejected.
- Each write flushes before it returns. If the operating system refuses the
  data, the writer closes the file and throws `std::runtime_error`, and every
  later write throws; a file that cannot be opened throws at construction.

Ordinary moment columns come in a fixed order: ascending for QMOM, row-major
for a tensor, the layout's `required_indices`, the state's `moment_indices()`,
or the order given for selected or requested moments. Generalized moment
columns, labeled `G[...]`, follow them, one per exponent vector in
`OutputOptions::generalized_moment_orders` or the `orders` argument, each
vector with one finite non-negative exponent per coordinate. The state writers
evaluate them with `evaluate_generalized_moment(...)` from the quadrature
passed with the row; they never enter the state. The other writers take the
values from the caller, or for QMOM from a `Quadrature1D`, and reject a row
without them.

The writers read only `generalized_moment_orders` from `OutputOptions`.
`directory`, `write_moments`, `write_quadrature` and `interval` are carried to
the host, which decides where, what and when to write.

Existing columns keep their order and labels, since readers of the files rely
on them; a new quantity becomes a new column after them or a new factory. A new
writer documents its format here and in the
[Reference](../user/7_reference.md#output-files) before it becomes public API;
a format needed by only one example or test stays local to it.

[`tests/unit/output.cpp`](../../tests/unit/output.cpp) covers the formats,
the factory pairing and the failure paths, and
[`tests/verification/output/`](../../tests/verification/output/) writes the
moment and quadrature files of a QMOM and a CQMOM case and reads them back.
