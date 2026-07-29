# Configuration Design

This page documents the current QMOMKit design for configuration handling. It
is a developer reference for maintainers, contributors, and AI-assisted
development workflows.

Configuration files describe reproducible numerical setups. They do not replace
the C++ API, own a simulation loop, define application physics, or execute
QMOMKit algorithms by themselves.

QMOMKit case configuration files use the user-facing `.qmom` extension.

For raw-moment, tensor, phase-space, conditioning-order, distribution, node, and
weight conventions, see `numerical_conventions.md`. For transport behavior, see
`transport_design.md`.

## Configuration Philosophy

QMOMKit configuration is a text form of the same numerical setup that users can
build through the C++ API. It exists so examples, verification drivers, Monte
Carlo workflows, and other projects that use QMOMKit can load repeatable setup
data from text files without recompiling C++.

The core design principle is separation:

- configuration files describe numerical setup data and selected policies;
- configuration loading validates that data and converts it to public case,
  request, or option objects;
- runtime algorithms consume those case, request, or option objects exactly as
  they do for programmatically constructed inputs.

A valid configuration is not a guarantee that a numerical reconstruction or
transport step will succeed. Configuration validation checks structure,
supported values, and basic consistency before execution. Algorithmic statuses
from QMOM, CQMOM, realizability checks, node advection, source updates, and
transport remain part of the runtime numerical algorithms.

Configuration is intentionally application agnostic. It does not define meshes,
boundary conditions, CFD operators, source callbacks, velocity callbacks, MPI
communication, units, or physical models.

## Configuration Architecture

The conceptual flow is:

```text
.qmom file
      |
      v
configuration loading
      |
      v
validated configuration objects
      |
      v
QMOM cases, CQMOM cases, and options
      |
      v
runtime algorithms
```

The public configuration API supports both file-based and string-based entry
points. Loading from disk reads exactly the requested file, then delegates to
the same parsing path used by string inputs.

The main parts of configuration handling are:

- parsing: reads the line-oriented `.qmom` text format into key/value entries,
  indexed tensor-moment entries, section names, source names, and line numbers;
- validation: generally rejects malformed syntax, duplicate entries, unknown
  sections, unknown keys, invalid values, missing required data, and
  inconsistent configuration combinations;
- initialization: turns explicit moments or analytical distribution
  specifications into `MomentSet`, CQMOM required moments, or explicit dense
  representation tensors;
- case construction: builds `QmomCase`, `CqmomCase`, or explicit
  representation `CqmomLeanCase` and `CqmomDenseCase` objects from initialized
  moments, phase-space hierarchy, and options;
- policy extraction: parses optional transport scheme, transport options,
  output options, and runtime realizability policy objects for drivers that
  choose to use them.

The parser's private data structures should stay hidden from users and
extension code. Successful public configuration calls return populated public
objects through `value()`, while failed calls return a formatted diagnostic
string through `error()` and no value.

## QMOM And CQMOM Mapping

Configuration maps user-facing setup concepts to the same public object types
that direct C++ callers use.

Configuration files define the phase-space problem. Phase-space structure
selects the configured formulation: one phase-space dimension maps to QMOM,
while multiple phase-space dimensions map to CQMOM. Coordinate names are
metadata and must not be used for formulation selection. Future quadrature
families should add explicit, documented selection rules when they are
implemented instead of overloading this phase-space rule.

`format_version` is optional. When omitted, configuration loading assumes the
current format version. When present, it is validated so future or unsupported
formats fail explicitly rather than being silently reinterpreted.

### QMOM

For one-dimensional QMOM, configuration maps to:

```text
MomentSet
node_count
QmomOptions
QmomCase
```

The moment input can be explicit raw moments or generated moments from a
one-dimensional analytical distribution or mixture. Explicit moments are stored
in ascending order:

```text
m0, m1, m2, ...
```

Configuration loading validates that the requested node count is positive and
that at least `2 * node_count` moments are available before creating a
`QmomCase`. Wheeler inversion success is still checked later by
`compute_qmom(...)`.

QMOM configuration requires exactly one phase-space node count. A file with
multiple node counts is not a QMOM case.

The `[qmom] adaptive` setting maps to `QmomOptions::adaptive`. It controls the
one-dimensional inversion behavior used by both ordinary QMOM and branch-local
CQMOM inversions.

### CQMOM

For multidimensional CQMOM, configuration maps to:

```text
PhaseSpaceHierarchy
QmomOptions
adaptive_branch_counts
CqmomDiagnosticLevel
VandermondeSolver
CqmomMomentStorage
CqmomCase
```

This is the recommended user-facing loading path exposed by
`load_cqmom(...)`. The returned `CqmomCase` stores the selected dense or lean
moment representation and reconstruction policies while avoiding any distinction
between analytical and explicit moment input in ordinary use. If
`cqmom.moment_storage` is omitted, the neutral parser selects lean storage.

The explicit lean representation workflow maps parsed CQMOM setup data to:

```text
InitialConditionSpecification
PhaseSpaceHierarchy
CqmomMomentLayout
RawMomentRequest
CqmomMomentData
```

For analytical independent-product initialization, configuration should generate
only the moments required by `create_cqmom_moment_layout(...)` and adapts the
selected result to `CqmomMomentData`. This path does not construct a dense
`RawMomentTensor` and does not require `initial.maximum_orders`.

The explicit dense representation workflow maps parsed CQMOM setup data to:

```text
RawMomentTensor
CqmomDenseCase
```

This path is exposed through `load_cqmom_dense(...)` for callers that
intentionally choose full rectangular tensor moment storage for full tensor
access, external data, validation, debugging, or interoperability. Analytical
independent-product configuration on this path requires
`initial.maximum_orders`, which defines the rectangular dense tensor extents.
`load_cqmom_lean(...)` exposes the explicit lean representation workflow and
returns `CqmomLeanCase` directly for callers that intentionally manage the
required CQMOM moment representation.

Phase-space node counts are attached to phase-space dimensions, not to
conditioning levels. If `conditioning_order` is omitted, identity order is used.
If it is provided, it must be a permutation of the phase-space dimensions.
CQMOM configuration requires at least two phase-space dimensions; a
one-dimensional phase space is interpreted as QMOM.

Explicit tensor moments are indexed in phase-space dimension order:

```text
M[i0,i1,...] = m_[i0,i1,...]
```

For explicit dense representation loading, the declared tensor extents define the
complete expected moment tensor. Configuration loading validates rank, index
range, completeness, and that each extent is large enough for the requested node
count in that dimension before creating a `CqmomDenseCase`. Recursive CQMOM
success is still checked later by `compute_cqmom(...)`.

The settings `phase_space.adaptive_branch_counts` and `[qmom] adaptive` are
parsed independently. `phase_space.adaptive_branch_counts` controls whether
CQMOM may accept branch-dependent quadrature counts and, when enabled,
automatically enables adaptive Wheeler for all CQMOM marginal and conditional
inversions. `[qmom] adaptive` still controls standalone one-dimensional QMOM
behavior.

The optional `[cqmom] diagnostic_level = detailed|basic|none` setting maps to
the CQMOM diagnostic policy stored in the loaded case. If omitted,
configuration preserves the public default `CqmomDiagnosticLevel::detailed`.

## Distribution Initialization Design

Configuration can reference analytical distributions as a compact way to define
initial moments. Distribution specifications are not runtime probability
objects consumed directly by QMOM or CQMOM. They are used during initialization
to generate raw moments of the number density function.

The one-dimensional paths are:

```text
distribution + number_density + maximum_order
      -> MomentSet

mixture components + component number densities + maximum_order
      -> MomentSet
```

The multidimensional analytical paths are:

```text
independent product marginals + global number_density + maximum_orders
      -> RawMomentTensor
```

for the explicit dense representation workflow, where `maximum_orders` is
required, and:

```text
independent product marginals + global number_density + CQMOM hierarchy
      -> CqmomMomentData
```

for the explicit lean representation workflow. The lean path derives the
requested raw moment multi-indices from the hierarchy and conditioning order
instead of from rectangular tensor extents. In this path, `maximum_orders` is
optional dense metadata and is not used to select moments.

Mixture component densities are absolute number-density contributions. An
independent-product initialization has one global number density; per-dimension
number densities are not part of the CQMOM configuration model.

Generated moments follow the conventions in `numerical_conventions.md`: `m0` is
the total number density, moments are not normalized away, and multidimensional
tensor entries follow phase-space dimension order.

## Transport And Policy Configuration

Configuration may include optional runtime settings that are useful to drivers:

- transport scheme selection;
- node-advection scheme selection and RK23 controls;
- CQMOM source-stage active branch-count policy;
- output options;
- runtime realizability policy.

These settings are parsed into public objects or enum values. They do not
install callbacks, run transport, write output, or define a physical model.
Case loaders validate recognized runtime settings when those settings are
present in the same `.qmom` file, but they still return only `QmomCase` or
`CqmomCase`.

| Settings | Validated by case loading | Stored in case/state | Must be loaded separately |
| --- | --- | --- | --- |
| QMOM/CQMOM inversion settings | yes | yes | no |
| Vandermonde and CQMOM storage settings | yes | yes | no |
| `time_integration.scheme` | yes | no | yes |
| `time_integration.node_advection.*` | yes | no | yes |
| `time_integration.realizability.*` | yes | no | yes |

Transport scheme configuration maps to `TransportScheme`. If omitted, the
parsed scheme defaults to `TransportScheme::strang`.

Node-advection configuration maps to `NodeAdvectionOptions`. The scheme lives
under `[time_integration.node_advection]` so it is not confused with
`[time_integration] scheme`, which controls operator splitting. RK23 controls
live under `[time_integration.node_advection.rk23]` and map directly to
`AdaptiveRk23Options`:

```text
relative_tolerance
absolute_tolerance
initial_step
max_step
min_step
max_steps
```

CQMOM source-stage active-count evolution is automatic. Continuous
reconstruction constrains each existing adaptive branch to its current active
count. After an accepted discontinuous source update, adaptive CQMOM passes the
current active topology and enables the internal one-node growth trial. New
descendants whose parent topology changed start from configured maximum counts
and adapt downward.

The parser default-constructs `NodeAdvectionOptions` and only overwrites keys
that are present, so configuration parsing does not duplicate
`AdaptiveRk23Options` defaults. If omitted, the parsed node-advection scheme
defaults to `NodeAdvectionScheme::ssp_rk3`. Applications can assign the parsed
object directly to transport options. This configuration choice affects only
continuous quadrature-node trajectories; it does not configure moment-source
integration, realizability policies, or transport splitting.

Runtime realizability configuration maps to `RealizabilityPolicy`. A parsed
policy is inert until a driver wires it into an acceptance predicate or another
runtime check. Support is explicit; it is not inferred from distribution family,
phase-space name, or coordinate label.

Output configuration maps to `OutputOptions`. Numerical algorithms do not write
files because output settings are present. Drivers choose whether and how to
connect those options to output helpers.

## Validation And Error Handling

Configuration APIs return:

```cpp
ConfigurationResult<T>
```

where `value()` is valid only on success and `error()` contains the public
diagnostic string on failure. Internally, parser failures are still assembled
from `ConfigurationError` records carrying:

- an error code;
- a source name;
- a line number when known;
- a human-readable message.

Configuration errors should be detected before numerical execution when the
information is available from the file. Examples include malformed syntax,
unknown sections, duplicate keys, missing required keys, invalid booleans,
non-finite values, unsupported distribution families, invalid conditioning
orders, incomplete tensor moment sets, insufficient moment counts, and invalid
output or policy settings.

Configuration validation should not silently reinterpret invalid input. A new
key or section must either be documented, parsed, validated, and tested, or be
reported as unknown.

Explicit moment initialization is the intentional precedence rule. When
explicit moments are provided, they are the source of truth for the initialized
moment representation. Distribution-based `[initial]` entries that would
otherwise generate moments are unused in that path, and the current
configuration behavior may ignore those unused initialization entries instead of
validating them as active inputs.

Runtime numerical failure remains separate. A configuration can be syntactically
and structurally valid while the configured moment sequence is not accepted by
QMOM or CQMOM reconstruction under the selected options. Callers must check both
configuration results and algorithm results.

## Public API Boundary

The public user-facing configuration API lives in:

```text
include/qmomkit/configuration.hpp
```

Public configuration capabilities include:

- `parse_qmom(...)` and `load_qmom(...)`;
- `parse_cqmom_dense(...)` and `load_cqmom_dense(...)`;
- `parse_cqmom_lean(...)` and `load_cqmom_lean(...)`;
- `parse_qmom_initial_condition(...)` and
  `parse_cqmom_initial_condition(...)`;
- `initialize_moments(...)` and `initialize_cqmom_dense_moments(...)`;
- `initialize_cqmom_moment_data(...)`;
- `create_qmom_case(...)`, `create_cqmom_dense_case(...)`, and
  `create_cqmom_lean_case(...)`;
- `parse_transport_scheme(...)`;
- `parse_node_advection_options(...)` and
  `load_node_advection_options(...)`;
- `parse_node_advection_scheme(...)` and
  `load_node_advection_scheme(...)` compatibility helpers;
- `parse_qmom_transport_options(...)`, `load_qmom_transport_options(...)`,
  `parse_cqmom_transport_options(...)`, and
  `load_cqmom_transport_options(...)`;
- `parse_realizability_policy(...)` and `load_realizability_policy(...)`;
- `parse_output_options(...)` and `load_output_options(...)`;
- `ConfigurationResult<T>`, `ConfigurationError`, and
  `ConfigurationErrorCode`.

Implementation details include:

- parser storage structures;
- line tokenization;
- section and key qualification;
- duplicate-entry tracking;
- unknown-key filtering;
- concrete schema lists;
- file reading code;
- helper functions that convert text tokens to typed values.

These details support the public API but are not extension points for projects
that use QMOMKit.

## Extension Principles

Configuration changes should preserve the existing separation between syntax and
runtime behavior.

Stable rules for adding configuration capabilities are:

- map configuration settings to existing public request, option, policy, or
  output objects when possible;
- keep parser internals behind `configuration.cpp`;
- validate all user input before constructing runtime objects;
- reject unknown active keys and sections instead of ignoring misspellings;
- keep phase-space dimensions, tensor ordering, conditioning order, and moment
  normalization consistent with `numerical_conventions.md`;
- do not infer physics from names, units, distribution families, or output
  labels;
- avoid coupling configuration to a single application, CFD framework, example,
  or validation driver;
- add focused tests for parsing, validation failures, default behavior, and the
  constructed runtime objects.

Configuration should remain a reproducible description of QMOMKit numerical
workflows, not an application-level solver or a hidden control plane for
application logic.
