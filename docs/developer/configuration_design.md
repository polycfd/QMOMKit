# Configuration design

This page documents the current QMOMKit design for configuration handling.

Configuration files describe reproducible numerical setups. They do not replace
the C++ API, own a simulation loop, define application physics, or execute
QMOMKit algorithms by themselves.

QMOMKit case configuration files use the user-facing `.qmom` extension.

For raw-moment, tensor, phase-space, conditioning-order, distribution, node,
and weight conventions, see
[numerical_conventions.md](numerical_conventions.md). For phase-space evolution
behavior, see
[phase_space_evolution_design.md](phase_space_evolution_design.md).

## Contents

- [Configuration Philosophy](#configuration-philosophy)
- [Configuration Architecture](#configuration-architecture)
- [QMOM And CQMOM Mapping](#qmom-and-cqmom-mapping)
- [Distribution Initialization Design](#distribution-initialization-design)
- [Phase-Space Evolution and Policy Configuration](#phase-space-evolution-and-policy-configuration)
- [Validation And Error Handling](#validation-and-error-handling)
- [Public API Boundary](#public-api-boundary)
- [Extension Principles](#extension-principles)

## Configuration philosophy

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
phase-space evolution step will succeed. Configuration validation checks
structure, supported values, and basic consistency before execution.
Algorithmic statuses from QMOM, CQMOM, realizability checks, node advection,
source updates, and phase-space evolution remain part of the runtime numerical
algorithms.

Configuration is intentionally application agnostic. It does not define meshes,
boundary conditions, CFD operators, source callbacks, velocity callbacks, MPI
communication, units, or physical models.

The `.qmom` format is deliberately narrower than the public C++ API. It maps
the supported case and runtime configuration surface; advanced lower-level
controls, such as reconstruction thresholds and source-integration acceptance
callbacks, remain API-only unless they acquire a documented configuration
contract.

## Configuration architecture

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
  specifications into `MomentSet`, CQMOM required moments, or explicit full
  representation tensors;
- case construction: builds `QmomCase`, `CqmomPreparedCase`, or an explicit
  `CqmomCase` from initialized moments, phase-space hierarchy, and options;
- policy extraction: parses an optional splitting scheme, evolution options,
  output options, and runtime realizability policy objects for drivers that
  choose to use them.

The parser's private data structures should stay hidden from users and
extension code. Successful public configuration calls return populated public
objects through `value()`, while failed calls return a formatted diagnostic
string through `error()` and no value.

## QMOM and CQMOM mapping

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
optional AbscissaBounds
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

Optional paired `[phase_space]` `support_lower` and `support_upper` entries
map to `QmomCase::abscissa_bounds`. QMOM requires exactly one entry in each
list. Both keys may be omitted, or both entries may be `none`, for unbounded
reconstruction. Numeric entries must be finite and satisfy `lower < upper`.
This reconstruction metadata is distinct from the bounded-interval
`MomentSupport` configuration used by realizability policies.

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
CqmomPreparedCase
```

This is the recommended user-facing loading path exposed by
`load_cqmom_prepared(...)`. The returned `CqmomPreparedCase` stores the moment
representation and reconstruction policies while avoiding any distinction
between analytical and explicit moment input in ordinary use.

Both CQMOM loaders initialize the moments the same way, mapping parsed CQMOM
setup data to:

```text
InitialConditionSpecification
PhaseSpaceHierarchy
CqmomMomentLayout
RawMomentRequest
CqmomMomentData
```

For analytical independent-product initialization, configuration requests
exactly `CqmomMomentLayout::required_indices` from
`create_cqmom_moment_layout(...)` and adapts the selected result to
`CqmomMomentData`. This path does not construct a dense
`RawMomentTensor`; `initial.maximum_orders` has no effect on the requested
moment layout.

Analytical independent-product configuration uses explicit
`initial.maximum_orders` when present; when omitted, the tensor extents come
from `PhaseSpaceHierarchy::moment_extents()`. The derived extent in dimension
`d` is `2 * node_counts[d]`, so the corresponding maximum raw order is
`2 * node_counts[d] - 1`; explicit orders may be larger but not smaller.
Reconstruction consumes only the `CqmomMomentLayout` required by the hierarchy
and conditioning order.

`load_cqmom(...)` returns the same configuration as an editable `CqmomCase`
definition, for callers that change settings in code before preparing it;
`load_cqmom_prepared(...)` returns that definition already prepared.

Phase-space node counts are attached to phase-space dimensions, not to
conditioning levels. If `conditioning_order` is omitted, identity order is
used. If it is provided, it must be a permutation of the phase-space
dimensions. CQMOM configuration requires at least two phase-space dimensions; a
one-dimensional phase space is interpreted as QMOM.

Explicit tensor moments are indexed in phase-space dimension order:

```text
M[i0,i1,...] = m_[i0,i1,...]
```

When moments are declared as an explicit tensor, the declared extents define
the complete expected moment tensor. Configuration loading validates rank,
index range and completeness, and selecting the required moments from it
rejects a tensor too small for the requested node counts.
Recursive CQMOM success is still checked later by `compute_cqmom(...)`.

`phase_space.adaptive_branch_counts` controls whether CQMOM may accept
branch-dependent quadrature counts and, when enabled, automatically enables
adaptive Wheeler for all CQMOM marginal and conditional inversions. A CQMOM
configuration with `[qmom] adaptive = true` must enable adaptive branch counts;
fixed-count CQMOM rejects that combination during configuration validation.
`[qmom] adaptive` remains the standalone one-dimensional QMOM control.

The optional `[cqmom] diagnostic_level = detailed|basic|none` setting maps to
the CQMOM diagnostic policy stored in the loaded case. If omitted,
configuration preserves the public default `CqmomDiagnosticLevel::none`.

### Immutable explicit-moment replacement

`config::with_initial_moments(const CqmomPreparedCase&, RawMomentTensor)` is
the programmatic boundary for replacing the explicit initial tensor of an
existing case. It keeps the base case's settings and prepared metadata and
replaces only its initial moments, deriving the fixed layout from the base
hierarchy and extracting its required values. The base case is unchanged. So is
the hierarchy, so every case derived from one prepared case shares one recursive
branch structure.

The overloads taking a `CqmomCase` hold the validation and the replacement. The
prepared overloads call them and then attach the base case's prepared metadata,
so both kinds of base case report the same errors.

This boundary validates tensor rank, positive non-overflowing extents, value
count, finiteness, and reconstruction coverage without parser
metadata. It deliberately does not invoke CQMOM inversion. Consequently,
programmatic and parsed explicit tensors share the same lifecycle: structural
errors are configuration failures, while non-realizability is reported by
`compute_cqmom(...)` or by `initial_state()`.

`config::with_initial_moments(const CqmomPreparedCase&, CqmomMomentData)` and
its `CqmomCase` counterpart take the required moment set instead of a tensor.
A caller whose moments already arrive as the required CQMOM moment set — from
`compute_raw_moments(...)` with the layout's `required_indices` adapted by
`create_cqmom_moment_data(...)`, or from any other producer that evaluates
exactly those indices — replaces the initial moments directly instead of
materializing a dense tensor whose additional entries the case never reads. For
a hierarchy whose dense tensor is much larger than its required set, that is the
difference between generating the required moments and generating the full
rectangle.

The overload taking the required moment set validates that the
supplied layout is the required CQMOM moment layout of the base hierarchy, that
there is one value per required index, and that the values are finite. It leaves
realizability to reconstruction exactly as the tensor overload does, and both
overloads produce the same case for the same moments.

## Distribution initialization design

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

for callers that build a full raw moment tensor themselves, where
`maximum_orders` sets the retained tensor extent, including when a sample-based
initialization evaluates raw moments from a sample population. Configured CQMOM
loading does not consult it. Its own path is:

```text
independent product marginals + global number_density + CQMOM hierarchy
      -> CqmomMomentData
```

for both CQMOM loaders. This path derives the requested raw moment
multi-indices from the hierarchy and conditioning order instead of from
rectangular tensor extents. In this path, `maximum_orders` is
optional tensor metadata and is not used to select moments.

Mixture component densities are absolute number-density contributions. An
independent-product initialization has one global number density; per-dimension
number densities are not part of the CQMOM configuration model.

Generated moments follow the conventions in
[numerical_conventions.md](numerical_conventions.md): $m_0$ is the total number
density, moments are not normalized away, and multidimensional tensor entries
follow phase-space dimension order.

## Phase-space evolution and policy configuration

Configuration may include optional runtime settings that are useful to drivers:

- phase-space evolution scheme selection;
- node-advection scheme selection and embedded-RK controls;
- output options;
- runtime realizability policy.

These settings are parsed into public objects or enum values. They do not
install callbacks, run phase-space evolution, write output, or define a
physical model. Case loaders validate recognized runtime settings when those
settings are present in the same `.qmom` file, but they still return only
`QmomCase` or `CqmomPreparedCase`.

| Settings | Validated by case loading | Stored in case/state | Must be loaded separately |
| --- | --- | --- | --- |
| QMOM/CQMOM inversion settings | yes | yes | no |
| Vandermonde and CQMOM storage settings | yes | yes | no |
| `time_integration.scheme` | yes | no | yes |
| `time_integration.node_advection.*` | yes | no | yes |
| `time_integration.realizability.*` | yes | no | yes |

Phase-space evolution scheme configuration maps to `PhaseSpaceEvolutionScheme`.
If omitted, the parsed scheme defaults to `PhaseSpaceEvolutionScheme::strang`.

Node-advection configuration maps to `NodeAdvectionOptions`. The scheme lives
under `[time_integration.node_advection]` so it is not confused with
`[time_integration] scheme`, which controls the operator-splitting scheme. The
available advection choices are `rk23` (the default Bogacki-Shampine 3(2)
pair), `rk45` (Dormand-Prince RK5(4)7M from Table 2), and the separate
fixed-step `ssp_rk3` path. Embedded-RK controls live under
`[time_integration.node_advection.adaptive]` and map directly to
`AdaptiveRkOptions`:

```text
relative_tolerance
absolute_tolerance
initial_step
max_step
min_step
max_steps
```

Drivers parse or load the complete `NodeAdvectionOptions` object with
`parse_node_advection_options(...)` or `load_node_advection_options(...)`.
When only the integration scheme is needed, they inspect its `.scheme` member.

CQMOM active-count evolution during moment-source stages is automatic.
Advection reconstruction constrains each existing adaptive branch to its
current active count. After an accepted moment-source update, adaptive CQMOM
passes the current active topology and enables the internal one-node growth
trial. New descendants whose parent topology changed start from configured
maximum counts and adapt downward.

The parser default-constructs `NodeAdvectionOptions` and only overwrites keys
that are present, so configuration parsing does not duplicate
`AdaptiveRkOptions` defaults. If omitted, the parsed node-advection scheme
defaults to `NodeAdvectionScheme::rk23`. Applications can assign the parsed
object directly to the evolution options. This configuration choice affects
only advection quadrature-node trajectories; it does not configure
moment-source integration, realizability policies, or operator splitting.

`[time_integration.realizability]` configures QMOM realizability checks during
moment-source SSP-RK3 evolution. For CQMOM,
`[time_integration.realizability.<coordinate>]` configures realizability
checks for individual phase-space coordinates; the zero-based coordinate index
must be smaller than the configured phase-space rank. When evolution
options are loaded from the configuration file, the corresponding QMOM or
CQMOM phase-space evolution path consumes these settings during its
moment-source stage. Support is explicit; it is not inferred from distribution
family, phase-space name, or coordinate label.

Output configuration maps to `OutputOptions`. Numerical algorithms do not write
files because output settings are present. Drivers choose whether and how to
connect those options to output helpers.

## Validation and error handling

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
configuration behavior may ignore those unused initialization entries instead
of validating them as active inputs.

Runtime numerical failure remains separate. A configuration can be
syntactically and structurally valid while the configured moment sequence is
not accepted by QMOM or CQMOM reconstruction under the selected options.
Callers must check both configuration results and algorithm results.

## Public API boundary

The public user-facing configuration API lives in:

```text
include/qmomkit/configuration.hpp
```

Public configuration capabilities include:

- `parse_qmom(...)` and `load_qmom(...)`;
- `parse_cqmom(...)` and `load_cqmom(...)`;
- `parse_qmom_initial_condition(...)` and
  `parse_cqmom_initial_condition(...)`;
- `initialize_qmom_moments(...)` and `initialize_cqmom_moments(...)`;
- `initialize_cqmom_moment_tensor(...)`;
- `create_qmom_case(...)` and `create_cqmom_case(...)`;
- `with_initial_moments(...)` for a `RawMomentTensor` or `CqmomMomentData`
  replacement;
- `parse_phase_space_evolution_scheme(...)`;
- `parse_node_advection_options(...)` and
  `load_node_advection_options(...)`;
- `parse_qmom_phase_space_evolution_options(...)`, `load_qmom_phase_space_evolution_options(...)`,
  `parse_cqmom_phase_space_evolution_options(...)`, and
  `load_cqmom_phase_space_evolution_options(...)`;
- `parse_realizability_policy(...)`, `load_realizability_policy(...)`, and
  `parse_coordinate_realizability_policies(...)`;
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

## Extension principles

Configuration changes should preserve the existing separation between syntax
and runtime behavior.

Stable rules for adding configuration capabilities are:

- map configuration settings to existing public request, option, policy, or
  output objects when possible;
- keep parser internals behind `configuration.cpp`;
- validate all user input before constructing runtime objects;
- reject unknown active keys and sections instead of ignoring misspellings;
- keep phase-space dimensions, tensor ordering, conditioning order, and moment
  normalization consistent with [numerical_conventions.md](numerical_conventions.md);
- do not infer physics from names, units, distribution families, or output
  labels;
- avoid coupling configuration to a single application, CFD framework, example,
  or validation driver;
- add focused tests for parsing, validation failures, default behavior, and the
  constructed runtime objects.

Configuration should remain a reproducible description of QMOMKit numerical
workflows, not an application-level solver or a hidden control plane for
application logic.
