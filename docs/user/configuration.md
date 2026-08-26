# Configuring QMOMKit Workflows

QMOMKit configuration files describe the repeatable numerical part of a
moment-method workflow in a small text format. They are optional: applications
can also create and configure QMOMKit cases directly through the C++ API.

## Configuration Files And Direct C++ Setup

`.qmom` files provide a convenient way to describe reproducible cases. They are
especially useful for examples, standalone workflows, testing, and studies that
vary initial distributions, node counts, conditioning order, output cadence, or
phase-space evolution scheme without recompiling C++.

Applications that already manage simulation setup in C++ can configure QMOMKit
directly instead. Both approaches describe the same numerical setup and lead to
public case objects that can be passed to the same reconstruction routines:

```text
configuration file:
    .qmom file
        -> load_qmom(...) or load_cqmom(...)
        -> QmomCase or CqmomCase

direct API:
    C++ setup
        -> QmomCase, CqmomDenseCase, or CqmomLeanCase
```

Neither approach is preferred universally. Use the one that fits the surrounding
application.

For CQMOM, `load_cqmom(...)` returns a `CqmomCase`: the configured problem
definition used to compute quadrature and, when needed, create an initial
phase-space evolution state. Direct C++ setup that bypasses configuration text uses the
public explicit CQMOM case types, such as `CqmomDenseCase` for full tensor
moments or `CqmomLeanCase` for required moment data. These explicit cases are
ordinary public workflows for applications that intentionally manage moment
storage themselves.

Think of a configuration file as the part of a case a researcher would want to
record, rerun, and compare. It captures the moment-method setup. The
application still supplies the surrounding model: physical callbacks, meshes,
boundary conditions, units, timestep loop, and any interpretation of
phase-space names.

Use configuration files when:

- a workflow should be easy to rerun with different distributions or node
  counts;
- a driver should load setup data from disk;
- examples or validation cases need a compact, readable input file.

Use direct C++ construction when:

- the setup is generated programmatically;
- moments come directly from another solver or data structure;
- you need full application-level control over options and callbacks.

Configuration files and direct C++ setup have the same boundary: they describe
the QMOMKit case setup. Physics callbacks, model coefficients, mesh data,
boundary conditions, units, timestep loops, and detailed phase-space evolution control remain
application code. Not every public C++ option has a `.qmom` key: `.qmom`
exposes the supported case and runtime configuration surface, while advanced
low-level controls, such as reconstruction thresholds and source-integration
acceptance callbacks, remain API-only.

## Direct C++ Setup At A User Level

Direct C++ setup supplies the same concepts that appear in a `.qmom` file, but
through public types and options rather than through text sections:

| Concept | Configuration file | Direct C++ setup |
|---|---|---|
| phase space | `[phase_space]` names, node counts, conditioning order, adaptive branch settings | construct `QmomCase` or a public CQMOM case with the desired coordinate count, node counts, hierarchy, and ordering |
| initial moments | `[moments]` explicit values | provide `MomentSet`, CQMOM moment data, or a dense tensor representation from application data |
| analytical initialization | `[initial]` distribution, mixture, or independent product sections | call distribution moment-generation utilities, then pass the generated moments into the case setup |
| numerical settings | `[qmom]`, `[cqmom]`, and `[vandermonde]` sections | set the corresponding QMOM, CQMOM, adaptive, diagnostic, and solver options in C++ |
| phase-space evolution options | `[time_integration]`, `[time_integration.node_advection]`, and optional realizability settings | set `QmomPhaseSpaceEvolutionOptions`, `CqmomPhaseSpaceEvolutionOptions`, `NodeAdvectionOptions`, source-update options, and realizability settings in C++ |
| output choices | `[output]` directory, enabled outputs, and cadence | create and call `MomentWriter` or `QuadratureWriter` from the surrounding driver |

This mapping is meant to orient application integration, not replace the public
API reference or examples. Use direct setup when the surrounding simulation
already owns configuration data, creates moments at runtime, or needs to choose
options programmatically.

## Example: Direct C++ Setup

Direct C++ setup is useful when the surrounding simulation already owns
configuration data, moments are generated at runtime, parameters are chosen
programmatically, or QMOMKit is embedded in a larger solver.

The two setup paths differ only in how the case data enters QMOMKit:

```text
configuration file workflow:
    .qmom file
        -> load_qmom(...) or load_cqmom(...)
        -> QmomCase or CqmomCase
        -> compute_qmom(...) or compute_cqmom(...)

direct C++ workflow:
    application data
        -> create/configure a public case object
        -> compute_qmom(...) or compute_cqmom(...)
```

For one-dimensional QMOM, direct setup can create `QmomCase` directly from the
moment vector, node count, and options:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = false;

    const qmomkit::QmomCase qmom_case{
        qmomkit::MomentSet({1.0, 0.8, 2.8, 4.4}),
        2,
        qmom_options};

    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmom_case);
    if (!qmom.ok())
    {
        std::cerr << "QMOM failed: " << qmom.message << '\n';
        return 1;
    }

    return 0;
}
```

For multidimensional CQMOM, direct setup can use an explicit public
representation. This example uses a dense raw moment tensor because the
application already owns full tensor moments:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {0, 1});

    const qmomkit::RawMomentTensor moments{
        {4, 4},
        {
            1.0, 0.44, 0.424, 0.5552,
            0.55, 0.242, 0.2332, 0.30536,
            1.2125, 0.5335, 0.5141, 0.67318,
            1.121875, 0.493625, 0.475675, 0.622865,
        }};

    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = false;

    const qmomkit::CqmomDenseCase cqmom_case{
        hierarchy,
        moments,
        {},
        qmom_options};

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
    if (!cqmom.ok())
    {
        std::cerr << "CQMOM failed: " << cqmom.message << '\n';
        return 1;
    }

    return 0;
}
```

Use `CqmomDenseCase` when the application has full rectangular tensor moments.
Use `CqmomLeanCase` when the application intentionally provides only the
required CQMOM moments. If you want QMOMKit to derive required CQMOM moments
from a configured independent-product initialization, `load_cqmom(...)` remains
the most concise route to the representation-neutral `CqmomCase`.

Phase-space evolution is still optional. QMOM phase-space evolution advances a `QmomState` created from
`QmomCase`. CQMOM phase-space evolution advances a `CqmomState` created from `CqmomCase`:

```cpp
qmomkit::QmomState qmom_state = qmom_case.initial_state();
const qmomkit::QmomPhaseSpaceEvolutionResult qmom_next =
    qmomkit::advance_qmom_phase_space(qmom_state, qmom_case.node_count, time, time_step, velocity, source, options);

qmomkit::CqmomState state = cqmom_case.initial_state();
const qmomkit::CqmomPhaseSpaceEvolutionResult next =
    qmomkit::advance_cqmom_phase_space(state, time, time_step, velocity, source, options);
```

The `velocity`, `source`, and phase-space evolution options are application code. The
[Phase-space evolution](phase_space_evolution.md) guide explains those callbacks and options in detail.

## File Structure

QMOMKit case configuration files use the `.qmom` extension.

Every QMOMKit case configuration defines a phase space:

```text
[phase_space]
names = xi
node_counts = 2
```

The phase-space structure determines the configured formulation. A
one-dimensional phase space uses QMOM. A phase space with more than one
dimension uses CQMOM. Coordinate names are labels only; QMOMKit does not infer
physics or method choice from names such as `d`, `u`, or `T`.

`format_version = 1` is optional. If omitted, QMOMKit assumes the current
format. If an unsupported version is written explicitly, loading fails with a
configuration error.

The current `.qmom` case format selects only between QMOM and CQMOM from
phase-space dimensionality.

The main sections are:

- `[phase_space]`: phase-space names, node counts, conditioning order, and
  adaptive-branch settings.
- `[moments]`: explicit raw moments supplied by the user.
- `[initial]`: distribution-based moment initialization.
- `[qmom]`: one-dimensional inversion options shared by QMOM and CQMOM.
- `[cqmom]`: CQMOM-specific reconstruction options.
- `[time_integration]`: phase-space evolution operator-ordering choices.
- `[time_integration.node_advection]`: optional advection node-trajectory
  integration scheme and embedded-RK controls.
- `[time_integration.realizability]`: optional QMOM moment-source
  realizability settings.
- `[time_integration.realizability.<coordinate>]`: optional CQMOM
  coordinate-specific moment-source realizability settings.
- `[output]`: optional output directory, enabled outputs, and cadence.

Configuration files do not define physics callbacks, meshes, boundary
conditions, or solver loops. Those remain application code.

## Loading A Configuration

Applications can load a configuration through the public configuration API:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    const auto configured =
        qmomkit::config::load_qmom("configuration.qmom");

    if (!configured)
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomResult result = qmomkit::compute_qmom(configured.value());
    if (!result.ok())
    {
        std::cerr << "QMOM failed: " << result.message << '\n';
        return 1;
    }
}
```

Use `load_qmom(...)` for one-dimensional phase-space configurations. For
ordinary multidimensional CQMOM configurations, use `load_cqmom(...)`:

```cpp
const auto configured =
    qmomkit::config::load_cqmom("configuration.qmom");

if (!configured)
{
    // Report configured.error(), then stop.
}

const qmomkit::CqmomCase cqmom_case = configured.value();

const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);
if (!result.ok())
{
    // Report result.message, then stop.
}
```

`CqmomCase` stores the setup loaded from the file and does not change during
phase-space evolution. `compute_cqmom(...)` returns the quadrature result for that setup; it
does not return a `CqmomState`. `CqmomCase::initial_state()` creates the
evolving `CqmomState` used by optional CQMOM phase-space evolution. QMOM phase-space evolution is also
optional and advances an evolving `QmomState` created from `QmomCase`. Use
`load_cqmom_dense(...)` only when an application deliberately chooses the
explicit dense representation workflow. Optional output settings are loaded
separately with `load_output_options(...)`; phase-space evolution controls are loaded with
`load_qmom_phase_space_evolution_options(...)` or `load_cqmom_phase_space_evolution_options(...)`.

## QMOM Configuration

Use a one-dimensional phase space for QMOM: a workflow with one phase-space
coordinate and a raw moment sequence `m0, m1, ...`. In this case, `node_counts`
contains exactly one positive integer: the requested number of QMOM nodes.

### Explicit Moments

Use explicit moments when another code, measurement pipeline, or data file
already provides raw moments:

```text
[phase_space]
names = xi
node_counts = 2

[moments]
values = 1.0, 0.8, 2.8, 4.4
```

`values` are raw moments in ascending order:

```text
m0, m1, m2, ...
```

For `N` QMOM nodes, provide at least `2 * N` moments. `m0` is the total number
density; QMOMKit does not normalize it away.

### Distribution Initialization

Use `[initial] type = distribution` when the initial condition is easier to
describe as an analytical one-dimensional distribution than as a table of raw
moments:

```text
[phase_space]
names = xi
node_counts = 2

[initial]
type = distribution
family = lognormal
number_density = 2.5
maximum_order = 3
mu_log = -0.2
sigma_log = 0.4
```

Generated moments use:

```text
m_k = number_density * E[xi^k]
m_0 = number_density
```

`maximum_order` is optional and inclusive. If omitted, QMOMKit generates the
minimum required order from the configured quadrature, `2 * node_count - 1`. A
larger explicit value remains available when additional initialization moments
are desired; this applies to analytical and sampled initialization workflows.

Supported families and configuration parameter names are:

```text
normal:       mean, standard_deviation
lognormal:   mu_log, sigma_log
beta:        alpha, beta
gamma:       shape, scale
chi:         degrees_of_freedom
chi-squared: degrees_of_freedom
weibull:     shape, scale
```

These canonical family names are preferred. The parser also accepts
`log-normal` and `log_normal` for `lognormal`, and `chi_squared` and
`chisquared` for `chi-squared`.

All supplied distribution parameters must be finite. The remaining family
constraints are:

- `normal`: `standard_deviation > 0`; `mean` may be any finite value.
- `lognormal`: `sigma_log > 0`; `mu_log` may be any finite value.
- `beta`: `alpha > 0` and `beta > 0`; its support is `[0, 1]`.
- `gamma` and `weibull`: `shape > 0` and `scale > 0`.
- `chi` and `chi-squared`: `degrees_of_freedom > 0`.

`number_density` for a single distribution must be finite and non-negative.

### Mixture Initialization

Use `[initial] type = mixture` when one analytical distribution is not enough
to describe the initial population, such as a population with multiple modes
or regimes:

```text
[phase_space]
names = xi
node_counts = 3

[initial]
type = mixture
maximum_order = 5

[initial.mixture.component.0]
family = gamma
number_density = 1.5
shape = 2.0
scale = 0.4

[initial.mixture.component.1]
family = lognormal
number_density = 0.8
mu_log = -0.2
sigma_log = 0.5
```

Mixture component densities are absolute:

```text
m_k = sum_r number_density_r * E_r[xi^k]
m_0 = sum_r number_density_r
```

Component sections must be contiguous starting at
`[initial.mixture.component.0]`. Each component `number_density` must be
finite and non-negative, and at least one component must have positive density.

### QMOM Options

The optional `[qmom]` section controls the one-dimensional inversion options
available through configuration:

```text
[qmom]
adaptive = true
```

For ordinary fixed-node QMOM, omit the section or set `adaptive = false`.

## CQMOM Configuration

Use a multidimensional phase space for CQMOM: a workflow with several
phase-space coordinates and a coupled quadrature representation. A file-based
CQMOM reconstruction path is:

```text
configuration
    -> load_cqmom(...)
    -> CqmomCase
    -> compute_cqmom(...)
    -> CqmomResult / CqmomQuadrature
```

After reconstruction, users may inspect the quadrature, evaluate derived
quantities, write output, or optionally evolve the distribution. Optional CQMOM
phase-space evolution creates a `CqmomState` with `CqmomCase::initial_state()` and advances
that state with `advance_cqmom_phase_space(...)`.

An application that loads a normal `CqmomCase` and later supplies complete raw
tensor moments may create a replacement case with
`qmomkit::config::with_initial_moments(base_case, tensor)`. This immutable
operation preserves the loaded hierarchy, support, storage choice, and
reconstruction policies. It performs tensor structural validation, while CQMOM
realizability remains a reconstruction-time concern. It is therefore equivalent
to using explicit moments in an otherwise identical CQMOM configuration without
serializing or reparsing configuration text.

The phase-space section defines the dimension count, node hierarchy, and
conditioning order. These settings connect the configuration file to the
physical model: every callback coordinate vector, output column, and moment
multi-index follows the same phase-space dimension order.

```text
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 0, 1
adaptive_branch_counts = false
```

The key ideas are:

- `names` are labels for humans and output. Distribution families and supports
  are selected explicitly in their own sections.
- `node_counts[d]` is the requested node count in phase-space dimension `d`.
- `conditioning_order` controls the CQMOM conditioning order. It does not
  change tensor storage order or the meaning of phase-space dimensions.
- `adaptive_branch_counts` controls whether adaptive CQMOM may accept fewer
  nodes on individual branches.

If `conditioning_order` is omitted, QMOMKit uses identity order:

```text
0, 1, 2, ...
```

### Independent-Product Initialization

Use `[initial] type = independent_product` when a multidimensional initial
condition can be described as independent one-dimensional marginals multiplied
by one global number density. `initial.dimension.N` always refers to
phase-space dimension `N`, not to CQMOM conditioning level `N`; dimension
sections are contiguous from zero and follow coordinate and tensor-index order.
Each marginal may be one supported distribution family or a `type = mixture`
marginal. Different dimensions may use different families:

```text
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 0, 1
adaptive_branch_counts = false

[initial]
type = independent_product
number_density = 2.5

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
type = mixture

[initial.dimension.1.mixture.component.0]
family = lognormal
number_density = 1.5
mu_log = -0.3
sigma_log = 0.25

[initial.dimension.1.mixture.component.1]
family = lognormal
number_density = 3.25
mu_log = 0.4
sigma_log = 0.45
```

The nested component sections must be numbered contiguously from
`initial.dimension.N.mixture.component.0`. Their `number_density` values are
absolute contributions to that marginal mixture, not fractions. QMOMKit first
normalizes each completed marginal through its zeroth raw moment, then combines
the normalized marginals as an independent product and applies the global
`initial.number_density` exactly once. Per-dimension number densities are not
part of this model.

This configuration represents independent products of marginals, including
marginal mixtures. It does not represent correlated joint distributions,
copulas, multidimensional mixture components, or mixtures of independent-product
populations.

For the normal `load_cqmom(...)` path, `node_counts` and
`conditioning_order` determine which raw moments are required. Analytical
independent-product initialization generates only those required moments and
does not construct a dense `RawMomentTensor`.

`maximum_orders[d]` optionally gives the largest moment order to generate in
dimension `d` when building a full tensor:

```text
extent[d] = maximum_orders[d] + 1
```

`maximum_orders` and `node_counts` describe different parts of the workflow.
Moment orders control how many raw tensor entries are generated for full tensor
workflows. Node counts control the requested CQMOM quadrature size. New users
following `load_cqmom(...)` do not need `maximum_orders`.

### Optional CQMOM Coordinate Support

CQMOM coordinates are unconstrained by default, which means real-line
reconstruction is used for every coordinate. When a phase-space coordinate has a
known closed bounded support, declare it with paired phase-space-order
`support_lower` and `support_upper` lists:

```text
[phase_space]
names = d, chi, u
node_counts = 2, 3, 2
conditioning_order = 0, 1, 2
support_lower = none, 0.0, none
support_upper = none, 1.0, none
```

Use `none` for unconstrained coordinates. Numeric lower and upper entries must be
finite values configured together, one entry per phase-space dimension, with
`lower < upper`. Supports are generic numerical bounds. They are not inferred
from coordinate names and do not assign physical units or particle semantics.

A shifted interval uses the same syntax:

```text
[phase_space]
names = x, y
node_counts = 2, 2
conditioning_order = 1, 0
support_lower = none, 2.0
support_upper = none, 5.0
```

Support is indexed by phase-space dimension, not by conditioning level. In the
example above, the second coordinate is constrained to `[2, 5]` wherever it is
reconstructed, even though it is the root coordinate in the conditioning order.
The same rule applies when a bounded coordinate is reconstructed at an
intermediate or final conditional level.

For a bounded coordinate, QMOMKit accepts a reconstructed CQMOM candidate only
when its returned abscissas lie inside the configured interval. The bounds are
carried by `CqmomCase` and `CqmomState`, so they apply to initial configuration
reconstruction and later phase-space evolution stage reconstructions. In
adaptive branch-count mode, a candidate with out-of-interval nodes retries
lower node counts down to one. In fixed-count mode, such a candidate fails with
diagnostics. Returned nodes are not clipped or projected. These bounds constrain
quadrature abscissas only; they do not configure moment realizability.

### Advanced: Explicit Dense Representation Workflow

Use explicit tensor moments when your application already has a full raw
moment tensor from external data, validation, debugging, or interoperability:

```text
[phase_space]
node_counts = 1, 1
conditioning_order = 0, 1

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.5
M[1,0] = 2.0
M[1,1] = 1.0
```

Tensor entries use phase-space dimension order. For a two-dimensional tensor,
`M[1,0]` means order `1` in dimension `0` and order `0` in dimension `1`.
Applications that need this representation can load it with
`load_cqmom_dense(...)` and use the explicit dense representation APIs.
Independent-product configurations loaded through this dense representation
workflow may omit `initial.maximum_orders`; QMOMKit then derives the minimum
dense tensor extents from the configured quadrature. Explicit larger extents
remain available for analytical and sampled initialization. Dense remains a
first-class supported expert workflow for callers that intentionally manage
full tensor moment storage.

### Advanced: CQMOM Diagnostics And Solver Options

Most applications can omit the `[cqmom]` and `[vandermonde]` sections. Add
them when you need more or less reconstruction diagnostics, or when a
validation study deliberately compares solver policies:

```text
[cqmom]
moment_storage = lean
diagnostic_level = detailed|basic|none

[vandermonde]
solver = automatic
```

`moment_storage` selects the internal raw-moment representation for
representation-neutral `load_cqmom(...)` workflows. Supported values are `lean`
and `dense`. If omitted, `lean` is used explicitly. The selected moment storage
is reported by `CqmomCase::moment_storage()` and `CqmomState::moment_storage()`.
This setting is distinct from `vandermonde.solver = dense`, which selects a
linear solver for conditional reconstruction rather than moment storage.

`diagnostic_level` controls reconstruction diagnostic retention and reporting
only. It does not change solver selection, residual acceptance, fallback
behavior, inversion results, phase-space evolution results, or success/failure status.
Supported values are `detailed`, `basic`, and `none`. If omitted, `detailed` is
used.

`vandermonde.solver` controls the branch-local conditional Vandermonde solve.
Supported values are:

| Value | Behavior |
|---|---|
| `automatic` | Default; currently resolves to Björck-Pereyra with eligible dense fallback. |
| `dense` | Dense Gaussian-elimination solve only. |
| `pivoted_qr` | Pivoted-QR solve only, with numerical-rank detection and no fallback. |
| `svd` | SVD solve only, with rank and convergence checks and no fallback. |
| `bjorck_pereyra` | Structured Björck-Pereyra solve only, with no fallback. |
| `bjorck_pereyra_with_dense_fallback` | Björck-Pereyra followed by eligible dense fallback. |

If omitted, `automatic` is used. The `pivoted_qr` and `svd` policies solve the
weighted conditional Vandermonde system directly and are selected only when
requested explicitly:

```text
[vandermonde]
solver = svd
```

Neither explicit direct policy falls back to the dense solver if it fails. A
rank-deficient conditional Vandermonde matrix makes the solve fail, and SVD
Jacobi nonconvergence is also reported as a failure. QMOMKit does not compute a
truncated, regularized, minimum-norm, or pseudoinverse solution in these paths.
SVD rank and Jacobi tolerances are internal defaults and are not configurable
through `.qmom` files. Loaded CQMOM cases and states expose the configured
policy through `vandermonde_solver()`.

### Fixed And Adaptive CQMOM

The default CQMOM setup is fixed, which gives a rectangular quadrature layout
and is usually the clearest starting point:

```text
adaptive_branch_counts = false

[qmom]
adaptive = false
```

Adaptive branch quadrature is requested explicitly when a workflow can accept
branch-local node reduction:

```text
[phase_space]
node_counts = 2, 2
conditioning_order = 0, 1
adaptive_branch_counts = true
```

Enabling adaptive CQMOM branch counts automatically enables adaptive Wheeler
for all marginal and conditional inversions. In a CQMOM configuration,
`[qmom] adaptive = true` therefore requires
`phase_space.adaptive_branch_counts = true`; it is rejected for fixed-count
CQMOM. The `[qmom] adaptive` setting remains valid for standalone
one-dimensional QMOM.

## Phase-Space Evolution Options

A `.qmom` file may contain both case settings and phase-space evolution settings.
`load_qmom(...)` and `load_cqmom(...)` return only the case configuration, but
they validate all recognized settings in the file. Phase-space evolution settings must
still be loaded into phase-space evolution options explicitly.

| Settings | Validated by case loading | Stored in case/state | Must be loaded separately |
| --- | --- | --- | --- |
| QMOM/CQMOM inversion settings | yes | yes | no |
| Vandermonde and CQMOM storage settings | yes | yes | no |
| `time_integration.scheme` | yes | no | yes |
| `time_integration.node_advection.*` | yes | no | yes |
| `time_integration.realizability.*` | yes | no | yes |

A configuration file may select the phase-space evolution operator-ordering scheme:

```text
[time_integration]
scheme = strang
```

Supported values are:

```text
strang
advsrc
srcadv
```

These choose how a phase-space evolution driver composes phase-space advection and
moment-source evolution. If `scheme` is omitted, the default is
`strang`, the usual second-order composition for cases where both parts are
active.

A configuration file may also select the advection node-advection scheme:

```text
[time_integration.node_advection]
scheme = rk23
```

Supported values are:

```text
rk23
rk45
ssp_rk3
```

`rk23` is the default. It is the adaptive Bogacki-Shampine 3(2) pair.
`rk45` is the adaptive Dormand-Prince RK5(4)7M pair from Table 2 of Dormand
and Prince (1980). `ssp_rk3` is a distinct fixed-stage method that advances
exactly over the timestep supplied by the application. The choice affects only
advection quadrature-node trajectories; it does not change moment-source integration,
realizability policy, or phase-space evolution splitting.

The adaptive controls apply to `rk23` and `rk45` and live under a shared nested
section. `initial_step` defaults to `1.0e-6 s` and may be overridden:

```text
[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-7
absolute_tolerance = 1.0e-10
initial_step = 1.0e-4
max_step = 1.0e-2
min_step = 1.0e-12
max_steps = 200000
```

Supported adaptive embedded-RK keys are:

```text
relative_tolerance
absolute_tolerance
initial_step
max_step
min_step
max_steps
```

If omitted, the tolerances retain the C++ `AdaptiveRkOptions` defaults:
`relative_tolerance = 1.0e-6`, `absolute_tolerance = 1.0e-9`,
`max_step = full integration interval`,
`min_step = machine-dependent automatic minimum`, and
`max_steps = 1000000000`, and `initial_step = 1.0e-6 s`. There is no automatic
initial-step estimator: QMOMKit does not derive `initial_step` from the state,
derivatives, or tolerances. Explicit tolerances and step sizes must be positive finite
numbers, and `min_step` must not exceed `max_step` when both are provided.
`max_steps` must be a positive integer. The higher-order member advances an
accepted state, and both embedded methods use FSAL derivative reuse internally.

Drivers can print `format_node_advection_options(...)` at startup to report
the effective scheme and controls. For fixed-step node advection, that summary
reports the selected scheme only; it does not imply that adaptive controls are
active.

### CQMOM Moment-Source Active Counts

When adaptive branch counts are enabled, QMOMKit automatically reduces
branch-local node counts when moment information collapses. Phase-space
advection cannot increase active counts. After an accepted
moment-source update, each existing branch may increase by at most one
node toward its configured maximum when the updated moments support that
additional node.

Initial reconstruction starts each branch from the configured maximum and lets
adaptive Wheeler determine the active count. Newly created descendant branches
also begin from the configured maximum and adapt downward. Users configure the
maximum node counts, whether CQMOM branch counts are adaptive, and adaptive
Wheeler tolerances; moment-source node-count evolution is automatic.

The configuration file stores runtime numerical choices and the initial
numerical setup. Your application still provides the advection velocity and
source callbacks, model coefficients, and timestep plan. Load QMOM phase-space evolution
controls with `load_qmom_phase_space_evolution_options(...)` and CQMOM phase-space evolution controls
with `load_cqmom_phase_space_evolution_options(...)`; these return complete phase-space evolution
option objects with default empty runtime callbacks.

## Output Options

Examples can use `[output]` to decide where and when to write moment or
quadrature histories:

```text
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
```

If `[output]` is absent, QMOMKit performs no file output by default. Output
settings are application-controlled; numerical algorithms do not write files
unless your driver wires these options into writers.

`qmomkit::MomentWriter` and `qmomkit::QuadratureWriter` write text tables with
a commented header line and whitespace-separated columns. QMOMKit-native
example output uses `.txt` filenames. Use CSV only for custom application data
that is not written through the QMOMKit output helpers.

CQMOM quadrature output uses one row per active terminal node:

```text
time node_index path_0 ... weight x0 ...
```

Multiple rows may share one time value, and adaptive CQMOM may write a different
number of rows at different output times. Group rows by `time` to recover each
saved terminal quadrature.

QMOM quadrature output also uses one row per active node:

```text
time node_index weight x0
```

Create `QuadratureWriter::for_qmom(path, maximum_node_count)` with the maximum
active node count the file should accept. Fixed and adaptive QMOM use the same
writer; inactive adaptive nodes are not padded or written as zero-weight nodes.

Moment output records native integer-order raw moments. Generalized moments
with real exponents are evaluated by caller code from a valid quadrature result
and are not part of the standard writer output.

## Runtime Realizability Settings

For QMOM, the realizability settings in
`[time_integration.realizability]` enable checks during moment-source SSP-RK3
evolution:

```text
[time_integration.realizability]
check_mode = final_only
support = nonnegative
accept_boundary = true
tolerance = 1e-12
relative_tolerance = true
```

Supported check modes are:

- `none`: do not perform a runtime realizability check;
- `final_only`: check the transported moments after the complete
  moment-source integration step;
- `every_stage`: check the transported moments after each SSP-RK3
  moment-source integration stage.

Supported supports are:

```text
real_line
nonnegative
bounded_interval
```

For bounded intervals, provide explicit bounds:

```text
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = 0.0
upper = 1.0
```

The QMOM settings apply only to moment-source evolution. They do not apply to
node advection, quadrature evolution, or reconstructed nodes and weights.

For CQMOM, use coordinate-specific realizability settings for the individual
phase-space coordinates that require checks. A zero-based coordinate suffix
identifies the coordinate and must be smaller than the configured phase-space
rank:

```text
[time_integration.realizability.0]
check_mode = final_only
support = bounded_interval
lower = 0.0
upper = 1.0

[time_integration.realizability.3]
check_mode = every_stage
support = real_line
```

Support is not inferred from a distribution family or a phase-space name.
For CQMOM, each configured coordinate setting applies only to the corresponding
one-dimensional marginal moment sequence:

\[
M[0,\ldots,k_d,\ldots,0].
\]

Coordinates without a setting are not checked. CQMOM realizability settings do
not check mixed moments or conditional moments, and they do not validate
quadrature reconstruction, adaptive topology, or branch counts. These are
intentionally separate from the coordinate-marginal realizability checks.

## Practical Guidance

Keep units consistent. Distribution parameters, moment values, and phase-space
coordinate names should refer to the same coordinate convention.

Keep moment normalization explicit. QMOMKit uses unnormalized raw moments of the
number density function. `m0` is the total number density. Dividing by `m0`
changes the setup to a unit-number-density problem.

Generate enough moments. For `N` QMOM nodes, an analytical initialization needs
`maximum_order >= 2 * N - 1`. For normal CQMOM configuration, `node_counts` and
`conditioning_order` define the reconstruction setup. Explicit full-tensor
workflows must choose tensor extents and `maximum_orders` large enough for the
requested CQMOM node counts.

Separate distribution settings from quadrature settings. Initial distributions
describe raw moments. `node_counts`, `conditioning_order`, and adaptive options
describe quadrature reconstruction.

Check configuration errors before computing quadrature. The public loading API
returns a diagnostic string through `error()` with source and line information
when available.

Check numerical statuses after loading. A valid configuration can still lead to
a failed QMOM or CQMOM reconstruction if the requested quadrature is not
accepted.

Remember precedence rules. Explicit `[moments]` values are the source of truth
when present. If explicit moments and `[initial]` are both present, the explicit
moments are used.

## Example Files

Small QMOMKit configuration examples live under `examples/configuration/`:

- `qmom_explicit_moments.qmom`
- `qmom_mixture_absolute_density.qmom`
- `realizability_policy_nonnegative.qmom`

Phase-space evolution examples with local configuration files include:

- `examples/qmom_1d_breakage_advection/configuration.qmom`
- `examples/cqmom_phase_space_evolution/configuration.qmom`

Monte Carlo companion examples live under `tools/monte_carlo/configuration/`:

- `qmom_distribution_lognormal.qmom`
- `cqmom_independent_product_2d.qmom`

Build and run examples as described in `installation.md`, then edit the
configuration files to explore different initial moments, node counts, and
output settings.
