# Configuring QMOMKit workflows

A `.qmom` file describes the repeatable numerical part of a case: the phase
space, the initial condition, and the numerical settings. Keeping that in a text
file means a case can be rerun, varied, and compared without recompiling, which
is why the examples and validation studies use one.

Configuration files are optional. Everything they express can also be built
directly in C++, which is the better choice when the surrounding application
already owns the setup data — see [Configuring without a file](#configuring-without-a-file).

## Your first configuration file

The smallest useful file declares a phase space and an initial condition. Save
this as `configuration.qmom`:

```ini
[phase_space]
names = xi
node_counts = 2

[moments]
values = 1.0, 0.8, 2.8, 4.4
```

That is a one-dimensional phase space with two quadrature nodes, initialised
from four raw moments. Load it and reconstruct:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    const auto configured = qmomkit::config::load_qmom("configuration.qmom");
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

    for (std::size_t node = 0; node < result.weights.size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << result.weights[node]
                  << ", abscissa=" << result.abscissas[node] << '\n';
    }
}
```

Two checks matter and they are separate. `configured` reports whether the file
was read and understood; `result.ok()` reports whether the requested quadrature
could actually be reconstructed from those moments. A perfectly valid file can
still describe moments that no two-node quadrature reproduces.

The rest of this guide builds outward from that file: a larger phase space, an
initial condition described analytically rather than by moments, and the
settings that control reconstruction, time evolution, and output.

## How a configuration file is organised

A `.qmom` file is a sequence of named sections. Only `[phase_space]` and an
initial condition are required; everything else has a default.

| Section | Purpose | Covered in |
| --- | --- | --- |
| `[phase_space]` | Coordinates, node counts, conditioning order, supports, scaling | [Describing the phase space](#describing-the-phase-space) |
| `[moments]` | Explicit raw moments | [Explicit moments](#explicit-moments) |
| `[initial]` | Moments generated from analytical distributions | [Setting the initial condition](#setting-the-initial-condition) |
| `[qmom]` | One-dimensional inversion options, used by QMOM and CQMOM alike | [Reconstruction settings](#reconstruction-settings) |
| `[cqmom]` | CQMOM diagnostics | [Advanced settings](#cqmom-diagnostics-and-solver-options) |
| `[vandermonde]` | Conditional solver policy | [Advanced settings](#cqmom-diagnostics-and-solver-options) |
| `[time_integration]` | Splitting scheme, node advection, realizability | [Time-evolution settings](#time-evolution-settings) |
| `[output]` | Output directory, selection, and cadence | [Output](#output) |

Sections may appear in any order. Unknown keys and sections are reported rather
than ignored, with the file and line number, so a typo surfaces as a load error
instead of a silently missing setting:

```text
configuration.qmom:4: Unknown key 'phase_space.node_countz'.
```

A configuration file stops at the numerical setup. Physics callbacks, model
coefficients, meshes, boundary conditions, units, and the timestep loop are
application code — see [Callbacks](callbacks.md).

## Describing the phase space

The `[phase_space]` section is the backbone of a case. Its dimension count
decides the method: one coordinate uses QMOM, several use CQMOM. Coordinate
names are labels for humans and output columns; QMOMKit infers no physics from
a name such as `d`, `u`, or `T`.

```ini
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 0, 1
adaptive_branch_counts = false
```

The keys are:

- `names` label the coordinates.
- `node_counts[d]` is the requested number of nodes in phase-space dimension
  `d`. For one-dimensional QMOM this is a single positive integer.
- `conditioning_order` sets the order in which CQMOM conditions the dimensions.
  It does not reorder coordinates, moments, or output; if omitted, dimensions
  are conditioned in their natural order `0, 1, 2, ...`.
- `adaptive_branch_counts` allows individual CQMOM branches to use fewer nodes
  than requested when their local moments cannot support more.

This section is what ties the file to your model: every callback coordinate
vector, output column, and moment multi-index follows the phase-space dimension
order declared here.

### Coordinate supports

Coordinates are unconstrained by default, so abscissas may be placed anywhere on
the real line. When a coordinate has a known closed bounded support, declare it
with paired `support_lower` and `support_upper` lists in phase-space dimension
order, using `none` where a coordinate is unconstrained:

```ini
[phase_space]
names = d, chi, u
node_counts = 2, 3, 2
conditioning_order = 0, 1, 2
support_lower = none, 0.0, none
support_upper = none, 1.0, none
```

A one-dimensional QMOM case uses the same keys with a single entry each:

```ini
[phase_space]
node_counts = 2
support_lower = 0.0
support_upper = 1.0
```

Numeric lower and upper entries must be finite, configured together, one per
phase-space dimension, with `lower < upper`. These are generic numerical bounds:
they are not inferred from coordinate names and carry no units or physical
meaning.

Support is indexed by phase-space dimension, not by conditioning level:

```ini
[phase_space]
names = x, y
node_counts = 2, 2
conditioning_order = 1, 0
support_lower = none, 2.0
support_upper = none, 5.0
```

Here the second coordinate is constrained to `[2, 5]` wherever it is
reconstructed, even though it is the root coordinate in the conditioning order.

QMOMKit accepts a reconstructed candidate only when its abscissas lie inside the
configured interval; nodes are never clipped or projected onto it. The bounds
travel with the case and its evolving state, so they apply to the initial
reconstruction and every later one. In adaptive mode a candidate with
out-of-interval nodes retries at lower node counts down to one; in fixed-count
mode it fails with diagnostics. These bounds constrain quadrature abscissas
only — they are not moment realizability settings, which are covered under
[runtime realizability](#runtime-realizability-settings).

### Automatic reconstruction scaling

Set `automatic_scaling = true` in `[phase_space]` to normalize coordinates
locally for a QMOM or CQMOM reconstruction:

```ini
[phase_space]
node_counts = 2
automatic_scaling = true
```

The default is `false`. Scaling is reconstruction-local: raw moments, source
callbacks, phase-space states, configured bounds, returned nodes, and bounded
support diagnostics remain in physical coordinates. QMOMKit selects one
positive scale per physical coordinate from the current raw moments, preferring
the RMS scale `sqrt(M[2e_d] / M[0])`, then the magnitude of the mean, and
finally one. For CQMOM, scales are selected once from the complete required
moment set before recursive conditioning.

Configured bounds are mapped to the normalized reconstruction only internally.
The bounded-reconstruction tolerance therefore applies in normalized
coordinates. If a safe transformation cannot be formed, QMOMKit follows the
ordinary unscaled reconstruction path.
## Setting the initial condition

A case needs initial raw moments. There are three ways to obtain them, and the
choice is usually made for you by where the data comes from:

- **`[moments]`** when another code, a measurement pipeline, or a data file
  already gives you raw moments;
- **`[initial]`** when the starting population is easier to describe as an
  analytical distribution and you want QMOMKit to generate the moments;
- **sampled moments** when you want the initial moments to come from a finite
  sampled population rather than from the analytical distribution itself.

If both `[moments]` and `[initial]` are present, the explicit `[moments]` values
win.

Sampling is the one option that is not a configuration-file setting. The
optional [Monte Carlo companion](monte_carlo.md) draws an ensemble from a
distribution and estimates its raw moments, which you then hand to QMOMKit from
C++ the same way you would any other explicit moments. Reach for it when a
deterministic run and a stochastic reference run should start from exactly the
same realization, so that any later difference between them comes from the
evolution rather than from a different starting population.

See [Sample and initialize from one realization](monte_carlo.md#sample-and-initialize-from-one-realization)
in the Monte Carlo guide for the API, and the `sampling.cpp` file in examples
such as [`cqmom_5d_particles_breakage_heat-transfer`](../../examples/cqmom_5d_particles_breakage_heat-transfer/README.md)
for a complete case built this way.

### Explicit moments

Use explicit moments when another code, measurement pipeline, or data file
already provides raw moments:

```ini
[phase_space]
names = xi
node_counts = 2

[moments]
values = 1.0, 0.8, 2.8, 4.4
```

`values` are raw moments in ascending order:

```ini
m0, m1, m2, ...
```

For `N` QMOM nodes, provide at least `2 * N` moments. `m0` is the total number
density; QMOMKit does not normalize it away.

### A single distribution

Use `[initial] type = distribution` when the initial condition is easier to
describe as an analytical one-dimensional distribution than as a table of raw
moments:

```ini
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

### A mixture of distributions

Use `[initial] type = mixture` when one analytical distribution is not enough
to describe the initial population, such as a population with multiple modes
or regimes:

```ini
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

### Independent marginals for CQMOM

Use `[initial] type = independent_product` when a multidimensional initial
condition can be described as independent one-dimensional marginals multiplied
by one global number density. `initial.dimension.N` always refers to
phase-space dimension `N`, not to CQMOM conditioning level `N`; dimension
sections are contiguous from zero and follow coordinate and tensor-index order.
Each marginal may be one supported distribution family or a `type = mixture`
marginal. Different dimensions may use different families:

```ini
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
copulas, multidimensional mixture components, or mixtures of
independent-product populations.

For the normal `load_cqmom_prepared(...)` path, `node_counts` and
`conditioning_order` determine which raw moments are required. Analytical
independent-product initialization generates only those required moments.

To retain more initial moments than the reconstruction needs, see
[retaining extra moments](#retaining-extra-moments).

## Reconstruction settings

The default setup is fixed: every branch uses exactly the requested node counts
and the quadrature layout is rectangular. That is usually the clearest starting
point, and it needs no configuration at all.

The default CQMOM setup is fixed, which gives a rectangular quadrature layout
and is usually the clearest starting point:

```ini
adaptive_branch_counts = false

[qmom]
adaptive = false
```

Adaptive branch quadrature is requested explicitly when a workflow can accept
branch-local node reduction:

```ini
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

The optional `[qmom]` section carries the one-dimensional inversion settings
shared by QMOM and CQMOM:

```ini
[qmom]
adaptive = true
```

For ordinary fixed-node QMOM, omit the section or set `adaptive = false`.
Diagnostic level and solver policy are covered under
[advanced settings](#advanced-settings).
## Time-evolution settings

A `.qmom` file may contain both case settings and phase-space evolution
settings. `load_qmom(...)` and `load_cqmom_prepared(...)` return only the case
configuration, but they validate all recognized settings in the file.
Phase-space evolution settings must still be loaded into phase-space evolution
options explicitly.

| Settings | Validated by case loading | Stored in case/state | Must be loaded separately |
| --- | --- | --- | --- |
| QMOM/CQMOM inversion settings | yes | yes | no |
| Vandermonde and CQMOM diagnostic settings | yes | yes | no |
| `time_integration.scheme` | yes | no | yes |
| `time_integration.node_advection.*` | yes | no | yes |
| `time_integration.realizability.*` | yes | no | yes |

A configuration file may select the phase-space evolution operator-splitting
scheme:

```ini
[time_integration]
scheme = strang
```

Supported values are:

```ini
strang
advsrc
srcadv
```

These choose how a phase-space evolution driver composes phase-space advection and
moment-source evolution. If `scheme` is omitted, the default is
`strang`, the usual second-order composition for cases where both parts are
active.

A configuration file may also select the node-advection scheme:

```ini
[time_integration.node_advection]
scheme = rk23
```

Supported values are:

```ini
rk23
rk45
ssp_rk3
```

`rk23` is embedded third-order RK scheme of [Bogacki & Shampine (1989)](https://www.doi.org/10.1016/0893-9659(89)90079-7) and chosen by default.
The embedded fifth-order RK scheme of [Dormand & Prince (1980)](https://www.doi.org/10.1016/0771-050X(80)90013-3), `rk45`, and the third-order SPP-RK 3 scheme of [Shu & Osher (1988)](https://www.doi.org/10.1016/0021-9991(88)90177-5), `ssp_rk3`, are available as alternatives.
The choice affects only advection quadrature-node trajectories; it does not change moment-source integration,
realizability policy, or phase-space evolution splitting.

The adaptive controls apply to `rk23` and `rk45` and live under a shared nested
section. `initial_step` defaults to `1.0e-6 s` and may be overridden:

```ini
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

```ini
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
numbers, and an explicit `min_step` must not exceed the effective `max_step` —
whether `max_step` is explicitly configured or automatically derived from the
integration interval. `min_step` takes precedence over the tolerances: when the
error controller would ask for a smaller step, a finite step at the floor is
accepted instead, so a successful run does not by itself certify that every step
met the requested tolerance. `max_steps` must be a positive integer. The higher-order member advances an
accepted state, and both embedded methods use FSAL derivative reuse internally.

Drivers can print `format_node_advection_options(...)` at startup to report
the effective scheme and controls. For fixed-step node advection, that summary
reports the selected scheme only; it does not imply that adaptive controls are
active.

The configuration file stores runtime numerical choices and the initial
numerical setup. Your application still provides the advection velocity and
source callbacks, model coefficients, and timestep plan. Load QMOM phase-space
evolution controls with `load_qmom_phase_space_evolution_options(...)` and
CQMOM phase-space evolution controls with
`load_cqmom_phase_space_evolution_options(...)`; these return complete
phase-space evolution option objects with default empty runtime callbacks.

## Output

Examples can use `[output]` to decide where and when to write moment or
quadrature histories:

```ini
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

```ini
time node_index path_0 ... weight x0 ...
```

Multiple rows may share one time value, and adaptive CQMOM may write a
different number of rows at different output times. Group rows by `time` to
recover each saved terminal quadrature.

QMOM quadrature output also uses one row per active node:

```ini
time node_index weight x0
```

Create `QuadratureWriter::for_qmom(path, maximum_node_count)` with the maximum
active node count the file should accept. Fixed and adaptive QMOM use the same
writer; inactive adaptive nodes are not padded or written as zero-weight nodes.

Moment output records native integer-order raw moments. The same
`output.generalized_moments` option can append quadrature-derived generalized
moments for QMOM or CQMOM when a driver passes the parsed `OutputOptions` to a
state-based writer. It is an ordered
pipe-separated list of comma-separated exponent vectors: QMOM uses rank-one
requests such as `0.5 | 1.5 | 2.5`; a two-dimensional CQMOM case can use
`0.5, 0 | 1, 1.5`. Exponents must be finite and non-negative, and vectors must
match the physical phase-space dimension.

Normal state-based output passes the parsed `OutputOptions`, current state, and
current quadrature to `MomentWriter`. Ordinary columns still come directly from
the state. Generalized columns are appended last in configured order, use
physical coordinate order, and do not change moment storage, reconstruction, or
evolution. An integer-valued generalized request remains an explicit appended
column even when the corresponding ordinary moment is already written.

## Runtime realizability settings

For QMOM, the realizability settings in
`[time_integration.realizability]` enable checks during moment-source SSP-RK3
evolution:

```ini
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

```ini
real_line
nonnegative
bounded_interval
```

For bounded intervals, provide explicit bounds:

```ini
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = 0.0
upper = 1.0
```

For CQMOM, use coordinate-specific realizability settings for the individual
phase-space coordinates that require checks. A zero-based coordinate suffix
identifies the coordinate and must be smaller than the configured phase-space
rank:

```ini
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
one-dimensional marginal moment sequence: $M[0,\ldots,k_d,\ldots,0]$.


Coordinates without a setting are not checked. CQMOM realizability settings do
not check mixed moments or conditional moments, and they do not validate
quadrature reconstruction, adaptive topology, or branch counts. These are
intentionally separate from the coordinate-marginal realizability checks.

## Loading a configuration from C++

Use `load_qmom(...)` for a one-dimensional configuration and
`load_cqmom_prepared(...)` for a multidimensional one:

```cpp
const auto configured = qmomkit::config::load_cqmom_prepared("configuration.qmom");
if (!configured)
{
    std::cerr << configured.error() << '\n';
    return 1;
}

const qmomkit::CqmomPreparedCase cqmom_prepared_case = configured.value();

const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_prepared_case);
if (!result.ok())
{
    std::cerr << "CQMOM failed: " << result.message << '\n';
    return 1;
}
```

On failure, `error()` gives a diagnostic string with source and line
information where available.

A loaded case is the fixed setup; it does not change as a simulation runs.
`compute_qmom(...)` and `compute_cqmom(...)` reconstruct a quadrature from it
directly. To evolve it in time instead, create the runtime state with
`QmomCase::initial_state()` or `CqmomPreparedCase::initial_state()` and advance
that `QmomState` or `CqmomState` with `advance_qmom_phase_space(...)` or
`advance_cqmom_phase_space(...)`. See [Phase-space
evolution](phase_space_evolution.md) for those calls and the callbacks they
need.

**Case settings and evolution settings load separately.** `load_qmom(...)` and
`load_cqmom_prepared(...)` validate every recognised key in the file, but they return
only the case. Time-evolution and output settings must be requested explicitly:

| To load | Call |
| --- | --- |
| the case | `load_qmom(...)` or `load_cqmom_prepared(...)` |
| evolution controls | `load_qmom_phase_space_evolution_options(...)` or `load_cqmom_phase_space_evolution_options(...)` |
| output settings | `load_output_options(...)` |

The evolution loaders return complete option objects with empty callbacks, which
your application then fills in with its physics.

`load_cqmom(...)` reads the same file into an editable `CqmomCase` definition
instead, for an application that changes settings in code before preparing the
case with `qmomkit::CqmomPreparedCase{cqmom_case}`.

An application that loads a case and later obtains a complete raw moment tensor
can build a replacement case with
`qmomkit::config::with_initial_moments(base_cqmom_prepared_case, tensor)`. This
keeps the loaded hierarchy, supports, and reconstruction settings, and
validates the tensor's structure — equivalent to having written those moments
into the file, without serialising and reparsing it.

If the moments already arrive as the required CQMOM moment set rather than a
full tensor, pass `CqmomMomentData` instead:
`qmomkit::config::with_initial_moments(base_cqmom_prepared_case, moment_data)`.
The result is the same case, and it avoids building a full tensor whose extra
entries the case never reads.

Every case derived from one prepared base case shares its recursive branch
structure, so derive the cases for many cells or samples from one base case
rather than loading or preparing each one. Both overloads also accept a
`CqmomCase` as the base, with the same checks and errors, and then return a
`CqmomCase`.

Use the `RawMomentTensor` overload when the moments come from outside QMOMKit,
as for instance the Monte Carlo companion's sampled moments do. It selects the moments the
hierarchy needs, and reports a tensor of the wrong rank or size, or one missing
a required moment, as a configuration error you can inspect. Extracting those
moments yourself with `extract_cqmom_moments(...)` signals the same problems by
throwing instead.
## Configuring without a file

Applications that already own their setup data can skip the file entirely and
construct cases directly. Both routes describe the same thing and produce case
objects the same reconstruction routines accept, so the choice is about where
the setup naturally lives rather than about capability.

Prefer a configuration file when a case should be easy to rerun with different
distributions or node counts, when a driver reads its setup from disk, or when
an example needs a compact readable input. Prefer direct construction when the
setup is generated programmatically, when moments arrive from another solver or
data structure, or when you want full control over options that the file format
does not expose.

The same concepts appear either way:

| Concept | `.qmom` file | Direct C++ setup |
|---|---|---|
| Concept | `.qmom` file | Direct C++ setup |
|---|---|---|
| phase space | `[phase_space]` names, node counts, conditioning order, adaptive branch settings | construct `QmomCase` or a CQMOM case with the desired coordinate count, node counts, hierarchy, and ordering |
| initial moments | `[moments]` explicit values | provide `MomentSet`, CQMOM moment data, or a full tensor of raw moments from application data |
| analytical initialization | `[initial]` distribution, mixture, or independent product sections | call distribution moment-generation utilities, then pass the generated moments into the case setup |
| numerical settings | `[qmom]`, `[cqmom]`, and `[vandermonde]` sections | set the corresponding QMOM, CQMOM, adaptive, diagnostic, and solver options in C++ |
| phase-space evolution options | `[time_integration]`, `[time_integration.node_advection]`, and optional realizability settings | set `QmomPhaseSpaceEvolutionOptions`, `CqmomPhaseSpaceEvolutionOptions`, `NodeAdvectionOptions`, source-update options, and realizability settings in C++ |
| output choices | `[output]` directory, enabled outputs, and cadence | create and call `MomentWriter` or `QuadratureWriter` from the surrounding driver |

Not every C++ option has a `.qmom` key. The file format covers the case and
runtime configuration surface; low-level controls such as reconstruction
thresholds and source-integration acceptance callbacks stay API-only.

For one-dimensional QMOM, construct the case from the moments, node count, and
options:

```cpp
#include <qmomkit/qmomkit.hpp>

#include <iostream>

int main()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = false;

    const qmomkit::QmomCase qmom_case{qmomkit::MomentSet({1.0, 0.8, 2.8, 4.4}), 2, qmom_options};

    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmom_case);
    if (!qmom.ok())
    {
        std::cerr << "QMOM failed: " << qmom.message << '\n';
        return 1;
    }

    return 0;
}
```

For multidimensional CQMOM, build a `CqmomCase` from the moments the
hierarchy requires, then construct a `CqmomPreparedCase` from it for
reconstruction and time evolution, as below. If you would rather have QMOMKit
derive the required moments from a configured independent-product
initialisation, `load_cqmom_prepared(...)` remains the simpler route.

This example starts from a full tensor the application already owns and selects
the required moments from it:

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

    const qmomkit::CqmomCase cqmom_case{
        qmomkit::extract_cqmom_moments(moments, qmomkit::create_cqmom_moment_layout(hierarchy)),
        {},
        qmom_options};

    const qmomkit::CqmomPreparedCase cqmom_prepared_case{cqmom_case};
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_prepared_case);
    if (!cqmom.ok())
    {
        std::cerr << "CQMOM failed: " << cqmom.message << '\n';
        return 1;
    }

    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    const qmomkit::CqmomQuadrature& quadrature = state.quadrature();

    return 0;
}
```

## Advanced settings

Most cases need none of the following.

### Retaining extra moments

`initial.maximum_orders` does not affect configured CQMOM loading. The moment
set QMOMKit builds follows from `node_counts` and `conditioning_order` alone, so
supplying the key changes nothing about the case a loader returns.

It is carried through `parse_cqmom_initial_condition(...)` for applications that
build a full raw moment tensor themselves with
`initialize_cqmom_moment_tensor(...)`, which requires it for an
independent-product specification. There, `maximum_orders[d]` gives the largest
raw-moment order retained in dimension `d`:

```text
extent[d] = maximum_orders[d] + 1
```

Parsing checks that the key has one entry per phase-space dimension, that each
order is representable, that the resulting extents do not overflow, and that the
entry count matches the configured marginals. It does not check the orders
against `node_counts`: a tensor too small to contain every moment the hierarchy
requires is rejected when those moments are extracted from it.

### Full moment tensor

Use explicit tensor moments when your application already has a full raw
moment tensor from external data, validation, debugging, or interoperability:

```ini
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
Configured CQMOM loading builds the moment set implied by `node_counts` and
`conditioning_order`; `initial.maximum_orders` does not change it. Sample-based
initialization can use larger extents when evaluating raw moments from its
sample population.

### CQMOM diagnostics and solver options

Most applications can omit the `[cqmom]` and `[vandermonde]` sections. Add
them when you need more or less reconstruction diagnostics, or when a
validation study deliberately compares solver policies:

```ini
[cqmom]
diagnostic_level = detailed|basic|none

[vandermonde]
solver = automatic
rcond = 1.0e-14
```

`diagnostic_level` controls optional reconstruction observability only. Required
numerical checks, solver acceptance, support checks, adaptive
branch selection, inversion results, and success/failure status are unchanged.
If omitted, `none` is used. Select `basic` for structural branch diagnostics or
`detailed` for residual and condition-estimate diagnostics; the latter may add
diagnostic work.

`vandermonde.solver` controls the branch-local conditional Vandermonde solve.
Supported values are:

| Value | Behavior |
|---|---|
| `automatic` | Default policy; uses Björck-Pereyra. |
| `dense` | Uses dense Gaussian elimination. |
| `pivoted_qr` | Uses pivoted QR with numerical-rank detection. |
| `svd` | Uses SVD with rank and convergence checks. |
| `moore_penrose` | Uses a truncated SVD Moore-Penrose pseudoinverse. |
| `bjorck_pereyra` | Uses Björck-Pereyra. |

If omitted, `automatic` is used. The `pivoted_qr`, `svd`, and `moore_penrose`
policies solve the weighted conditional Vandermonde system directly and are
selected only when requested explicitly:

```ini
[vandermonde]
solver = svd
```

`vandermonde.rcond` controls the explicit `moore_penrose` policy. It defaults
to `1.0e-14` and must be finite with `0 <= rcond < 1`, regardless of the
selected solver. The pseudoinverse retains precisely the singular directions
with `sigma > rcond * sigma_max`; discarded directions contribute zero. Its
value is ignored by every other solver. For example:

```ini
[vandermonde]
solver = moore_penrose
rcond = 0.25
```

Rank-deficient systems and SVD nonconvergence remain failures for the regular
`svd` policy. `moore_penrose` accepts a converged rank-deficient SVD and returns
the finite truncated minimum-norm solution. SVD rank and Jacobi tolerances are
internal defaults. Loaded CQMOM cases and states expose the configured policy
through `vandermonde_solver()` and `rcond` through `conditional_options()`.

## Practical guidance

Keep units consistent. Distribution parameters, moment values, and phase-space
coordinate names should refer to the same coordinate convention.

Keep moment normalization explicit. QMOMKit uses unnormalized raw moments of
the number density function. `m0` is the total number density. Dividing by `m0`
changes the setup to a unit-number-density problem.

Generate enough moments. For `N` QMOM nodes, an analytical initialization needs
`maximum_order >= 2 * N - 1`. For normal CQMOM configuration, `node_counts` and
`conditioning_order` define the reconstruction setup. A full tensor of raw
moments, whether given in `[moments]` or built with
`initialize_cqmom_moment_tensor(...)` (which requires `initial.maximum_orders`),
must contain every moment the hierarchy requires; a tensor that is too small is
rejected when those moments are extracted from it.

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
when present. If explicit moments and `[initial]` are both present, the
explicit moments are used.

## Example files

The inline examples in this guide cover explicit moments, mixture
initialization, and runtime realizability settings. Executable examples keep
their complete configurations alongside the corresponding workflow.

Phase-space evolution examples with local configuration files include:

- [examples/qmom_mixture_initialization/configuration.qmom](../../examples/qmom_mixture_initialization/configuration.qmom)
- [examples/qmom_1d_breakage_advection/configuration.qmom](../../examples/qmom_1d_breakage_advection/configuration.qmom)
- [examples/cqmom_2d_bubbles_cavitation/configuration.qmom](../../examples/cqmom_2d_bubbles_cavitation/configuration.qmom)

Monte Carlo companion examples live under
[tools/monte_carlo/configuration/](../../tools/monte_carlo/configuration/):

- [qmom_distribution_lognormal.qmom](../../tools/monte_carlo/configuration/qmom_distribution_lognormal.qmom)
- [cqmom_independent_product_2d.qmom](../../tools/monte_carlo/configuration/cqmom_independent_product_2d.qmom)

Build and run examples as described in the [Installation](installation.md)
guide, then edit the
configuration files to explore different initial moments, node counts, and
output settings.
