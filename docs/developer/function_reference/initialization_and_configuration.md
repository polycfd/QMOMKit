# Initialization, configuration, and policy ownership

Distribution arithmetic is in
[distributions.cpp](../../../src/distributions.cpp), declared in
[distributions.hpp](../../../include/qmomkit/distributions.hpp). Configuration
definitions are in [configuration.cpp](../../../src/configuration.cpp),
declared in [configuration.hpp](../../../include/qmomkit/configuration.hpp).
Public distribution functions belong to `qmomkit`; configuration entry points
to `qmomkit::config`; parsing/formula helpers are local. See [distribution
design](../distributions_design.md), [configuration
design](../configuration_design.md), and [user
configuration](../../user/configuration.md).

## Single-distribution moments — Level A

`compute_raw_moments` dispatches a typed distribution variant through visitors
to local `compute_normal`, `compute_log_normal`, `compute_beta`,
`compute_gamma`, `compute_chi`, `compute_chi_squared`, or `compute_weibull`.
Selected-order variants use the corresponding `*_selected` helper. Inputs
include physical total number density, distribution parameters, and either a
maximum integer order or selected raw indices. Results retain raw density:
`m0=N`, not one unless `N=1` was requested.

`normal_expectations` builds unit-density moments with `E0=1`, `E1=mean`, and
`Ek=mean*E(k-1)+(k-1)*variance*E(k-2)`. Its ascending-order loop depends on the
previous two orders. `compute_normal` multiplies these by `N` afterward and
checks finiteness. Selected normal moments still compute recurrence through the
largest requested order, then gather requested entries. An intermediate
unit-density moment may overflow even when subsequent small-density
multiplication would make the exact raw moment representable.

Lognormal and Weibull use `evaluate_log_expectation` with the local
`evaluate_*_moment` functions:

| Family | Log of unit-density order-`k` moment |
|---|---|
| Lognormal | `k*mu + k*k*sigma*sigma/2` for normal log-coordinate parameters |
| Weibull | `k*log(scale)+lgamma(1+k/shape)` |

The shared evaluator adds `log(N)`, compares with the logs of maximum finite
double and minimum denormal, then exponentiates and checks the result. Order
zero returns `N` directly. `compute_from_log_moments` loops over positive
orders until completion or the first failure.

Beta, gamma, chi-squared and chi build each moment as a product of one factor
per order in a local `ScaledProduct`, which starts from `N` and keeps the
running product as a mantissa and a binary exponent:

| Family | Factor from order `k-1` to `k` |
|---|---|
| Beta | `(alpha+k-1)/(alpha+beta+k-1)` (`multiply_beta_step`) |
| Gamma | `shape+k-1`, then `scale` (`multiply_gamma_step`) |
| Chi-squared | `nu+2(k-1)` (`multiply_chi_squared_step`) |
| Chi | `nu+k-2` from order `k-2` (`multiply_chi_step`); order one is `sqrt(2)*half_step_gamma_ratio(nu/2)` (`chi_first_moment`) |

`ScaledProduct::value` reports overflow when the exponent exceeds the double
range and underflow when the value rounds to zero, so both are decided exactly.
`compute_from_products` extends one product order by order; the chi
maximum-order path keeps separate even and odd products. Selected variants
evaluate each requested order with `evaluate_product_moment` or
`evaluate_chi_moment`, multiplying the same factors in the same order, and
retain their order and duplicates. `half_step_gamma_ratio` evaluates
`Gamma(x+1/2)/Gamma(x)` from an eight-term asymptotic expansion from `x = 10`
on and raises smaller `x` by the recurrence `r(x) = r(x+1) x/(x+1/2)`.

Checked preconditions include finite nonnegative density, positive finite
shape/ scale/standard-deviation parameters where required, and capacity-safe
maximum order. Zero density produces zeros once the relevant validation path
has accepted its request. Failures distinguish invalid input, overflow,
underflow, and nonfinite arithmetic, with failed order/component/dimension
where available. Normal arithmetic and logarithmic arithmetic do not have
identical underflow-detection paths.

Evidence: [distribution unit tests](../../../tests/unit/distributions.cpp) use
analytical moments, density scaling, selected/full parity, invalid parameters,
and overflow/underflow assertions. [Distribution verification](../../../tests/verification/distributions/distributions_analytical.cpp)
compares the product families with error-free references for parameters up to
`1e8` and odd chi moments with high-precision values; these cases cover orders
up to ten and do not establish a uniform envelope in order or density.

## Mixtures and independent products — Level A

The mixture `compute_raw_moments` overload validates component densities and
requires at least one positive component. It skips zero-density components,
including their distribution-parameter evaluation. For each remaining
component, it obtains already density-weighted moments and accumulates each
requested order. Local `add_compensated` and `add_selected_compensated` use a
Neumaier-style correction chosen by relative magnitude of the accumulated sum
and new term. Both intermediate accumulation and final `sum+correction` are
checked. Failure identifies the contributing component when available.
Compensation reduces some summation error; it does not prevent overflow or
guarantee relative accuracy near cancellation of odd moments.

The independent-product overloads implement
`M_alpha = N * product_d (marginal_m[alpha[d]] / marginal_m[0])`.
This **explicit normalization of marginal shapes** is mathematically distinct
from normalizing away the global number density. Each marginal must have finite
positive zeroth moment. `normalize_marginal_moment` checks division, including
nonzero underflow. The per-index coordinate loop checks multiplication overflow
and nonzero-to-zero underflow, starting from global `N` exactly once.

Full-tensor generation checks maximum-order rank, constructs checked extents, then
loops over all row-major tensor indices. Selected generation gathers the unique
orders needed in each coordinate (including zero for normalization), computes
only those marginal requests, and forms each requested mixed moment in caller
order. Helpers `contains_order`, `order_position`, and
`first_requested_index_with_order` preserve request-to-failure correspondence.
The independence assumption belongs to the specified distribution; no
correlated model is inferred. Zero global density returns a zero result before
marginal numerical evaluation, so bad unused marginal parameters can be
bypassed on this branch.

Evidence: mixture and product assertions in [distribution
tests](../../../tests/unit/distributions.cpp) compare component sums,
normalized marginal products, tensor/selected outputs, and failure attribution.
They do not characterize many-component cancellation or extreme unequal
marginal scales. The [minimal CQMOM
example](../../../examples/cqmom_2d_minimal/README.md) uses independent-product
initialization through the public configuration route; it is workflow evidence,
not an automated test.

## Parsing structure and precedence — Levels A/B

Local `parse_file` is a line-oriented parser. It strips `#` comments, trims,
recognizes sections, then stores qualified scalar keys or explicit `M[...]`
entries. Duplicate scalar/moment keys become errors with source/line
information. It is not a quoted-string language: a `#` is treated as a comment
delimiter. Section occurrences are recorded, but repetition is rejected only
where subsequent validators require it. `parse_size` uses unsigned conversion
with a negative-sign check; `parse_double` requires full finite conversion.
List parsers reject empty entries/trailing commas; optional-bound lists
additionally recognize `none`.

`parse_qmom`, `parse_cqmom`, and `parse_cqmom_prepared` each
call `parse_file` exactly once and share the resulting `ParsedFile` with the
initial-condition step, instead of re-parsing the same text a second time.
`parse_qmom_initial_condition_from_parsed` and
`parse_cqmom_initial_condition_from_parsed` are internal-linkage siblings of
the public `parse_qmom_initial_condition` and `parse_cqmom_initial_condition`
that take an already-parsed `ParsedFile` instead of raw text; the public
functions remain independently callable and still parse exactly once each. The
CQMOM case parsers pass a copy into `parse_cqmom_case_settings` so its
validation cannot affect the `ParsedFile` handed to the initial-condition step
afterward; `parse_qmom` instead moves its own already-validated, still
error-free `ParsedFile` directly, which is equivalent to a fresh parse at that
point.

`parse_qmom_initial_condition` and `parse_cqmom_initial_condition` choose
explicit moments when present. That branch intentionally takes precedence over
analytical initialization sections; selected otherwise-unused initialization
keys are allowed rather than overriding explicit data. Syntax and duplicate-key
errors are still errors. Without explicit data, `parse_distribution`,
`parse_mixture`, and `parse_independent_product` build typed specs. Numbered
mixture components and product dimensions must be contiguous from zero.
Family-specific allowed-key checks prevent parameters for one family silently
configuring another.

`parse_explicit_tensor` validates rank/extents, checks extent-product overflow,
allocates a supplied-entry mask, and loops over explicit moments to validate
rank, bounds, finiteness, and complete coverage. Missing entries are errors,
not zeros. Explicit tensor moments contain all tensor entries; extraction
later selects the canonical required subset.

Evidence: [configuration tests](../../../tests/unit/configuration.cpp) exercise
unknown keys/sections, duplicates, explicit-data precedence, incomplete
tensors, family parsing, and invalid shapes. Parser acceptance alone does not
establish realizability of the initialized moments.

## Mapping numerical policy — Level A/B

`parse_cqmom_case_settings` resolves counts, names, conditioning permutation,
reconstruction policies, and paired support intervals. When omitted, diagnostics
default to none, the solver to automatic, and the conditioning order to
identity. Adaptive branches enable the recursive adaptive workflow; a standalone
adaptive Wheeler request without adaptive branches is rejected. Names are
metadata and do not establish units or support.

`parse_qmom_options_from_parsed` maps `qmom.adaptive`; other Wheeler thresholds
are API options, not arbitrary accepted `.qmom` keys.
`parse_automatic_scaling_from_parsed` separately maps
`phase_space.automatic_scaling`. Case parsers and
`parse_qmom_phase_space_evolution_options` validates the shared configuration,
but only the case parser places reconstruction policy into the `QmomCase`.
`QmomCase::initial_state()` transfers that policy and the configured maximum to
the evolving `QmomState`; phase-space evolution options contain only per-call
controls. CQMOM evolving state likewise retains its case reconstruction
options.

`parse_coordinate_support` pairs lower/upper lists in phase-space order: each
coordinate has both finite bounds with lower below upper, or neither. Bounds
are node admissibility controls, separate from support-specific realizability
policy. `parse_realizability_policy_from_parsed` maps support, tolerance,
boundary acceptance, and check mode.
`parse_coordinate_realizability_policies_from_parsed` constructs optional
per-coordinate entries; case-level rank validation rejects policies outside the
hierarchy. These policies still only check marginals.

`parse_node_advection_options_from_parsed` maps scheme and embedded controls,
requires positive configured controls, and checks minimum/maximum consistency
when supplied. API zero-as-default controls are not equivalent to explicitly
configured zero values. `parse_phase_space_evolution_scheme_from_parsed` and
`parse_*_phase_space_evolution_options` package operator order and numerical
options; they never install physics callbacks or own a simulation loop.
`parse_output_options` returns writer preferences, not an open writer.

Evidence: configuration unit tests inspect defaults, explicit overrides, enum
rejection, support pairing, and coordinate policy rank. The
[QMOM advection example](../../../examples/qmom_1d_breakage_advection/main.cpp)
loads case and evolution options separately. The configuration unit regression
covers true, false, and omitted
automatic-scaling values across the QMOM evolution loader.

## Initialization, construction, and replacement — Levels A/B

`initialize_qmom_moments` copies an explicit sequence or computes analytical
moments; explicit copying does not independently invert or classify it.
`initialize_cqmom_moment_tensor` copies an explicit tensor or evaluates the
full independent product, with maximum orders required to set the tensor
extents. `initialize_cqmom_moments` first obtains the hierarchy's
required layout, then extracts explicit data or requests selected analytical
moments. Configured case construction uses required indices rather than
computing an unnecessary full tensor, so it neither needs nor consults
`initial.maximum_orders`.

`create_qmom_case` and `create_cqmom_case` validate structural
suitability and package initialized data with numerical settings. They do not
guarantee that `initial_state` reconstruction will succeed. The QMOM's
large-count guard protects its own requirement; CQMOM moment extraction rejects
a tensor that does not contain every required moment. `with_initial_moments` validates
replacement tensor structure/finiteness and required hierarchy coverage,
preserves case settings/prepared metadata, and defers realizability/inversion.
It extracts the moments the hierarchy requires from the supplied tensor. Its
`CqmomMomentData` overload takes that required moment set directly, validating
that the supplied layout is the hierarchy's required layout with one finite
value per index. Both overloads yield the same case for the same moments; the
tensor overload exists for moment data arriving from outside QMOMKit, where
reporting a malformed tensor as a configuration error is more useful than
throwing. The `CqmomCase` overloads hold the validation. The
`CqmomPreparedCase` overloads call them, pass a failure on unchanged, and give
the replaced case the base case's prepared metadata through
`with_prepared_metadata_of`, so both kinds of base case report the same errors.

Public `parse_qmom`, `parse_cqmom`, and `parse_cqmom_prepared` compose parser,
initializer, and constructor; `parse_cqmom_prepared` additionally prepares the
case. Each `load_*` reads a file then
invokes its corresponding parser with the filename as diagnostic source. The
reader does not enforce a filename extension; `.qmom` is the documented
configuration convention.

Evidence: [configuration tests](../../../tests/unit/configuration.cpp) and
[case-replacement
regression](../../../tests/component_regression/cqmom_prepared_case_moment_replacement.cpp)
cover preserved settings, layout compatibility, deferred inversion, and
replacement through both initial-moment overloads. [Explicit QMOM
example](../../../examples/qmom_explicit_moments/README.md) shows the
distinction between specified raw moments and reconstructed nodes.

## Helper families — Level C

| Family | Purpose and important semantics | Evidence |
|---|---|---|
| Distribution `check_*`, `is_positive_finite`, `finish_*`, `zero_selected_moments`, result `ok` | Validate family/request domain and package completed raw values; zero shortcuts differ by layer as described above. | Distribution invalid/zero tests |
| Distribution `make_*failure`, `invalid_*`, `failure_at_order`, `finite_failure_status`, `selected_failure_from_moment` | Preserve order/component/dimension attribution; distinguish underflow, overflow, invalid input, and nonfinite arithmetic. | Failure-attribution tests |
| Distribution visitors and `compute_independent_product_marginal_moments` | Dispatch typed single/mixture inputs to the same checked evaluators. | Family/product parity tests |
| Parser `require_entry`, `find_entry`, `parse_required_*`, `parse_optional_*` | Required-key diagnostics and default-preserving optional conversion; no physical defaults inferred. | Missing/default/invalid configuration tests |
| Parser key/section recognizers, allowed-key builders, `reject_unknown_*`, `validate_optional_*`, enum parsers | Define accepted configuration surface and prevent unrelated keys changing runtime state. | Unknown-key and enum tests |
| `trim`, `parse_moment_key`, `mixture_component_prefix`, `validate_names`, phase-space/version validators | Text/index syntax and numerical rank metadata checks. | Parser syntax/rank tests |
| `checked_product`, `validate_replacement_tensor_structure`, `format_moment_index` | Checked tensor allocation/shape boundary and precise missing-index diagnostics. | Tensor/replacement tests |
| `add_error`, `format_configuration_error(s)`, `failure_result`, `file_failure`, `read_configuration_file` | Accumulate source/line errors and convert file access failure; do not partially return a valid case. | File/parse failure tests |
| `ConfigurationResult` constructors, boolean/status/value/error accessors | Own an optional value and error string; `value()` without a value throws `std::bad_optional_access`, while `error()` always returns the stored string. | Configuration tests and installed-header coverage |

The [inventory](coverage_inventory.md) lists exact helper spellings and overloads;
shared rows here describe their common responsibility, not a claim of identical
validation across different input domains.
