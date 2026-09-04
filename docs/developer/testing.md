# Testing

This page documents the current QMOMKit testing philosophy and test
architecture for maintainers, contributors, and AI-assisted development
workflows. It describes stable correctness expectations and the current root
CTest inventory, not temporary CI setup.

For package layout and consumer boundaries, see `architecture.md`. For moment,
tensor, and phase-space conventions that tests must protect, see
`numerical_conventions.md`. For realizability status behavior and runtime
acceptance policies, see `realizability.md`. For Monte Carlo companion package
validation guidance, see `monte_carlo.md`.

## Debug Build Requirement

The unit tests use C/C++ `assert()`. `assert()` is disabled when `NDEBUG` is
defined, and normal Release builds define `NDEBUG`; therefore a Release CTest
run does not validate unit assertions. The authoritative functional test
configuration is Debug:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DQMOMKIT_BUILD_TESTS=ON

cmake --build build --parallel

ctest --test-dir build --output-on-failure
```

Run one registered unit contract with:

```bash
ctest --test-dir build -R '^unit_wheeler$' --output-on-failure
```

Keep the existing assertion-based unit tests intact. Do not use a Release CTest
result as evidence that the functional unit suite passed.

## Test-suite Organization

QMOMKit is a numerical C++ library, so correctness requires more than exercising
code paths. A passing test should usually demonstrate a mathematical invariant,
a known analytical result, a manufactured solution, or a documented failure
contract.

Unit tests alone are insufficient for moment-method software. They are good at
protecting local behavior, option validation, indexing helpers, and failure
paths, but they do not by themselves establish that an inversion, conditional
reconstruction, phase-space evolution update, or realizability check matches the intended
mathematics.

The core testing principle is:

```text
QMOMKit testing

                    Test purpose
                         |
        -----------------------------------
        |                 |               |
      unit        component_regression   verification
        |                 |               |
 isolated correctness   implementation   mathematical and
 local APIs             stability        workflow correctness
```

Analytical and manufactured solutions are especially important because they make
the expected answer independent of the implementation under test. When a
verification compares recovered quadrature, evolved moments, or reconstructed
raw moments, the reference values should come from the mathematical case being
tested, not from another call into the same production algorithm.

When adding a new test, choose the directory by purpose:

1. If it verifies an isolated function, class, or local API contract, use
   `tests/unit/`.
2. If it protects specific implementation behavior or a regression-prone
   component interaction, use `tests/component_regression/`.
3. If it verifies a mathematical property or user workflow against an
   analytical, manufactured, convergence, or otherwise independent reference,
   use `tests/verification/`.

### Test Project Structure

The primary test project lives under:

```text
tests/
```

It is included by the root CMake project when `QMOMKIT_BUILD_TESTS=ON`. Normal
unit, component regression, verification, and public-header tests link the
current build-tree `qmomkit` target and compile against the current source-tree
public headers. Installation is not required for the normal correctness suite.

### Running Tests

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DQMOMKIT_BUILD_TESTS=ON

cmake --build build --parallel

ctest --test-dir build --output-on-failure
```

### CTest Naming and Current Count

CTest names identify the tested contract. Several executables accept a local
selector so independently meaningful groups can be registered separately. The
individual CTest registrations are authoritative; running such an executable
without a selector is a developer convenience, not an aggregate CTest contract.
The current root suite has **113** registered CTest tests: 32 unit, 7
component-regression, 73 verification, and 1 supporting public-header test.
Obtain the current count from `ctest -N` (and confirm structured registration
with `ctest --show-only=json-v1`), rather than maintaining it by arithmetic
alone.

The test project currently registers CTest executables through these
subdirectories:

```text
tests/
    unit/
    component_regression/
    verification/
    public_headers/
```

Some repository tests include `src/` as a private include path for focused
internal validation. That is permitted inside the repository test project, but
it does not change the public API boundary. Downstream users should consume only
installed headers and exported package targets.

Package consumption is also tested explicitly by:

```text
tests/package_consumer/
```

That directory is a minimal external CMake project. It validates installed
package discovery, exported target linkage, and installed public-header usage
without repository source include paths or `add_subdirectory`.

## Test Categories

### Unit Tests

Unit tests verify isolated functions, classes, data structures, options, and
local API contracts. They should be focused, fast, deterministic, and specific
about the contract being protected.

They intentionally use `assert()` and therefore require a Debug build.

Unit coverage is appropriate for:

- request and option validation;
- status values, diagnostics, and `ok()` behavior;
- indexing and storage helpers;
- local algorithm invariants;
- error handling and state preservation;
- boundary cases that should fail cleanly.

Representative unit areas include Wheeler inversion, recursive CQMOM support
structures, realizability checks and policies, configuration mapping,
distribution helpers, output helpers, node advection, moment-source evolution,
raw-moment reconstruction, and phase-space evolution components.

### Component Regression Tests

Component regression tests protect important implementation behavior that is
larger than an isolated unit contract but narrower than mathematical
verification. They are appropriate for direct component composition,
diagnostics, status propagation, policy mapping, and previously fragile
implementation paths that should remain visible during maintenance.

Component regression coverage belongs under:

```text
tests/component_regression/
```

These tests may use manufactured inputs or simple analytical expectations, but
their primary purpose is regression protection. They are not failed
verification tests, not a substitute for unit tests, and not replacements for
verification tests when the behavior is a public mathematical workflow with an
independent reference.

Independent numerical verification requires an analytical, manufactured,
convergence, conservation, or otherwise independent reference. Dense/lean or
batch/scalar parity, round trips, captured fixtures, and comparisons between
related production paths remain valuable consistency/regression evidence, but
must not be presented as independent algorithm verification.

### Verification Tests

Verification tests provide mathematical evidence that the numerical algorithms
behave as intended. They compare against analytical solutions, manufactured
solutions, reconstruction identities, convergence behavior, or other
independently known results.

Verification tests serve as executable documentation of mathematical
correctness and public workflows.

Verification coverage belongs around:

- QMOM and Wheeler reconstruction;
- CQMOM conditional moments and recursive quadrature;
- distribution-generated raw moments;
- support-aware realizability checks;
- Runge-Kutta integration;
- node advection;
- moment-source evolution;
- phase-space evolution composition.

Verification tests should preserve QMOMKit conventions: raw moments are moments
of the number density function, tensor indices follow phase-space dimension
order, conditioning order does not redefine phase-space order, and quadrature
weights are not normalized unless the case explicitly defines unit total number
density.

### Public-Header Tests

Public-header tests are compile-time interface checks. Each translation unit in
`tests/public_headers/` includes one current public header from
`include/qmomkit/` through the source-tree `qmomkit` target.

These tests protect the public API boundary by verifying that headers are
self-contained and do not depend on accidental include order. The linked CTest
executable does not need numerical runtime assertions; successful compilation
and linkage are the meaningful checks.

When a public header is added, renamed, or given new direct dependencies, the
public-header coverage should be updated so downstream include behavior remains
visible.

### Package Consumer Tests

Package consumer tests validate the installed-package experience rather than the
in-tree build. The consumer project should behave like an external application:

- discover QMOMKit with `find_package(QMOMKit CONFIG REQUIRED)`;
- link only `qmomkit::qmomkit`;
- include installed public headers;
- avoid repository source paths, internal headers, and test utilities.

`tests/package_consumer/` is intentionally separate from the normal CTest unit,
verification, and public-header suite. It represents an external consumer
workflow and requires QMOMKit to be installed before the consumer project is
configured. Its purpose is to validate installed package discovery and exported
CMake target behavior, not correctness or regression behavior inside the normal
`tests/` project.

This category is most important when package exports, installation rules,
public headers, or CMake package metadata change.

## Numerical Validation Principles

Analytical solutions are preferred when the exact answer is known in closed
form. They are useful for Wheeler inversion, low-dimensional CQMOM cases,
distribution moment generation, simple ODE integration, moment-source evolution, and
phase-space evolution cases with closed-form operator behavior.

Manufactured solutions are useful when the algorithm is more complex than the
available closed-form examples. They should define moments, nodes, weights,
conditioning paths, source terms, or phase-space evolution behavior from an independent
construction and then verify that QMOMKit recovers the expected mathematical
object.

Convergence studies are useful for time integration and phase-space evolution. They should
test numerical order or monotonic error reduction using deterministic cases
whose exact or high-confidence reference behavior is known.

A convergence registration must identify its independent reference, refinement
variable and sequence, error quantity, expected order or trend, and pass
threshold. Tolerance refinement for adaptive methods is accuracy evidence; it
is not a formal fixed-step order proof.

Deterministic regression cases are useful for boundary behavior, diagnostics,
and previously fragile numerical situations. They should preserve the meaningful
status, tolerance, ordering, or reconstruction property rather than simply
locking in implementation-shaped output.

High-dimensional CQMOM needs special care. Tests should separate phase-space
dimension order from conditioning order, verify terminal quadrature weights and
coordinates, and reconstruct raw tensor moments from the final quadrature when
that is the mathematical guarantee being tested.

Realizability tests should cover accepted, boundary, rejected, invalid-input,
and numerically indeterminate cases. They should also verify that policies reject
candidate moment states without implying automatic correction.

## Adding New Tests

When adding functionality, add tests at the level where the new behavior is
owned.

Add unit tests for local behavior: request validation, option defaults, status
values, diagnostics, indexing, state preservation, and explicit failure paths.

Add verification tests for mathematical guarantees: moment reconstruction,
closed-form evolution, manufactured quadrature recovery, support-aware
realizability, convergence, or invariant preservation.

Add public-header coverage when public APIs change. Installed headers should
compile when included directly and should not rely on a caller including another
header first.

Add or update package-consumer coverage when installation, export, package
configuration, or public target behavior changes.

The optional Monte Carlo companion package has its own CMake project and tests
under `tools/monte_carlo/`. Configure that project with `BUILD_TESTING=ON` when
running its tests, then run them explicitly with CTest.

Update examples only when the user-facing workflow changes. Examples are not a
substitute for automated correctness tests.

Prefer deterministic tests. Stochastic or Monte Carlo workflows are valuable for
independent validation, but normal correctness tests should not depend on random
sampling noise unless the uncertainty model is explicit and controlled.

## Relationship With Benchmarks

Tests establish correctness. Benchmarks measure performance.

A benchmark should not be the only evidence that an algorithm is correct, and a
correctness test should not depend on timing-sensitive behavior. Performance
work should preserve the same mathematical contracts checked by unit and
verification tests.

Benchmark structure and interpretation belong in `benchmarks.md`.

## Registered Root CTest Inventory

This inventory lists every root CTest registration exactly once. It is current
at the documented count of 114 tests.

### Unit Tests (32)

| Unit test | What it tests | How correctness is checked | Important cases covered |
|---|---|---|---|
| `unit_wheeler` | One-dimensional Wheeler inversion and adaptive node reduction | Uses hard-coded quadrature-derived moments; compares recovered nodes, weights, ordering, and reconstructed moments, plus explicit statuses | One to three nodes, scaling, degenerate/non-realizable sequences, adaptive reduction, negative mass/variance, non-finite input, invalid options |
| `unit_runge_kutta` | RK23, RK45, and SSP-RK3 contracts | Checks published tableau coefficients, independently evaluated scalar steps, callback/stage traces, and exact status/counter results | FSAL, forward/backward integration, retries, minimum/maximum step, maximum steps, invalid options, stage predicates |
| `unit_recursive_indexing_layout` | CQMOM indexing, hierarchy, conditioning paths, and branch-local layouts | Compares hand-specified multi-indices, flattened indices, offsets, paths, and failure statuses | 1D--3D, unequal counts, identity/non-identity order, invalid hierarchy/index, branch diagnostics |
| `unit_recursive_lean_reconstruction` | Lean CQMOM moment storage and recursive reconstruction interface | Checks required layouts and compares lean reconstruction/state mapping with dense representation; this establishes representation consistency, not an independent mathematical oracle | 2D/3D, non-identity order, unequal/adaptive counts, missing/non-finite/incompatible data |
| `unit_recursive_bounded_support` | Bounded-support recursive CQMOM reconstruction | Uses explicit feasible/infeasible intervals; checks returned nodes remain within bounds without clipping and verifies fixed/adaptive statuses | Shifted intervals, coordinate-specific support, fixed/adaptive reduction and growth, dense/lean routes |
| `unit_recursive_solver_policy` | Conditional Vandermonde policy selection inside recursive CQMOM | Uses manufactured reconstructions and compares quadrature/conditional values and diagnostics across solver policies | Dense, Björck--Pereyra, fallback, pivoted QR, SVD, automatic selection, failure diagnostics |
| `unit_conditional_moments` | Conditional-moment sequence recovery | Hand-composes mixed raw moments from known parent nodes, weights, and conditional moments before numerical recovery | Non-identity conditioning, unused coordinates, singular/nearly singular systems, warnings/failures |
| `unit_moment_rates` | Generic QMOM/CQMOM moment-rate interfaces | Uses callbacks with known rates and checks callback ordering, layout/index forwarding, status propagation, and output values | QMOM/CQMOM, callback errors, diagnostics, source contexts |
| `unit_advection_moment_rates` | Quadrature-derived advection moment rates | Independently sums terminal-node velocity contributions and compares returned rates | QMOM/CQMOM, non-identity order, zero/linear velocity, adaptive topology, callback failures |
| `unit_moment_source_rates` | Direct source-rate evaluation | Uses explicit callback outputs and invocation expectations; checks rates, layout/index forwarding, and error/status propagation | Malformed sizes, non-finite values, throwing callbacks, QMOM/CQMOM contexts |
| `unit_bjorck_pereyra` | Transposed-Vandermonde Björck--Pereyra solve | Solves explicit manufactured systems with known coefficient matrices | Multiple/zero RHS, clustered nodes, repeated-node rejection, invalid dimensions/finite values/tolerances |
| `unit_householder_qr` | Pivoted Householder QR solve | Checks known solutions, direct residual and normal-equation conditions, factor reuse, and pivot restoration | Square/tall/inconsistent systems, rank deficiency, invalid input |
| `unit_jacobi_svd` | Jacobi SVD factorization and solve | Checks singular values, orthogonality, matrix reconstruction, rank, and known solutions | Rank deficiency, scaling, ill-conditioned systems, invalid matrices |
| `unit_vandermonde` | Weighted conditional Vandermonde dispatch and policies | Uses manufactured coefficients, componentwise solution checks, scale-normalized residuals, and selected-solver/fallback diagnostics | All policies, fallbacks, clustered nodes, zero RHS, larger systems, invalid input; the large-dynamic-range case uses resolvable small coefficients with relative-plus-near-zero checks and an `A`/`x`/`b`-scaled residual |
| `unit_public_api` | Public QMOM/CQMOM data, layout, and reconstruction contracts | Checks exact layout/index/value expectations, dense/lean and 1D-CQMOM/QMOM mapping, and invalid-input failures | Missing/extra/non-finite data, adaptive constraints, diagnostics |
| `unit_configuration_initialization` | Configuration-driven QMOM/CQMOM initialization | Compares parsed cases and generated moments with direct/programmatic construction and analytical initialization APIs | Distributions, mixtures, independent products, dense/lean CQMOM, support settings |
| `unit_configuration_cqmom_case_state` | `CqmomCase` and `CqmomState` contracts | Checks stored settings, moment data, prepared metadata, state compatibility, and replacement behavior | Dense/lean storage, support, immutable metadata, state compatibility |
| `unit_configuration_runtime_options` | Runtime `.qmom` option mapping | Checks exact defaults and parsed enum/numeric values plus explicit invalid-option rejection | Diagnostics, solver, realizability, splitting, node-advection schemes, adaptive RK, output |
| `unit_configuration_failure_paths` | Malformed and invalid configuration behavior | Uses intentionally malformed strings/files and explicit failure expectations | Missing/duplicate/unknown fields, non-finite values, incompatible options, invalid dimensions/ranks, ambiguous/missing files |
| `unit_output` | Moment and quadrature output writers | Reads generated temporary files and compares headers, columns, ordering, and numerical values | Selected/duplicate requests, adaptive histories, CQMOM layouts, malformed input |
| `unit_distribution_selected_results` | Analytical distribution raw-moment requests | Compares selected results with known values and dense generation; checks diagnostic/failure contracts | Request ordering, duplicate orders, rank, invalid parameters, overflow/underflow, zero density |
| `unit_distribution_mixtures` | Mixture raw-moment generation | Compares explicit component-weighted sums with generated moments and checks component-index diagnostics | Heterogeneous families, zero-density components, empty/invalid mixtures, accumulation failures |
| `unit_distribution_independent_products` | Multidimensional independent-product moment generation | Checks product moments, tensor flattening, selected/dense consistency, and normalization | Rank/order validation, zero global density, mixture marginals, per-dimension diagnostics |
| `unit_realizability` | Raw-moment realizability diagnostics | Checks exact status, failed matrix/order, pivot metadata, and boundary classification for constructed sequences | Real-line/non-negative/bounded support, invalid input, mass policy, boundary/indeterminate cases |
| `unit_realizability_policy` | Runtime realizability-policy adapters | Checks predicate/converter invocation counts and accept/reject/retry behavior | Disabled, final-only, every-stage, support forwarding, bounded interval, converter failures |
| `unit_node_advection` | QMOM/CQMOM node advection and integrator selection | Uses analytically predictable velocity fields, expected node updates, preserved weights/topology, callback counts, and diagnostics | QMOM/CQMOM; RK23/RK45/SSP-RK3; adaptive topology; callback failures; option propagation. Direct CQMOM RK45 and SSP-RK3 cases use `xi(t + dt) = xi(t) + v dt` and check terminal paths and hierarchy counts |
| `unit_moment_source` | Moment-source evolution and source-context contracts | Checks known source updates, retries/rollback, callback context access, and policy/error expectations | Generalized moments, boundary contributions, realizability policies, invalid values, callback exceptions |
| `unit_moment_reconstruction` | Raw/generalized moment reconstruction from QMOM/CQMOM quadratures | Uses manufactured quadratures and independently summed moments, a finite-difference rate check, and explicit domain failures | Selected/duplicate moments, fixed/adaptive CQMOM, non-identity order, invalid powers/layouts; a negative CQMOM abscissa with fractional exponent throws `std::domain_error`, while the same coordinate with integer exponent remains valid |
| `unit_qmom_phase_space_evolution` | QMOM phase-space evolution orchestration | Compares staged operations with composed public operations and checks state, diagnostics, and rollback; composition parity tests orchestration, not convergence | Advection, source, boundary, Lie/Strang composition, adaptive counts, policies/failures |
| `unit_cqmom_phase_space_evolution` | CQMOM phase-space evolution orchestration | Checks dense/lean state equivalence, topology transitions, source contexts, support/policy behavior, and rollback | Fixed/adaptive branches, support bounds, solver policy, source/boundary callbacks, failures |
| `unit_phase_space_evolution_batch` | Batch phase-space evolution APIs | Compares each batch result with its scalar API and separately checks ordering, isolation, and callback indexing; scalar parity is a dispatch oracle, not an independent numerical reference | Empty batches, indexed callbacks, normal failures, escaping exceptions, state isolation |
| `unit_performance_instrumentation` | Optional performance instrumentation | Checks exact counter values and enabled/disabled semantics while verifying instrumentation preserves numerical routes/results | Solver/source/advection counters, metadata/layout reuse, fixed/dynamic state paths, counter isolation |

### Component-Regression Tests (7)

| Component-regression test | What it verifies | How correctness is established | Important cases / criterion |
|---|---|---|---|
| `component_regression_moment_source` | Moment-source component interactions | Deterministic source-update, status, and diagnostic regression contracts | Callback errors, retry/rollback, and propagated diagnostics |
| `component_regression_node_advection` | Node-advection component interactions | Deterministic manufactured quadratures exercise node advancement, reconstruction, and failure paths | QMOM/CQMOM routes, adaptive topology, and callback failure propagation |
| `component_regression_phase_space_evolution` | Operator-split phase-space-evolution composition | Deterministic state, stage-diagnostic, topology, and rollback contracts across composed workflows | Lie/Strang staging, source/advection interaction, dense/lean state handling |
| `component_regression_cqmom_compute_failures` | CQMOM failure/status propagation | Deliberately invalid or fragile dense/lean inputs are checked for the documented status, diagnostics, and preserved failure behavior | Invalid tensors, conditional reconstruction failures, and adaptive/fixed policy failures |
| `component_regression_cqmom_case_moment_replacement` | Immutable CQMOM case moment replacement | Replaces configured initial tensors and checks retained policies, storage, and state construction contracts | Dense/lean cases, support/options retention, invalid replacement data |
| `component_regression_cqmom_dense_lean_parity_3d` | Dense and lean CQMOM representation/workflow consistency | Compares layouts, reconstructed states, quadratures, and diagnostics across routes sharing the same reconstruction machinery | Consistency/regression evidence only; it is not an independent mathematical CQMOM oracle |
| `component_regression_cqmom_physical_transport_probe` | Historical adaptive CQMOM ragged reinversion/state regression | Replays captured fixture data whose original generation provenance is not recoverable; it is therefore non-independent regression evidence | Four expected terminal paths recover, each final `u` branch reduces from two requested nodes to one, and active counts persist through two zero-rate state advances |

### Verification Tests (74)

| Verification test | What it verifies | How correctness is established | Important cases / criterion |
|---|---|---|---|
| `verification_qmom_wheeler_analytical` | Wheeler reconstruction of prescribed one-, two-, and three-node rules | Independently specified discrete rules generate raw moments; recovered ordered nodes, weights, and reconstructed moments are compared with those rules | One--three nodes, input ordering, mass and moment recovery |
| `verification_qmom_adaptive_reconstruction` | Public adaptive QMOM reduction to the supported node count | Six moments are independently summed from `(-0.75, 0.4)` and `(1.25, 0.6)`; a three-node request must recover that two-node rule | Nodes, weights, mass, and all six moments agree to about `1e-12` |
| `verification_qmom_compute_failures` | Public QMOM failure/status behavior | Deliberately impossible or invalid moment sequences are required to return the documented status/message contract | Failure/status regression contract; not independent numerical verification |
| `verification_qmom_analytical_single_dirac` | One-node QMOM for a Dirac measure | Prescribed node and weight define exact raw moments and the known one-node rule | Node, weight, and moments recover within `1e-12` |
| `verification_qmom_analytical_two_node_discrete` | Two-node discrete QMOM | An independently prescribed two-node rule supplies moments and the expected quadrature | Ordered nodes, weights, and moments recover within `1e-12` |
| `verification_qmom_analytical_three_node_discrete` | Three-node discrete QMOM | An independently prescribed three-node rule supplies moments and the expected quadrature | Ordered nodes, weights, and moments recover within `1e-12` |
| `verification_qmom_analytical_uniform_n2` | Two-node quadrature for a uniform density | Mapped Gauss--Legendre nodes/weights and closed-form uniform moments are the reference | Known rule and moments recover within `1e-12` |
| `verification_qmom_analytical_uniform_n3` | Three-node quadrature for a uniform density | Mapped three-point Gauss--Legendre rule and closed-form uniform moments are the reference | Known rule and moments recover within `1e-12` |
| `verification_qmom_analytical_exponential_n2` | Two-node quadrature for an exponential density | Gauss--Laguerre rule and exponential raw moments are independently known | Nodes, weights, and moments recover within `1e-12` |
| `verification_qmom_analytical_exponential_n3` | Three-node quadrature for an exponential density | Three-point Gauss--Laguerre rule and exponential raw moments are independently known | Nodes, weights, and moments recover within `1e-12` |
| `verification_qmom_analytical_gaussian_n2` | Two-node quadrature for a standard Gaussian | Probabilists' Gauss--Hermite rule and Gaussian moments are the reference | Symmetric nodes/weights, vanishing odd moments, and recovery within `1e-12` |
| `verification_qmom_analytical_gaussian_n3` | Three-node quadrature for a standard Gaussian | Three-point probabilists' Gauss--Hermite rule and Gaussian moments are the reference | Symmetry, zero odd moments, nodes, weights, and moments recover within `1e-12` |
| `verification_cqmom_recursive_quadrature_manufactured` | Fixed and ragged recursive CQMOM construction | Test-local terminal rules independently generate raw tensors before recursive inversion | Topology, paths, local/global weights, abscissas, and reconstructed moments are compared with the prescribed rules |
| `verification_cqmom_quadrature_bounds` | Bounded CQMOM reconstruction, including multidimensional recursion | Feasible/infeasible 1D rules are explicit. A 2D four-terminal rule uses conditioning `{1,0}`, x support `[1,4]`, and y support `[-2,0]` | 2D raw moments reproduce to `1e-12`; every terminal obeys its physical-coordinate support and branch diagnostics report y bounds at the root and x bounds on child branches |
| `verification_cqmom_adaptive_branch_growth` | Adaptive CQMOM branch growth/reinjection through native source evolution | A constant independently prescribed moment-rate changes a `{2,2}`-maximum state from paths `{0,0},{1,0}` to a known three-terminal rule `{0,0},{1,0},{1,1}` | Checks growth transition, active counts, paths, nodes, weights, and target moments to about `5e-11`; an increased count alone is not sufficient |
| `verification_cqmom_recursive_validation` | Fixed and adaptive recursive quadrature recovery | Independently built fixed and ragged terminal trees generate the input moments | Branch diagnostics, topology/paths, weights, abscissas, mass closure, and raw-moment errors meet case tolerances |
| `verification_cqmom_conditional_moments_analytical` | Conditional raw-moment recovery | Parent nodes/weights and per-parent conditional sequences are explicitly composed into mixed moments | Conditional orders 0--3 and residual norms recover to `1e-13`/`1e-14` |
| `verification_cqmom_moment_reconstruction_selected` | Selected raw-moment reconstruction from CQMOM quadrature | An independent-product configuration supplies analytical tensor moments and selected multi-index requests | Selected values, request order/duplicates, terminal quadrature, and reconstruction agree with the known product rule |
| `verification_cqmom_public_workflows` | Dense, lean, and default `CqmomCase` workflows | An independently specified 2D product quadrature supplies moments, terminal nodes, and required lean layout | Each public route recovers the reference quadrature/moments and reports the expected layout and diagnostics |
| `verification_cqmom_independent_product` | Fixed 2D independent-product CQMOM | Prescribed marginal rules define the product terminal quadrature and raw tensor | Parent/conditional nodes, local/global weights, paths, diagnostics, and reconstructed moments match the rule |
| `verification_cqmom_deterministic_affine` | Fixed 2D affine conditional CQMOM | A prescribed x rule and deterministic affine y(x) conditional construction define terminal nodes and moments | Fixed topology, conditional nodes/weights, paths, and reconstructed moments match the construction |
| `verification_cqmom_branch_dependent` | Adaptive 2D branch-dependent CQMOM | Prescribed conditional rules give different supported child counts on parent branches | Adaptive counts/reduction flags, ragged paths, local/global weights, nodes, and moments match the known tree |
| `verification_cqmom_affine_chain_5d` | Five-dimensional affine-chain recursive CQMOM | A terminal rule with specified affine coordinate relations is summed independently into raw moments | Configured nontrivial conditioning, terminal paths/weights/nodes, local weights, affine relations, and moments are checked |
| `verification_cqmom_independent_product_3d` | Three-dimensional product CQMOM | Independently prescribed one-dimensional rules are combined in phase-space order | Fixed terminal paths, per-level local rules, global weights, coordinates, and moments recover |
| `verification_cqmom_independent_product_4d` | Four-dimensional product CQMOM | Independently prescribed one-dimensional rules are combined in phase-space order | Fixed topology and known rule/moment recovery, not dimension count alone |
| `verification_cqmom_independent_product_5d` | Five-dimensional product CQMOM | Independently prescribed one-dimensional rules are combined in phase-space order | Fixed topology and known rule/moment recovery, including conditioning-order diagnostics |
| `verification_cqmom_independent_product_6d` | Six-dimensional product CQMOM | Independently prescribed one-dimensional rules are combined in phase-space order | Fixed topology, paths, local/global weights, coordinates, and moments match the known rule |
| `verification_cqmom_mixed_affine_independent_4d` | Four-dimensional mixed affine/independent CQMOM | A prescribed terminal construction combines affine relations with independent rules | Expected paths, branch counts, local/global weights, coordinates, affine relations, and moments are checked |
| `verification_cqmom_automatic_scaling_reference`, `verification_cqmom_automatic_scaling_rescaled` | Adaptive 3D automatic-scaling CQMOM covariance | A prescribed five-terminal ragged rule independently generates two complete configured tensors in units `(1,1,1)` and `(1e6,1e-6,1e3)` | Non-identity conditioning `{2,0,1}`, dense/lean routes, paths, active counts, physical coordinates, weights, mixed moments, and required-layout invariance |
| `verification_cqmom_automatic_scaling_feature_sensitive` | Scaling-sensitive 3D adaptive CQMOM reconstruction | A prescribed five-terminal ragged physical rule in units `(1e8,1e-8,1e4)` independently generates the complete raw tensor | Scaling on recovers paths, topology, physical nodes, weights, and moments; scaling off must lose the prescribed five-terminal topology |
| `verification_distribution_families_analytical` | Raw moments of normal, lognormal, beta, gamma, chi-square, chi, and Weibull families | Closed-form family formulas are evaluated independently of the implementation | Multiple parameters/orders, normalization, and invalid-parameter contracts; the Weibull--gamma identity is an additional consistency check |
| `verification_distribution_mixtures_analytical` | Mixture raw-moment generation | Explicit component-density-weighted sums form the oracle | Bimodal/heterogeneous/multimodal mixtures, permutation invariance, odd-moment cancellation, zero density, and failures |
| `verification_distribution_independent_products_analytical` | Multidimensional independent-product moments | Factorized marginal moments form each tensor entry independently | One- through three-dimensional products, signed odd moments, zero density, and overflow cases |
| `verification_distribution_configuration_initialization` | Configuration-driven distribution initialization | Parsed configurations are compared with analytical moments, selected requests, and independently factorized product tensors | Families, mixtures, dimension sweep, non-unit lognormal-mixture marginal, dense/lean states, and public CQMOM reconstruction |
| `verification_realizability_analytical` | Hamburger, Stieltjes, and Hausdorff realizability classification | Known normal, gamma/lognormal, beta/uniform, zero-measure, and non-realizable sequences have known support classifications | Strict interior, boundary, negative-mass/variance, localizing, and bounded-support failures with diagnostic matrix/order checks |
| `verification_vandermonde_solver_selection` | Conditional Vandermonde solver policies on an ordinary product case | Every policy reconstructs a known prescribed product quadrature and its independently generated moments | Dense, Björck--Pereyra, pivoted QR, SVD, automatic/fallback; judged against the known rule, not solver-to-solver parity alone |
| `verification_vandermonde_solver_policy_clustered` | Conditional Vandermonde policies on a clustered, branch-dependent 2D CQMOM rule | A prescribed positive six-terminal rule independently generates the `6 x 4` mixed-moment tensor; the root weighted `3 x 3` system has `kappa_2` about `2.33e5` | Every public policy recovers the fixed topology, nodes, weights, and moments; diagnostics show requested/selected policy and no fallback on this solvable case |
| `verification_runge_kutta_analytical` | Public RK23, RK45, and SSP-RK3 solution/error behavior | Constant, exponential-decay, and oscillatory ODEs have closed-form solutions | Numerical states, embedded estimates, and method/status behavior meet case tolerances |
| `verification_ssp_rk3_convergence` | SSP-RK3 global convergence | Integrates `y'=-y` to its exact exponential solution over five halved step sizes | Consecutive errors decrease; final two observed orders exceed about `2.85` |
| `verification_rk23_adaptive_convergence` | Adaptive RK23 tolerance-refinement accuracy | Integrates `y'=-y` to the exact exponential solution under progressively tighter tolerances | Error decreases, accepted steps do not decrease, and final average-step observed trends exceed `2.2` |
| `verification_rk23_fixed_step_convergence` | Fixed-step RK23 and embedded-estimate order | Four halved steps integrate `y'=-y` against the exact solution and one-step embedded difference | Solution and embedded behavior are both consistent with third order |
| `verification_rk45_adaptive_convergence` | Adaptive RK45 tolerance-refinement accuracy | Integrates `y'=-y`, `y(0)=1`, to `T=2` against `exp(-2)` at tolerances `1e-4`, `1e-5`, `1e-6`, `1e-7` | Error decreases, accepted work does not decrease, and final average-step trends exceed `4.2`; this is not a formal fixed-step order proof |
| `verification_rk45_fixed_step_convergence` | Fixed-step RK45 and embedded-estimate order | Four halved steps integrate `y'=-y` against the exact solution and one-step embedded difference | Solution and embedded behavior are both consistent with fifth order |
| `verification_node_advection_constant_velocity` | QMOM constant-velocity node transport | Closed-form translation advances each prescribed node; weights are unchanged | Nodes, weights, reconstructed moments, and topology agree with the translated rule |
| `verification_node_advection_linear_growth` | QMOM linear-growth node transport | The scalar exponential trajectory is the exact node reference | Nodes/moments recover the exponential solution; tighter adaptive controls reduce error |
| `verification_node_advection_affine` | QMOM affine node transport | Closed-form affine ODE trajectory advances prescribed nodes | Nodes, weights, and reconstructed moments agree with the exact trajectory |
| `verification_node_advection_linear_system_2d` | CQMOM coupled 2D linear node transport | Matrix-exponential trajectories advance an independently prescribed terminal rule | Terminal nodes/weights, paths, topology, and reconstructed tensor moments match the exact solution |
| `verification_moment_source_order_dependent_relaxation` | QMOM order-dependent relaxation | The source has a closed-form exponential relaxation of each raw moment | Public source evolution is compared with the exact moment vector over configured time plans |
| `verification_moment_source_common_rate_relaxation` | CQMOM common-rate relaxation | An independently specified target terminal rule defines the exact common-rate raw-moment relaxation | CQMOM moment evolution and reconstructed quadrature agree with the analytical target law |
| `verification_moment_source_linear_coupled` | QMOM coupled linear source integration | The coupled linear system has an exact solution | Four halved plans give decreasing error; final two SSP-RK3 observed orders exceed `2.8` |
| `verification_phase_space_evolution_commuting_qmom` | Commuting QMOM advection/source splitting | Scaling advection and diagonal source rates have an exact commuting moment solution | Public split state and reconstructed quadrature match the exact law; refinement checks the configured high-order source trend |
| `verification_phase_space_evolution_qmom_boundary_contribution` | QMOM boundary contribution integration | Prescribed constant source/boundary rates give exact raw-moment increments | Boundary-only and combined increments, reconstructed moments, and failure contracts are checked |
| `verification_phase_space_evolution_qmom_phase_space_evolution_failures` | QMOM evolution failure propagation | Deliberately failing velocity/RK cases have prescribed status, failed-stage, rollback, and diagnostic outcomes | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_qmom_phase_space_evolution_inversion_failures` | QMOM post-source inversion failure propagation | A deliberately non-realizable post-source state has a prescribed failure time/status/diagnostic outcome | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_noncommuting_convergence` | QMOM noncommuting operator-splitting order | An independently derived exact linear moment-system solution is compared over five halved refinements | Strang approaches second order; `advsrc` and `srcadv` approach first order and differ at finite step |
| `verification_phase_space_evolution_noncommuting_cqmom_convergence` | CQMOM noncommuting operator-splitting order | An independently derived CQMOM moment reference is compared over five halved refinements | Strang approaches second order; `advsrc`/`srcadv` approach first order while topology/moments remain valid |
| `verification_phase_space_evolution_commuting_cqmom_2d` | Commuting 2D CQMOM advection/source splitting | Prescribed terminal rule and commuting analytical rates provide exact moments/quadrature | Dense/lean workflow results, terminal topology, and reconstructed moments meet analytical tolerances |
| `verification_phase_space_evolution_cqmom_boundary_contribution` | CQMOM boundary contribution integration | Prescribed CQMOM boundary/source rates give exact moment increments | Final stored/reconstructed moments, terminal quadrature, and boundary failure contracts are checked |
| `verification_phase_space_evolution_commuting_cqmom_state_adaptive_2d` | Adaptive-state commuting CQMOM evolution | A commuting analytical rule is evolved from a known adaptive state | Active counts/topology are preserved as required and final moments/quadrature match the exact state |
| `verification_phase_space_evolution_automatic_scaling_cqmom_reference`, `verification_phase_space_evolution_automatic_scaling_cqmom_rescaled` | Automatic scaling through adaptive CQMOM Strang evolution | Constant physical node translation and the SSP-RK3 polynomial for `dM_alpha/dt=0.1M_alpha` define the manufactured result | Configured physical units, source-context coordinates, ragged paths, final weights/nodes, and full moments match the oracle |
| `verification_phase_space_evolution_cqmom_native_source` | Native CQMOM source callback evolution | A closed-form source law supplies independently calculated required-moment values | Callback context, final time, state moments, and final quadrature meet analytical SSP-RK3 tolerance |
| `verification_phase_space_evolution_cqmom_combined_native_source` | Combined CQMOM advection and native source evolution | Independently calculated split-reference moments combine the prescribed advection and source laws | Callback observes the post-advection context; final moments/quadrature and stage diagnostics meet tolerance |
| `verification_phase_space_evolution_cqmom_native_source_failures` | Native CQMOM source failure propagation | Wrong-size, non-finite, and throwing callbacks have prescribed failure statuses/diagnostics | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_cqmom_phase_space_evolution_failures` | CQMOM evolution failure propagation | Deliberately failing velocity/RK/source cases specify status, failed stage, rollback, and diagnostic expectations | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_cqmom_phase_space_evolution_inversion_failures` | CQMOM post-source inversion failure propagation | A deliberately non-realizable post-source state specifies the failure status/time/diagnostics | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_cqmom_dense_phase_space_evolution_failures` | Dense CQMOM evolution failure propagation | The dense workflow replays deliberate advection/source failures with prescribed status and state contracts | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_cqmom_dense_phase_space_evolution_inversion_failures` | Dense CQMOM inversion-failure propagation | The dense workflow replays a prescribed non-realizable post-source state | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_cqmom_lean_phase_space_evolution_failures` | Lean CQMOM evolution failure propagation | The lean workflow replays deliberate advection/source failures with prescribed status and state contracts | Workflow/regression contract; not independent numerical verification |
| `verification_phase_space_evolution_cqmom_lean_phase_space_evolution_inversion_failures` | Lean CQMOM inversion-failure propagation | The lean workflow replays a prescribed non-realizable post-source state | Workflow/regression contract; not independent numerical verification |
| `verification_output_writer_workflows` | QMOM/CQMOM output writer workflow | Parses generated output and compares headers, ordering, requested fields, and numerical columns with known quadrature/moment data | Workflow regression; it does not independently verify a numerical algorithm |
| `verification_configuration_failure_paths` | Configuration parser failure behavior | A corpus of malformed/invalid configurations is required to produce specific parse/validation failures | Parser/failure contract; not independent numerical verification |

### Supporting/Public-Header Test (1)

| Test | Source | What | Method |
|---|---|---|---|
| `public_headers` | `tests/public_headers/*.cpp` | Self-contained installed public headers | One include per translation unit; compile and link |

### External Installed-Package Consumer Check

`tests/package_consumer/` is outside the root CTest count. After installing
QMOMKit, configure and build it as an external project to verify
`find_package(QMOMKit CONFIG REQUIRED)`, `qmomkit::qmomkit`, and installed
headers without repository include paths.

## Test-suite Maintenance

When adding, removing, renaming, or splitting a root CTest registration, update
this inventory in the same change and give every entry a distinct contract and
method. Re-run `ctest -N` and `ctest --show-only=json-v1` to record the actual
count. Keep functional validation in Debug because unit assertions are disabled
by `NDEBUG`; do not replace the assertion mechanism without an explicit design
request.
