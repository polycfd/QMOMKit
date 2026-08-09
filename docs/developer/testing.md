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

CTest names identify the tested contract. Some executables accept a small local
selector so independently meaningful groups can be registered separately;
running such an executable with no selector remains a developer convenience and
is not registered as an aggregate test. The current root suite has **103**
registered CTest tests: 32 unit, 5 component-regression, 65 verification, and
1 supporting public-header test. Obtain the current count from `ctest -N` (and
confirm structured registration with `ctest --show-only=json-v1`), rather than
maintaining it by arithmetic alone.

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
at the documented count of 103 tests.

### Unit Tests (32)

| Test | Source | What | Method |
|---|---|---|---|
| `unit_wheeler` | `tests/unit/wheeler.cpp` | Wheeler inversion and failures | Debug `assert()` contracts |
| `unit_runge_kutta` | `tests/unit/runge_kutta.cpp` | RK step APIs and option failures | Debug `assert()` contracts |
| `unit_recursive_indexing_layout` | `tests/unit/recursive_quadrature.cpp` | Recursive paths, multi-indices, and layouts | Debug `assert()` contracts |
| `unit_recursive_lean_reconstruction` | `tests/unit/recursive_quadrature.cpp` | Lean moments and dense/lean reconstruction parity | Debug `assert()` contracts |
| `unit_recursive_bounded_support` | `tests/unit/recursive_quadrature.cpp` | Bounded-support recursive CQMOM behavior | Debug `assert()` contracts |
| `unit_recursive_solver_policy` | `tests/unit/recursive_quadrature.cpp` | Recursive solver selection and diagnostics | Debug `assert()` contracts |
| `unit_conditional_moments` | `tests/unit/conditional_moments.cpp` | Conditional-moment calculations and validation | Debug `assert()` contracts |
| `unit_moment_rates` | `tests/unit/moment_rates.cpp` | Generic moment-rate evaluation | Debug `assert()` contracts |
| `unit_advection_moment_rates` | `tests/unit/advection_moment_rates.cpp` | Advection moment-rate construction | Debug `assert()` contracts |
| `unit_moment_source_rates` | `tests/unit/moment_source_rates.cpp` | Moment-source rate construction | Debug `assert()` contracts |
| `unit_bjorck_pereyra` | `tests/unit/bjorck_pereyra.cpp` | Björck-Pereyra linear solve | Debug `assert()` contracts |
| `unit_householder_qr` | `tests/unit/householder_qr.cpp` | Householder QR linear solve | Debug `assert()` contracts |
| `unit_jacobi_svd` | `tests/unit/jacobi_svd.cpp` | Jacobi SVD linear solve | Debug `assert()` contracts |
| `unit_vandermonde` | `tests/unit/vandermonde.cpp` | Vandermonde solver behavior | Debug `assert()` contracts |
| `unit_public_api` | `tests/unit/public_api.cpp` | Public QMOM/CQMOM workflow contracts | Debug `assert()` contracts |
| `unit_configuration_initialization` | `tests/unit/configuration.cpp` | Parsed distribution and initial-moment construction | Debug `assert()` contracts |
| `unit_configuration_cqmom_case_state` | `tests/unit/configuration.cpp` | CQMOM case, state, and ownership behavior | Debug `assert()` contracts |
| `unit_configuration_runtime_options` | `tests/unit/configuration.cpp` | Runtime solver, RK, evolution, and output options | Debug `assert()` contracts |
| `unit_configuration_failure_paths` | `tests/unit/configuration.cpp` | Parser and file failure behavior | Debug `assert()` contracts |
| `unit_output` | `tests/unit/output.cpp` | Output requests and quadrature serialization | Debug `assert()` contracts |
| `unit_distribution_selected_results` | `tests/unit/distributions.cpp` | Dense/selected raw-moment requests and errors | Debug `assert()` contracts |
| `unit_distribution_mixtures` | `tests/unit/distributions.cpp` | Mixture masses, selected moments, and failures | Debug `assert()` contracts |
| `unit_distribution_independent_products` | `tests/unit/distributions.cpp` | Independent products, marginal mixtures, and failures | Debug `assert()` contracts |
| `unit_realizability` | `tests/unit/realizability.cpp` | Moment realizability diagnostics | Debug `assert()` contracts |
| `unit_realizability_policy` | `tests/unit/realizability_policy.cpp` | Realizability-policy decisions | Debug `assert()` contracts |
| `unit_node_advection` | `tests/unit/node_advection.cpp` | Node-advection stages and options | Debug `assert()` contracts |
| `unit_moment_source` | `tests/unit/moment_source.cpp` | Moment-source update stages and diagnostics | Debug `assert()` contracts |
| `unit_moment_reconstruction` | `tests/unit/moment_reconstruction.cpp` | Raw-moment reconstruction APIs | Debug `assert()` contracts |
| `unit_qmom_phase_space_evolution` | `tests/unit/phase_space_evolution.cpp` | QMOM state, stages, policies, and rollback | Debug `assert()` contracts |
| `unit_cqmom_phase_space_evolution` | `tests/unit/phase_space_evolution.cpp` | CQMOM topology, stages, policies, and rollback | Debug `assert()` contracts |
| `unit_phase_space_evolution_batch` | `tests/unit/phase_space_evolution_batch.cpp` | Batched phase-space evolution | Debug `assert()` contracts |
| `unit_performance_instrumentation` | `tests/unit/performance_instrumentation.cpp` | Optional performance-instrumentation counters | Debug `assert()` contracts |

### Component-Regression Tests (5)

| Test | Source | What | Method |
|---|---|---|---|
| `component_regression_moment_source` | `tests/component_regression/moment_source.cpp` | Moment-source component interactions | Deterministic regression contracts |
| `component_regression_node_advection` | `tests/component_regression/node_advection.cpp` | Node-advection component interactions | Manufactured deterministic cases |
| `component_regression_phase_space_evolution` | `tests/component_regression/phase_space_evolution.cpp` | Operator-split evolution interactions | Deterministic regression contracts |
| `component_regression_cqmom_compute_failures` | `tests/component_regression/cqmom_compute_failures.cpp` | CQMOM failure/status propagation | Deliberate invalid and fragile cases |
| `component_regression_cqmom_case_moment_replacement` | `tests/component_regression/cqmom_case_moment_replacement.cpp` | Immutable CQMOM moment replacement | Deterministic API contracts |

### Verification Tests (65)

| Test | Source | What | Method |
|---|---|---|---|
| `verification_qmom_wheeler_analytical` | `tests/verification/qmom/wheeler_analytical.cpp` | Wheeler reconstruction | Closed-form/discrete reference |
| `verification_qmom_compute_failures` | `tests/verification/qmom/qmom_compute_failures.cpp` | QMOM public failure behavior | Manufactured failure cases |
| `verification_qmom_analytical_single_dirac` | `tests/verification/qmom/qmom_analytical.cpp` | One-node QMOM | Exact Dirac moments |
| `verification_qmom_analytical_two_node_discrete` | `tests/verification/qmom/qmom_analytical.cpp` | Two-node QMOM | Exact discrete quadrature |
| `verification_qmom_analytical_three_node_discrete` | `tests/verification/qmom/qmom_analytical.cpp` | Three-node QMOM | Exact discrete quadrature |
| `verification_qmom_analytical_uniform_n2` | `tests/verification/qmom/qmom_analytical.cpp` | Uniform two-node reconstruction | Analytical moments |
| `verification_qmom_analytical_uniform_n3` | `tests/verification/qmom/qmom_analytical.cpp` | Uniform three-node reconstruction | Analytical moments |
| `verification_qmom_analytical_exponential_n2` | `tests/verification/qmom/qmom_analytical.cpp` | Exponential two-node reconstruction | Analytical moments |
| `verification_qmom_analytical_exponential_n3` | `tests/verification/qmom/qmom_analytical.cpp` | Exponential three-node reconstruction | Analytical moments |
| `verification_qmom_analytical_gaussian_n2` | `tests/verification/qmom/qmom_analytical.cpp` | Gaussian two-node reconstruction | Analytical moments |
| `verification_qmom_analytical_gaussian_n3` | `tests/verification/qmom/qmom_analytical.cpp` | Gaussian three-node reconstruction | Analytical moments |
| `verification_cqmom_recursive_quadrature_manufactured` | `tests/verification/cqmom/recursive_quadrature.cpp` | Recursive CQMOM reconstruction | Manufactured quadrature |
| `verification_cqmom_recursive_validation` | `tests/verification/cqmom/recursive_cqmom.cpp` | Fixed and adaptive recursive quadrature | Independent terminal rules |
| `verification_cqmom_conditional_moments_analytical` | `tests/verification/cqmom/conditional_moments.cpp` | Conditional moments | Analytical formulas |
| `verification_cqmom_dense_lean_parity_3d` | `tests/verification/cqmom/dense_lean_parity_3d.cpp` | Three-dimensional dense/lean parity | Manufactured moments and reconstruction |
| `verification_cqmom_physical_transport_probe` | `tests/verification/cqmom/physical_transport_probe.cpp` | Captured CQMOM transport behavior | Deterministic regression probe |
| `verification_cqmom_moment_reconstruction_selected` | `tests/verification/cqmom/moment_reconstruction.cpp` | Selected CQMOM moment reconstruction | Independent-product reference |
| `verification_cqmom_public_workflows` | `tests/verification/cqmom/cqmom_public_workflows.cpp` | Dense, lean, and default workflows | Analytical terminal quadrature |
| `verification_cqmom_independent_product` | `tests/verification/cqmom/cqmom_2d_analytical.cpp` | Two-dimensional independent product | Analytical moments/nodes |
| `verification_cqmom_deterministic_affine` | `tests/verification/cqmom/cqmom_2d_analytical.cpp` | Two-dimensional affine conditionals | Manufactured quadrature |
| `verification_cqmom_branch_dependent` | `tests/verification/cqmom/cqmom_2d_analytical.cpp` | Branch-dependent CQMOM | Manufactured quadrature |
| `verification_cqmom_affine_chain_5d` | `tests/verification/cqmom/cqmom_high_dimensional_analytical.cpp` | Five-dimensional affine chain | Analytical recursion |
| `verification_cqmom_independent_product_3d` | `tests/verification/cqmom/cqmom_high_dimensional_analytical.cpp` | Three-dimensional independent product | Analytical moments |
| `verification_cqmom_independent_product_4d` | `tests/verification/cqmom/cqmom_high_dimensional_analytical.cpp` | Four-dimensional independent product | Analytical moments |
| `verification_cqmom_independent_product_5d` | `tests/verification/cqmom/cqmom_high_dimensional_analytical.cpp` | Five-dimensional independent product | Analytical moments |
| `verification_cqmom_independent_product_6d` | `tests/verification/cqmom/cqmom_high_dimensional_analytical.cpp` | Six-dimensional independent product | Analytical moments |
| `verification_cqmom_mixed_affine_independent_4d` | `tests/verification/cqmom/cqmom_high_dimensional_analytical.cpp` | Mixed four-dimensional CQMOM | Manufactured analytical case |
| `verification_distribution_families_analytical` | `tests/verification/distributions/distributions_analytical.cpp` | Normal, lognormal, beta, gamma, chi, and Weibull moments | Closed-form raw moments |
| `verification_distribution_mixtures_analytical` | `tests/verification/distributions/distributions_analytical.cpp` | Mixture weights, permutations, and failures | Independent weighted sums |
| `verification_distribution_independent_products_analytical` | `tests/verification/distributions/distributions_analytical.cpp` | Direct independent-product moments | Factorized analytical moments |
| `verification_distribution_configuration_initialization` | `tests/verification/distributions/distributions_analytical.cpp` | Configuration initialization, including a normal × non-unit lognormal-mixture marginal | Analytical factorization, selected moments, and public CQMOM reconstruction |
| `verification_realizability_analytical` | `tests/verification/realizability/realizability_analytical.cpp` | Realizability conditions | Analytical boundary and interior cases |
| `verification_vandermonde_solver_selection` | `tests/verification/vandermonde/solver_selection.cpp` | Solver-policy selection | Independent quadrature/reconstruction |
| `verification_runge_kutta_analytical` | `tests/verification/runge_kutta/runge_kutta_analytical.cpp` | RK solution/error behavior | Closed-form ODEs |
| `verification_ssp_rk3_convergence` | `tests/verification/runge_kutta/runge_kutta_convergence.cpp` | SSP-RK3 global order | Exponential-decay refinement |
| `verification_rk23_adaptive_convergence` | `tests/verification/runge_kutta/runge_kutta_convergence.cpp` | Adaptive RK23 tolerance trend | Exponential decay with accepted-step study |
| `verification_rk23_fixed_step_convergence` | `tests/verification/runge_kutta/runge_kutta_convergence.cpp` | Fixed-step RK23 order | Exponential-decay refinement |
| `verification_rk45_fixed_step_convergence` | `tests/verification/runge_kutta/runge_kutta_convergence.cpp` | Fixed-step RK45 order | Exponential-decay refinement |
| `verification_node_advection_constant_velocity` | `tests/verification/node_advection/node_advection_analytical.cpp` | Constant node velocity | Closed-form translation |
| `verification_node_advection_linear_growth` | `tests/verification/node_advection/node_advection_analytical.cpp` | Linear node growth | Closed-form scalar ODE |
| `verification_node_advection_affine` | `tests/verification/node_advection/node_advection_analytical.cpp` | Affine node evolution | Closed-form affine ODE |
| `verification_node_advection_linear_system_2d` | `tests/verification/node_advection/node_advection_analytical.cpp` | Two-dimensional linear node system | Matrix-exponential reference |
| `verification_moment_source_order_dependent_relaxation` | `tests/verification/moment_source/moment_source_analytical.cpp` | QMOM order-dependent relaxation | Closed-form moment source |
| `verification_moment_source_common_rate_relaxation` | `tests/verification/moment_source/moment_source_analytical.cpp` | CQMOM common-rate relaxation | Closed-form moment source |
| `verification_moment_source_linear_coupled` | `tests/verification/moment_source/moment_source_analytical.cpp` | Coupled QMOM source wrapper | Closed-form solution, four refinements, and SSP-RK3 third-order check |
| `verification_phase_space_evolution_commuting_qmom` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Commuting QMOM operators | Exact split solution |
| `verification_phase_space_evolution_qmom_boundary_contribution` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | QMOM boundary contribution | Manufactured source term |
| `verification_phase_space_evolution_qmom_phase_space_evolution_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | QMOM evolution failures | Deliberate failure propagation |
| `verification_phase_space_evolution_qmom_phase_space_evolution_inversion_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | QMOM inversion failures | Deliberate failure propagation |
| `verification_phase_space_evolution_noncommuting_convergence` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Noncommuting QMOM splitting order | Temporal refinement |
| `verification_phase_space_evolution_noncommuting_cqmom_convergence` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Noncommuting CQMOM splitting order | Temporal refinement |
| `verification_phase_space_evolution_commuting_cqmom_2d` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Commuting two-dimensional CQMOM | Exact split solution |
| `verification_phase_space_evolution_cqmom_boundary_contribution` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | CQMOM boundary contribution | Manufactured source term |
| `verification_phase_space_evolution_commuting_cqmom_state_adaptive_2d` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Adaptive CQMOM state evolution | Exact commuting case |
| `verification_phase_space_evolution_cqmom_native_source` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Native CQMOM source callback | Closed-form source update |
| `verification_phase_space_evolution_cqmom_combined_native_source` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Combined native source callback | Closed-form source update |
| `verification_phase_space_evolution_cqmom_native_source_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Native CQMOM source failures | Deliberate failure propagation |
| `verification_phase_space_evolution_cqmom_phase_space_evolution_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | CQMOM evolution failures | Deliberate failure propagation |
| `verification_phase_space_evolution_cqmom_phase_space_evolution_inversion_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | CQMOM inversion failures | Deliberate failure propagation |
| `verification_phase_space_evolution_cqmom_dense_phase_space_evolution_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Dense CQMOM evolution failures | Deliberate failure propagation |
| `verification_phase_space_evolution_cqmom_dense_phase_space_evolution_inversion_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Dense CQMOM inversion failures | Deliberate failure propagation |
| `verification_phase_space_evolution_cqmom_lean_phase_space_evolution_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Lean CQMOM evolution failures | Deliberate failure propagation |
| `verification_phase_space_evolution_cqmom_lean_phase_space_evolution_inversion_failures` | `tests/verification/phase_space_evolution/phase_space_evolution_analytical.cpp` | Lean CQMOM inversion failures | Deliberate failure propagation |
| `verification_output_writer_workflows` | `tests/verification/output/output_workflows.cpp` | QMOM and CQMOM output workflows | Parsed output against known quadrature |
| `verification_configuration_failure_paths` | `tests/verification/configuration/configuration_failures.cpp` | Invalid configuration files | Explicit parser failure corpus |

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
