# Architecture

QMOMKit is a serial C++17 library that reconstructs quadratures from raw
moments of a number density and advances them in phase space. The core library
owns the numerics: QMOM and CQMOM reconstruction, moments of analytical
distributions, realizability checks, node advection, moment-source integration,
operator splitting, the `.qmom` parser, and output writers. Physics enters only
through callbacks. Everything else belongs to the host application, as the
[user guide](../user_new/README.md) describes.

## Package

The root CMake project builds and installs one target, `qmomkit::qmomkit`,
together with the public headers, `QMOMKitConfig.cmake`, its version file
(`SameMinorVersion` compatibility), and the exported targets. The package has no
dependencies.

Everything else is a standalone consumer, and the core never depends on one:

| Consumer | Links against |
| --- | --- |
| `tests/` (unit, component regression, verification, public headers) | The build-tree target `qmomkit`, added by the root project with `QMOMKIT_BUILD_TESTS=ON`. The tests also have `src/` on their include path. |
| `tests/package_consumer/` | The installed package, through `find_package(QMOMKit CONFIG REQUIRED)`. |
| `tests/benchmarks/` | The installed package and Google Benchmark. |
| `examples/*` | The installed package, or `qmomkit::qmomkit` when the enclosing build already defines it; the Monte Carlo examples find `QMOMKitMonteCarlo` instead. |
| `tools/monte_carlo/` | The installed package; it exports `qmomkit_monte_carlo::qmomkit_monte_carlo`, which links `qmomkit::qmomkit` publicly. |

Configuring `tests/` on its own stops with an error; it is enabled from the
root with `QMOMKIT_BUILD_TESTS=ON`.

## Data flow

A configured workflow passes through the library in this order:

1. A case holds moments, node counts and reconstruction options. Its entry
   points validate them once.
2. Reconstruction turns the moments into a quadrature: Wheeler inversion for
   QMOM ([QMOM](3_qmom.md)), and for CQMOM a recursion that alternates Wheeler
   inversions with conditional Vandermonde solves ([CQMOM](4_cqmom.md),
   [Linear algebra](5_linear_algebra.md)).
3. A state stores the accepted quadrature and the moments summed from it.
4. A step moves the nodes without reconstructing, integrates the moments with
   the source, and reconstructs once after it; a source may also request a
   quadrature of its stage moments
   ([Phase-space evolution](8_phase_space_evolution.md)).

Realizability checks ([Realizability](6_realizability.md)), analytical moments
([Distributions](7_distributions.md)) and the `.qmom` parser and writers
([Configuration and output](9_configuration_and_output.md)) sit beside this
path.

## Source map

| File | Owns | Implements |
| --- | --- | --- |
| [`src/public_api.cpp`](../../src/public_api.cpp) | `compute_qmom`, `compute_cqmom`, cases, states, state candidates, `CqmomQuadrature`, conversion of internal recursive results into `CqmomResult` | `qmom.hpp`, `cqmom.hpp` |
| [`src/wheeler.cpp`](../../src/wheeler.cpp) | Wheeler inversion, adaptive reduction, bounded-support selection | `wheeler.hpp`, `internal/bounded_quadrature.hpp` |
| [`src/recursive_quadrature.cpp`](../../src/recursive_quadrature.cpp) | Recursive CQMOM traversal, prepared branch structure, required-moment layout, `extract_cqmom_moments`, multi-index flattening, `PhaseSpaceHierarchy` | `internal/recursive_quadrature.hpp`, `multi_index.hpp`, `phase_space_hierarchy.hpp`, the layout functions of `cqmom.hpp` |
| [`src/conditioning_path.cpp`](../../src/conditioning_path.cpp) | Conditioning paths, path tables, path-to-node-index conversion | `conditioning_path.hpp` |
| [`src/coordinate_scales.cpp`](../../src/coordinate_scales.cpp) | Coordinate scales for adaptive conditional inversions | `internal/coordinate_scales.hpp` |
| [`src/moment_set.cpp`](../../src/moment_set.cpp) | `MomentSet` | `moment_set.hpp` |
| [`src/moment_reconstruction.cpp`](../../src/moment_reconstruction.cpp) | Raw and generalized moments from a quadrature | `moment_reconstruction.hpp`, `internal/moment_reconstruction.hpp` |
| [`src/realizability.cpp`](../../src/realizability.cpp) | Realizability classification, policy validation | `realizability.hpp`, `internal/realizability.hpp` |
| [`src/distributions.cpp`](../../src/distributions.cpp) | Moments of distributions, mixtures and independent products | `distributions.hpp` |
| [`src/runge_kutta.cpp`](../../src/runge_kutta.cpp) | SSP-RK3 step, embedded Runge-Kutta driver | `runge_kutta.hpp` |
| [`src/node_advection.cpp`](../../src/node_advection.cpp) | Node advection, `make_phase_space_velocity` | `node_advection.hpp` |
| [`src/moment_source.cpp`](../../src/moment_source.cpp) | Source contexts, moment-source update | `moment_source.hpp`, `internal/moment_source.hpp` |
| [`src/phase_space_evolution.cpp`](../../src/phase_space_evolution.cpp) | Operator-split steps, advection and source stages, node-count transitions, batches | `phase_space_evolution.hpp` |
| [`src/configuration.cpp`](../../src/configuration.cpp) | `.qmom` parser, loaders, `config::with_initial_moments` | `configuration.hpp` |
| [`src/output.cpp`](../../src/output.cpp) | `MomentWriter`, `QuadratureWriter` | `output.hpp` |
| [`src/internal/linear_algebra.cpp`](../../src/internal/linear_algebra.cpp) | `DenseMatrix`, dense solve with partial pivoting, symmetric Jacobi eigensolver, condition estimate | `internal/linear_algebra.hpp` |
| [`src/internal/vandermonde.cpp`](../../src/internal/vandermonde.cpp) | Conditional Vandermonde systems and solver selection | `internal/vandermonde.hpp` |
| [`src/internal/bjorck_pereyra.cpp`](../../src/internal/bjorck_pereyra.cpp), [`householder_qr.cpp`](../../src/internal/householder_qr.cpp), [`jacobi_svd.cpp`](../../src/internal/jacobi_svd.cpp) | Björck–Pereyra, pivoted QR, and SVD kernels | their headers |
| [`src/internal/embedded_runge_kutta.cpp`](../../src/internal/embedded_runge_kutta.cpp) | Embedded tableaus and step attempts | `internal/embedded_runge_kutta.hpp` |

Header-only internals: `cqmom_case.hpp` (case successors, marginal extraction,
construction access), `cqmom_validation.hpp` (layout and storage checks),
`qmom_state_access.hpp`, `numeric_utils.hpp`, `checked_size.hpp`,
`ode_rhs_validation.hpp`, `phase_space_velocity_validation.hpp`, and
`performance.hpp`. The public `realizability_policy.hpp` is header-only.

## Public and internal boundary

Installed headers live in `include/qmomkit/`; `qmomkit.hpp` includes them all.
They never include `src/internal/`. Implementation helpers live in
`src/internal/`, in namespace `qmomkit::internal`, and are not installed. Tests
may include them; that does not make them public. Configuration functions live
in `qmomkit::config`.

Public headers expose requests, options, results, statuses, moment containers,
quadratures and states. They do not expose solver workspaces, branch layouts or
parser internals. States, prepared cases and `CqmomQuadrature` have private
constructors; the implementation builds them through friend access classes,
such as `internal::QmomStateAccess` and `internal::CqmomConstructionAccess`,
which public headers only forward-declare. So only the library creates states,
and a state's moments are always those of its accepted quadrature.
`CqmomPreparedCase` and `CqmomState` hold their data behind
`std::unique_ptr<Data>`, so the prepared branch structure stays internal.

A new public header is included by `qmomkit.hpp` and gets an
`include_<header>.cpp` in `tests/public_headers/`; a change to installed
headers or exports also runs the package consumer
([Testing and validation](10_testing_and_validation.md)).
[AGENTS.md](../../AGENTS.md#public-and-internal-apis) has the full rules for
public API changes.

## Errors and statuses

An operation whose failure is a possible outcome of valid use returns a result
with a status, `ok()`, and a message: reconstruction, moment generation,
realizability checks, state candidates, every evolution operator, and
configuration loading. Internal helpers may throw; a function that promises a
status catches what they throw and maps it to a status. Callback exceptions
are mapped the same way
([Phase-space evolution](8_phase_space_evolution.md#callbacks)).

An exception means that the calling code broke a contract, that I/O failed, or
that a setup cannot start:

- malformed structure while constructing objects: `PhaseSpaceHierarchy` with
  an invalid conditioning order or a zero node count, and `CqmomPreparedCase`
  from a case whose layout is not the required layout of its hierarchy
  (`std::invalid_argument`);
- `initial_state()` when the moments of the case cannot be reconstructed
  (`std::runtime_error`);
- moment assembly and reconstruction helpers: `extract_cqmom_moments`,
  `create_cqmom_moment_data`, `reconstruct_raw_moments`,
  `reconstruct_cqmom_moments`, `evaluate_generalized_moment`, and
  `std::overflow_error` when a value leaves the range of a double;
- rectangular accessors such as `CqmomQuadrature::weight(node_index)` on ragged
  storage (`std::logic_error`);
- `ssp_rk3_step` with empty or non-finite input;
- output writers.

Data returned with a non-success status is diagnostic. A failed Wheeler
inversion can carry a one-node rule, and a failed step returns its input state;
neither is a result. The [user reference](../user_new/7_reference.md#statuses)
lists the statuses and the throwing functions a user meets.

## Thread safety

The library has no global mutable state. There are two thread-local objects:
the pool of coordinate scratch buffers in `src/node_advection.cpp`, and the
active counter session in `src/internal/performance.hpp`. Every other
workspace, such as `EmbeddedRkWorkspace`, is local to its call. Calls on
different states, cases and results may therefore run concurrently, also for
states of one prepared case, which share its branch structure read-only
through `std::shared_ptr<const ...>`. A source context caches lazily
(`mutable` members) and is valid only during the callback that receives it.
Code added to the library keeps these properties: no mutable `static` data,
and `thread_local` for scratch reused across calls.

## Coding rules

[AGENTS.md](../../AGENTS.md) holds the general rules. Three more apply to the
core:

- **Runtime checks versus assertions.** Runtime checks protect what a caller
  controls: inputs, options, shapes, non-finite values, and numerical failure.
  They map to a status or an exception by the rule above. `assert` states an
  internal invariant that earlier validation or construction guarantees, such
  as the index bounds of `DenseMatrix`. `NDEBUG` removes it, so it must never
  guard anything a caller can cause. Do not assert properties of the final
  output on intermediate values of an algorithm.
- **Validate once, at the boundary.** Public entry points validate. Internal
  `*_trusted_*` entry points, such as
  `compute_recursive_quadrature_from_trusted_moment_data` and
  `advance_trusted_qmom_moments_with_source`, skip checks that their caller
  already made, such as the layout of a prepared case. To speed up a check on
  a hot path, move it to such a boundary rather than remove it.
- **Sizes.** Size arithmetic that can overflow, such as `2 * n` moments or
  products of extents, goes through
  [`checked_size.hpp`](../../src/internal/checked_size.hpp).

Building, formatting and CI are in
[Testing and validation](10_testing_and_validation.md); installing and
embedding are in [Installation](../user_new/1_installation.md).
