# QMOMKit developer documentation

These pages are for people and AI agents who change the QMOMKit core library.
They explain what the code does not say itself: why the algorithms are built
as they are, which conventions and invariants they keep, and how a change is
tested. For everything else, the source is the reference. How to use QMOMKit
is in the [user documentation](../user_new/README.md); the rules for working on
the repository are in [AGENTS.md](../../AGENTS.md) and
[CONTRIBUTING.md](../../CONTRIBUTING.md).

## Pages

Read the first two pages before any change; they hold the structure and the
conventions every other page assumes.

1. [Architecture](1_architecture.md): the package, the source files and what
   they own, the public and internal boundary, errors and statuses, and thread
   safety.
2. [Numerical conventions](2_numerical_conventions.md): moments, indices,
   orderings, quadrature storage, tolerances, supports and bounds.

Reconstruction, from moments to a quadrature:

3. [QMOM](3_qmom.md): the Wheeler inversion, adaptive reduction and abscissa
   bounds.
4. [CQMOM](4_cqmom.md): the recursive reconstruction, its moment layout,
   branches and prepared cases.
5. [Linear algebra](5_linear_algebra.md): the internal dense kernels and the
   conditional Vandermonde solve.

Moments:

6. [Realizability](6_realizability.md): the moment checks per support, their
   numerical classification, and realizability policies.
7. [Distributions](7_distributions.md): moments of analytical distributions,
   mixtures and independent products.

Evolution, input and output:

8. [Phase-space evolution](8_phase_space_evolution.md): states, operator
   splitting, node advection, the moment-source step and node counts.
9. [Configuration and output](9_configuration_and_output.md): the `.qmom`
   parser and the output writers.

Evidence:

10. [Testing and validation](10_testing_and_validation.md): the test tree, CI,
    and the evidence numerical changes need.
11. [Benchmarks and performance](11_benchmarks_and_performance.md): the
    benchmark project and the performance counters.

## Making a change

| Changing | Code | Read first |
| --- | --- | --- |
| QMOM inversion, adaptive reduction, bounds | `src/wheeler.cpp` | [QMOM](3_qmom.md) |
| CQMOM reconstruction, layouts, branch counts | `src/recursive_quadrature.cpp`, `src/public_api.cpp` | [CQMOM](4_cqmom.md) |
| Moments and generalized moments of a quadrature | `src/moment_reconstruction.cpp` | [CQMOM](4_cqmom.md#moments-from-a-quadrature) |
| Dense solvers, Vandermonde systems | `src/internal/` | [Linear algebra](5_linear_algebra.md) |
| Realizability checks and policies | `src/realizability.cpp` | [Realizability](6_realizability.md) |
| Distribution families | `src/distributions.cpp` | [Distributions](7_distributions.md) |
| Time stepping, sources, node counts in evolution | `src/phase_space_evolution.cpp`, `src/node_advection.cpp`, `src/moment_source.cpp`, `src/runge_kutta.cpp` | [Phase-space evolution](8_phase_space_evolution.md) |
| Configuration keys, output files | `src/configuration.cpp`, `src/output.cpp` | [Configuration and output](9_configuration_and_output.md) |
| A public header | `include/qmomkit/` | [Architecture](1_architecture.md#public-and-internal-boundary) |
| A numerical convention | anywhere | [Numerical conventions](2_numerical_conventions.md) |

Every change then follows the same steps:

1. Add or update tests at the level that owns the behavior, and give a
   numerical change evidence beyond agreement with the old code
   ([Testing and validation](10_testing_and_validation.md)).
2. Build and run the suite in Debug, and format the changed files.
3. Update the page that owns the topic, and the
   [user documentation](../user_new/README.md) where public behavior or a
   configuration key changes.
4. Run the examples a change affects and check their numbers.

Tolerances, fallbacks, node ordering and status mapping change only with
explicit approval ([AGENTS.md](../../AGENTS.md#numerical-contracts)).
