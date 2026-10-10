# QMOMKit developer documentation

These pages are for contributors to the QMOMKit core library. They explain the
design, the numerical conventions and the invariants the code relies on, and
how changes are tested. How to use QMOMKit is described in the
[user documentation](../user_new/README.md); the rules for working on the
repository are in [AGENTS.md](../../AGENTS.md) and
[CONTRIBUTING.md](../../CONTRIBUTING.md).

The pages explain what the code does not say itself: why an algorithm is built
as it is, which invariants it keeps, and which page or file owns a topic. For
everything else, the source is the reference.

## Pages

1. [Architecture](1_architecture.md): the package, the source files and what
   they own, the public and internal boundary, errors and statuses, and thread
   safety.
2. [Numerical conventions](2_numerical_conventions.md): moments, indices,
   orderings, quadrature storage, tolerances, supports and bounds.
3. [QMOM](3_qmom.md): the Wheeler inversion, adaptive reduction and abscissa
   bounds.
4. [CQMOM](4_cqmom.md): the recursive reconstruction, its moment layout,
   branches and prepared cases.
5. [Linear algebra](5_linear_algebra.md): the internal dense kernels and the
   conditional Vandermonde solve.
6. [Realizability](6_realizability.md): the moment checks per support, their
   numerical classification, and realizability policies.
7. [Distributions](7_distributions.md): moments of analytical distributions,
   mixtures and independent products.
8. [Phase-space evolution](8_phase_space_evolution.md): states, operator
   splitting, node advection, the moment-source step and node counts.
9. [Configuration and output](9_configuration_and_output.md): the `.qmom`
   parser and the output writers.
10. [Testing and validation](10_testing_and_validation.md): the test tree, CI,
    and the evidence numerical changes need.
11. [Benchmarks and performance](11_benchmarks_and_performance.md): the
    benchmark project and the performance counters.

## Where to make a change

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
| A public header | `include/qmomkit/` | [Architecture](1_architecture.md), [Testing](10_testing_and_validation.md) |
| A numerical convention | anywhere | [Numerical conventions](2_numerical_conventions.md) |

A change to public behavior also updates the
[user documentation](../user_new/README.md) where it describes that behavior,
and a new configuration key appears in its
[Reference](../user_new/7_reference.md).
