# QMOMKit developer guide

The developer guide is for QMOMKit maintainers, contributors, advanced
developers, and AI-assisted development workflows. It documents stable
engineering knowledge needed to maintain and extend the core numerical library.

The CTest suite contains 871 individually registered logical runtime tests; see
[Testing](testing.md) for guidance and selection commands, and the [Test
inventory](test_inventory.md) for the complete logical-test reference.

This guide describes the current design of QMOMKit:

- architecture;
- numerical conventions;
- algorithms;
- testing;
- validation;
- benchmarks;
- coding practices.

It is distinct from the [user docs](../user/README.md), which explains how to
use QMOMKit from a downstream application, example, or validation workflow. It
is also distinct from source comments, which explain local implementation
details close to the code they clarify.

Developer documentation preserves stable conventions and design intent. It does
not preserve development stages, task history, roadmaps, or agent-specific
operating instructions.

## Table of contents

- [Implementation and Function Reference](function_reference/README.md): source-to-intent
  mapping and function contracts for the current implementation.
- [Developer Glossary](glossary.md): concise navigation guide for common
  QMOMKit developer symbols, functions, and concepts.
- [Architecture](architecture.md): package layout, public/internal boundaries,
  repository organization, consumer model, and design principles.
- [Numerical Conventions](numerical_conventions.md): raw NDF moment
  conventions, tensor ordering, phase-space ordering, conditioning order,
  quadrature terminology, and distribution initialization assumptions.
- [QMOM Wheeler Inversion](qmom_wheeler.md): one-dimensional QMOM inversion,
  Wheeler algorithm design, adaptive node behavior, numerical assumptions, and
  public API relationship.
- [Recursive CQMOM](cqmom_recursive.md): conditional moment representation,
  recursive construction, terminal quadrature storage, adaptive branch counts,
  failure handling, and public/internal API boundaries.
- [Internal Linear Algebra](linear_algebra.md): internal dense numerical
  kernels, the conditional Vandermonde adapter, solver policies, and
  numerical-convention ownership.
- [Distribution Design](distributions_design.md): analytical distribution
  specifications, moment generation, mixtures, independent products, numerical
  failure behavior, and extension principles.
- [Phase-space evolution Design](phase_space_evolution_design.md): node advection, moment-source
  evolution, SSP-RK3, adaptive embedded RK23/RK45, Strang splitting, coupled-direction
  schemes, status propagation, and phase-space evolution API boundaries.
- [Configuration Design](configuration_design.md): configuration philosophy,
  parsing and validation boundaries, mapping to QMOM cases, CQMOM cases,
  distribution initialization, phase-space evolution policies, output options, and extension
  principles.
- [Output Design](output_design.md): output architecture, writer
  responsibilities, moment and quadrature history contracts, format
  conventions, and configuration interaction.
- [Realizability](realizability.md): supported moment domains, realizability
  checks, diagnostics, runtime acceptance policies, phase-space evolution interaction, and
  failure classification.
- [Testing](testing.md): testing philosophy, standalone test architecture, unit
tests, verification tests, public-header tests, package-consumer tests, and
contributor expectations.
- [Test Inventory](test_inventory.md): current CTest names, categories,
  implementation sources, selectors, and claims for every logical runtime test.
- [Validation](validation.md): numerical validation strategy, analytical and
  manufactured references, Monte Carlo validation, error interpretation, and
  reproducible validation case design.
- [Monte Carlo Companion Package](monte_carlo.md): companion package boundary,
  stochastic reference construction, uncertainty, output conventions, and
  validation case guidance.
- [Benchmarks](benchmarks.md): benchmark purpose, standalone benchmark
  architecture, Google Benchmark dependency handling, benchmark categories,
  reproducibility principles, and interpretation.
- [Performance Instrumentation](performance_instrumentation.md): performance
  philosophy, internal workload counters, diagnostic separation, and profiling
  boundaries.
- [Coding Guidelines](coding_guidelines.md): QMOMKit-specific development
  principles, public API discipline, numerical code practices, documentation
  expectations, dependency philosophy, and code quality expectations.

## Guide boundaries

The developer guide describes design rules and implementation contracts that
help maintainers reason about the library. It stays focused on the core
numerical library rather than application-specific physics, CFD framework
integration, or one-off contributor procedures.
