# QMOMKit User Guide

QMOMKit is a C++ library for the quadrature method of moments (QMOM) and
conditional quadrature method of moments (CQMOM). It provides the numerical
tools needed to represent and evolve a population in phase space using a finite
set of raw moments.

QMOMKit can reconstruct one- and multidimensional quadratures from moments,
evolve their abscissas under phase-space dynamics, and account for population
changes through moment sources and phase-space boundary contributions. It can
be used for standalone population models or embedded in a larger simulation
code, while the application remains responsible for the physical model and,
where applicable, transport and discretization in physical space.

A typical QMOMKit workflow starts from a set of raw moments, reconstructs the
corresponding quadrature, and then either uses that quadrature directly to
evaluate population-dependent quantities or advances the moment state through
phase-space evolution.

QMOMKit provides the reusable numerical pieces for that workflow:

- one-dimensional quadrature reconstruction from raw moment sets using the quadrature method of moments (QMOM) [[McGraw, 1997](https://www.doi.org/10.1080/02786829708965471)];
- multidimensional quadrature reconstruction from raw moment sets using the conditional quadrature method of moments (CQMOM) [[Yuan & Fox, 2011](https://www.doi.org/10.1016/j.jcp.2011.07.020)];
- analytical moment initialization for supported distributions and mixtures;
- optional realizability diagnostics for checking whether moments are
  mathematically consistent with a chosen support;
- phase-space evolution utilities that combine advection node motion with source-like
  raw-moment updates;
- output helpers and optional Monte Carlo companion tooling for sample
  realizations and validation references.

QMOMKit handles moment indexing, quadrature reconstruction, configured case
setup, diagnostics, and local phase-space evolution. The application defines
the physical meaning and units of the phase-space variables and supplies the
model-specific velocity, source, and boundary callbacks. In coupled
simulations, the application also remains responsible for physical-space
transport, spatial discretization, and the surrounding solver infrastructure.

## Recommended Reading Path

For a first pass through QMOMKit, read:

```text
installation
  -> quick_start
  -> qmom / cqmom
  -> distributions / configuration
  -> phase-space evolution / realizability
  -> examples / monte_carlo
```

Start with the installation guide and quick start if you are using QMOMKit
from another project. The QMOM guide covers one-dimensional populations,
whereas the CQMOM guide introduces multidimensional phase spaces and
conditional quadrature reconstruction.

Use the distribution-initialization guide when initial conditions are more
naturally specified through analytical distributions than explicit moment
tables, and the configuration guide when you want reproducible cases defined
through `.qmom` files. The phase-space evolution and realizability guides cover
time-dependent population states, source and boundary contributions, and
support-aware realizability checks.

The examples and Monte Carlo documentation provide complete workflows and
reference calculations for verification and comparison.

If you already have a concrete modelling workflow, these entry points are often
the most direct:

- known one-dimensional moments: [One-Dimensional QMOM](qmom.md);
- multidimensional phase space with configured marginals:
  [Multidimensional CQMOM](cqmom.md);
- analytical initial distributions:
  [Distributions](distributions.md);
- repeatable input files: [Configuration](configuration.md);
- time evolution through velocity and source models:
  [Phase-space evolution](phase_space_evolution.md);
- checking support consistency: [Realizability](realizability.md).

## Output At A Glance

QMOMKit-native output helpers write whitespace-delimited text tables.
`MomentWriter` records native integer-order raw moments. `QuadratureWriter`
records reconstructed node weights and abscissas. Generalized moments with real
exponents are derived quantities evaluated separately from a valid quadrature
result, so they are not part of the standard moment writer output.

## User Extension Points

Most user extensions happen outside the core numerical algorithms:

- use configuration files or direct C++ setup to describe phase space, initial
  moments or distributions, numerical settings, phase-space evolution scheme, and output
  cadence;
- provide physics callbacks for phase-space velocity, raw-moment source terms,
  and optional boundary contributions;
- evaluate derived quantities from reconstructed quadrature in application code;
- use the [developer documentation](../developer/README.md) when adding new numerical algorithms or changing
  QMOMKit internals.

## User Documents

- [Glossary](glossary.md): concise explanations of scientific terms used
  throughout the user guide.
- [Installation](installation.md): building, installing, and consuming
  QMOMKit.
- [Quick Start](quick_start.md): minimal QMOM and CQMOM examples.
- [One-Dimensional QMOM](qmom.md): one-dimensional QMOM inputs, outputs,
  options, and statuses.
- [Multidimensional CQMOM](cqmom.md): multidimensional CQMOM hierarchy,
  conditioning order, and quadrature output.
- [Distributions](distributions.md): analytical raw moments, mixtures, and
  independent products.
- [Configuration](configuration.md): configuration-file format and examples.
- [Phase-space evolution](phase_space_evolution.md): node advection, moment-source updates, and
  supported phase-space evolution schemes, including physics callback requirements.
- [Realizability](realizability.md): optional one-dimensional realizability
  checks and runtime policies.
- [Examples](examples.md): recommended path through example applications.
- [Monte Carlo Companion Package](monte_carlo.md): sample realizations and
  stochastic validation-reference workflows.

## Examples

For a more hands-on introduction, see the examples in [`examples/`](../../examples/). They provide complete, runnable QMOM and CQMOM cases of increasing complexity and show how configuration, quadrature reconstruction, phase-space evolution, sources, and boundary contributions fit together in practice.