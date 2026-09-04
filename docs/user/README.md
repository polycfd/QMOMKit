# QMOMKit User Guide

This guide is for researchers who want to use QMOMKit as a numerical C++
library inside their own modelling, simulation, analysis, or validation
workflow.

QMOMKit helps when a model is naturally expressed through moments of a number
density function, but a later step needs a compact quadrature representation.
For example, a population-balance or reduced-order model may evolve raw
moments, then need weighted nodes for evaluating closures, evolving
phase-space coordinates, writing diagnostics, or comparing against a reference
calculation.

QMOMKit provides the reusable numerical pieces for that workflow:

- one-dimensional QMOM reconstruction from raw moment sequences;
- multidimensional CQMOM reconstruction from tensor moments or configured
  phase-space descriptions;
- analytical moment initialization for supported distributions and mixtures;
- optional realizability diagnostics for checking whether moments are
  mathematically consistent with a chosen support;
- phase-space evolution utilities that combine advection node motion with source-like
  raw-moment updates;
- output helpers and optional Monte Carlo companion tooling for sample
  realizations and validation references.

QMOMKit keeps moment ordering consistent, reconstructs quadrature rules,
initializes configured cases, reports failures, and updates phase-space
evolution state. Your application provides the coordinate meanings, units,
velocity and source callbacks, timestep policy, and application physics.

## Intended Audience

Use this guide if you want to:

- install and consume the QMOMKit CMake package;
- decide whether QMOM or CQMOM matches a modelling problem;
- construct one-dimensional QMOM quadratures;
- construct multidimensional CQMOM quadratures;
- generate raw moments from supported analytical distributions;
- configure QMOMKit examples or workflows;
- use phase-space evolution and source-evolution components;
- run examples or Monte Carlo companion tools as a user.

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

Start with installation and the quick start if you are using QMOMKit from
another project. Use the QMOM page when the population is represented by one
coordinate, and the CQMOM page when the state depends on several phase-space
coordinates. Move to distribution initialization when initial conditions are
easier to describe analytically than as moment tables, and to configuration
once you want repeatable cases that can be run by examples or drivers.
Phase-space evolution and realizability cover evolved states and support-aware checks.
Examples and Monte Carlo documentation show complete workflows and
validation-oriented usage.

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
- use `docs/developer/` when adding new numerical algorithms or changing
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

## Documentation Boundary

Use this guide for installing QMOMKit, choosing QMOM or CQMOM, writing
configuration files, running examples, and connecting phase-space evolution callbacks. For
implementation rationale, test design, validation policy, benchmark guidance,
and maintenance rules, use `docs/developer/`.
