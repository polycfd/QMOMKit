# QMOMKit Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `compute_wheeler_recurrence`, which returns the Wheeler recurrence
  coefficients, the diagonal of the Wheeler table and its rounding bounds as
  every QMOM reconstruction computes them.
- HyQMOM (Fox & Laurent, 2022) as the `hyperbolic` QMOM family, selected
  with `QmomOptions::family` or `[phase_space] qmom_family`. With node count
  N, it reconstructs up to N + 1 nodes from 2N + 1 moments.
  `compute_hyqmom_characteristic_speeds` returns the characteristic speeds of
  a HyQMOM state. Two new examples show both.
- `QmomResult::base_node_count` and `QmomState::base_node_count()`, the
  number of nodes of the classical quadrature on which a QMOM quadrature
  rests. For the N + 1 nodes of `hyperbolic`, it is N.
- GQMOM (Fox, Laurent & Passalacqua, 2023) with the `normal` continuation as
  the `generalized` QMOM family, selected with `QmomOptions::family` or
  `[phase_space] qmom_family`. With node count N, it reconstructs up to N_G > N
  nodes from 2N + 1 moments. A new example shows the family.

### Changed

- Abscissa bounds accept a node outside them by at most 1e-12 times the
  largest node magnitude, instead of 1e-12 times the larger bound magnitude. A
  distant or very large bound therefore no longer widens the margin. Where the
  nodes are much smaller than the bounds, a node slightly outside a bound may
  now be rejected.

## [1.0.0] - 2026-10-28

First public release.

### Added

- Quadrature reconstruction from raw moments: QMOM with the Wheeler
  algorithm for one phase-space coordinate, and CQMOM for several, with a
  configurable conditioning order.
- Adaptive node counts, and abscissa bounds on an interval.
- Phase-space evolution by Strang or Lie splitting of node advection
  (RK23, RK45, SSP-RK3) and moment sources (SSP-RK3), including exchange
  across phase-space boundaries.
- Realizability checks of moment sequences on the real line, the half-line
  and bounded intervals, and realizability policies for evolution.
- States rebuilt from moments that the host has transported in physical
  space.
- Case setup in C++ or from `.qmom` configuration files, with moments given
  explicitly or generated from analytical distributions, mixtures and
  independent products.
- Writers for moment and quadrature histories.
- Examples, user and developer documentation, and the CMake package
  `QMOMKit` with the target `qmomkit::qmomkit`, with no third-party
  dependencies.
- The optional Monte Carlo companion package `QMOMKitMonteCarlo`, with the
  target `qmomkit_monte_carlo::qmomkit_monte_carlo`, which computes
  finite-sample reference solutions to compare QMOMKit results against.

[Unreleased]: https://github.com/polycfd/QMOMKit/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/polycfd/QMOMKit/releases/tag/v1.0.0
