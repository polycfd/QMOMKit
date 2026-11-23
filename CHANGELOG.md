# QMOMKit Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2026-11-23

### Added

- HyQMOM (Fox & Laurent, 2022) as the `hyperbolic` QMOM family, selected
  with `QmomCase::family` or `[phase_space] qmom_family`. With node count
  N, it reconstructs up to N + 1 nodes from 2N + 1 moments.
  `compute_hyqmom_characteristic_speeds` returns the characteristic speeds of
  a HyQMOM state, the weight of each speed for kinetic fluxes, and the
  moment of order 2N + 1 that closes the transport equations. Two new
  examples show the family and its speeds.
- GQMOM (Fox, Laurent & Passalacqua, 2023) as the `generalized` QMOM family,
  selected with `QmomCase::family` or `[phase_space] qmom_family`, with its
  settings in `QmomCase::generalized_options` or the `generalized_*` keys.
  With node count N, it reconstructs up to N_G > N nodes from 2N + 1
  moments. Its `normal` continuation is for a number density on the real
  line, its `gamma` and `lognormal` continuations for one on [0, inf), and
  its `beta` continuation for one on [0, 1]. For these three, abscissa bounds
  replace the interval. Three new examples show the family.
- `QmomState::family()` and `QmomState::generalized_options()`, the family
  and the settings of the generalized family of a state's case.
- `QmomResult::base_node_count` and `QmomState::base_node_count()`, the
  number of nodes of the classical quadrature on which a QMOM quadrature
  rests. For the N + 1 nodes of `hyperbolic`, it is N.
- `compute_wheeler_recurrence`, which returns the Wheeler recurrence
  coefficients, the diagonal of the Wheeler table and its rounding bounds as
  every QMOM reconstruction computes them.

### Changed

- `QmomResult` has the new last member `base_node_count`, and `QmomCase` the
  new last members `family` and `generalized_options`. Code that unpacks
  either struct with a structured binding must add names for them.
- Abscissa bounds accept a node outside them by at most 1e-12 times the
  largest node magnitude, instead of 1e-12 times the larger bound magnitude. A
  distant or very large bound therefore no longer widens the margin. Where the
  nodes are much smaller than the bounds, a node slightly outside a bound may
  now be rejected.
- An adaptive reconstruction with abscissa bounds is faster when the moments
  support fewer nodes than requested. Before, it rejected those fewer nodes
  and repeated the inversion with one node fewer each time, which returned
  the same nodes, until it asked for exactly that many. It now accepts them
  at once if they lie within the bounds. The nodes are unchanged. A CQMOM
  branch in this situation now reports the node count it started from as
  `requested_node_count`, as a branch without bounds does.
- `[output] generalized_moments` with the wrong number of exponents for the
  phase space is now reported when the file is loaded. Before, the error only
  appeared when the moment writer was created.
- The Monte Carlo companion requires a QMOMKit of at least its own version
  and of the same major version. Its installed package requires exactly the
  QMOMKit it was built against.

### Fixed

- The static libraries of QMOMKit and of the Monte Carlo companion are now
  compiled as position-independent code, so a host can link them into a
  shared library of its own on Linux.
- An optional lookup of the installed Monte Carlo companion now reports it as
  not found when QMOMKit is missing, instead of stopping with an error.

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

[Unreleased]: https://github.com/polycfd/QMOMKit/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/polycfd/QMOMKit/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/polycfd/QMOMKit/releases/tag/v1.0.0
