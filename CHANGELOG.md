# QMOMKit Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
