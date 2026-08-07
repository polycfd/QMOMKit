# Monte Carlo Companion Package

QMOMKit includes the optional `QMOMKitMonteCarlo` companion package under
`tools/monte_carlo/`. It provides finite ensemble sampling, sample-derived raw
moments, advection sample evolution, and moment-history comparison utilities.
It is useful when you want a sampled reference beside a QMOM or CQMOM
calculation.

The companion package is not needed for ordinary QMOM or CQMOM reconstruction.
It is not a production particle solver, a source-term model, or a replacement
for QMOMKit phase-space evolution.

## When To Use It

Use the Monte Carlo package for two common tasks.

Sample realization:

```text
analytical distribution
    -> finite phase-space ensemble
    -> sample-derived raw moments or sample-level analysis
```

Stochastic validation reference:

```text
finite ensemble
    -> Monte Carlo sample evolution
    -> Monte Carlo raw-moment history

sample-derived raw moments
    -> QMOM or CQMOM evolution
    -> reconstructed raw-moment history
```

For direct comparisons, initialize both paths from the same finite realization
when practical. That keeps initial sampling fluctuation common to the Monte
Carlo and deterministic paths.

## Inject A Local Population Birth

`inject_population(...)` adds a generic birth contribution to an existing
Monte Carlo ensemble. The caller supplies the integrated birth amount
`Delta n` for its chosen interval; the companion package does not own physical
rate laws, time stepping, or operator splitting.

For an ensemble of `N` samples and positive number density `n0`, every sample
has representative weight `w = n0 / N`. The operation draws
`K ~ Poisson(Delta n / w)` and invokes the coordinate callback once per born
sample. That callback represents the normalized birth distribution and may
ignore its supplied random engine for deterministic injection. The returned
number density is `w * (N + K)`, so existing and new samples retain the same
representative weight.

The operation owns one locally seeded `std::mt19937_64` for the Poisson draw
and all coordinate samples. Its reproducibility limits are the same as the
other sampling APIs. A zero birth amount is a successful no-op. Positive birth
requires a non-empty ensemble with positive finite number density; birth into
an empty or zero-density population is not supported by this initial API.

## Install QMOMKit Core

Build and install the core QMOMKit package from the repository root first:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build
```

This writes the core package to `install/` in the repository root.

## Build The Monte Carlo Package

Configure the companion package against the installed QMOMKit package:

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

This builds the Monte Carlo library and the user-facing
`monte_carlo_sample_based_initialization_example` executable, then installs
the companion package into the same repository-local `install/` prefix as
QMOMKit. It does not run tests.

## Run The Example Driver

The sample-based initialization driver reads an analytical QMOM or CQMOM
configuration, samples the initial distribution once, and writes aligned
deterministic and Monte Carlo moment histories.

One-dimensional QMOM example:

```bash
./tools/monte_carlo/build/monte_carlo_sample_based_initialization_example \
    tools/monte_carlo/configuration/qmom_distribution_lognormal.qmom \
    tools/monte_carlo/build/sample_based_qmom
```

Multidimensional CQMOM example:

```bash
./tools/monte_carlo/build/monte_carlo_sample_based_initialization_example \
    tools/monte_carlo/configuration/cqmom_independent_product_2d.qmom \
    tools/monte_carlo/build/sample_based_cqmom
```

The QMOM run writes:

```text
sample_based_qmom_moments.txt
sample_based_mc_moments.txt
sample_based_mc_statistics.txt
```

The CQMOM run writes:

```text
sample_based_cqmom_moments.txt
sample_based_mc_moments.txt
sample_based_mc_statistics.txt
```

The deterministic and Monte Carlo moment-history files use the same time grid
and moment columns. The statistics file records sample count, seed, number
density, and standard-error columns.

## Compare Moment Histories

Use the comparison script on two aligned moment-history text tables:

```bash
python3 tools/monte_carlo/scripts/compare_moments.py \
    tools/monte_carlo/build/sample_based_qmom/sample_based_mc_moments.txt \
    tools/monte_carlo/build/sample_based_qmom/sample_based_qmom_moments.txt \
    --mc-statistics tools/monte_carlo/build/sample_based_qmom/sample_based_mc_statistics.txt \
    --output-dir tools/monte_carlo/build/sample_based_qmom/comparison
```

The script writes:

```text
moment_comparison.txt
```

and can write *.pdf error plots with `--plots` when matplotlib is available.

For direct QMOMKit-versus-Monte-Carlo validation, pass the Monte Carlo moment
history as the reference table and the deterministic QMOMKit history as the
computed table. Then:

```text
signed_error = computed - reference
absolute_error = abs(signed_error)
relative_error = absolute_error / (abs(reference) + epsilon)
z_score = signed_error / Monte Carlo standard error
```

Large z-scores are diagnostic signals. They are not automatic proof that either
calculation is wrong.

## Estimate Selected Multidimensional Moments

When only a requested set of multidimensional raw moments is needed, use the
selective estimator instead of constructing a dense tensor:

```cpp
#include <monte_carlo/moments.hpp>

const qmomkit::RawMomentRequest request{{{2, 0}, {0, 1}, {2, 0}}};
const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate =
    qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, request);
```

On success, `estimate.indices`, `estimate.values`, and
`estimate.standard_errors` are aligned in the original request order, including
duplicate indices. The estimator applies the same total-number-density
normalization as the dense API, computes each requested monomial directly, and
does not allocate a dense tensor for omitted moments. Check `estimate.ok()` and
use `estimate.message` for a failed validation or numerical evaluation.

## Output Format

Normal Monte Carlo output follows the same text-table convention as QMOMKit
example output:

- filenames use `.txt`;
- columns are whitespace-delimited;
- the data header is a commented line beginning with `#`;
- metadata lines are also comments;
- floating-point values are written with enough precision for reproducible
  scientific checks.

Example moment-history table:

```text
# time sample_count M0 M1 M2
0.0 100000 1.0 2.0 5.0
```

Moment columns use one-dimensional labels such as `M0`, `M1`, and `M2`, or
multidimensional labels such as `M00`, `M01`, and `M10` in QMOMKit
phase-space order. Each row in a moment-history file represents one output
time. Each row in `moment_comparison.txt` represents one time and one compared
moment.

## Run Tests Separately

Configure separate builds with tests enabled. Run both Debug and Release so the
same validation checks are exercised with and without optimization.

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build-tests-debug \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DBUILD_TESTING=ON

cmake --build tools/monte_carlo/build-tests-debug
ctest --test-dir tools/monte_carlo/build-tests-debug --output-on-failure
```

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build-tests-release \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DBUILD_TESTING=ON

cmake --build tools/monte_carlo/build-tests-release
ctest --test-dir tools/monte_carlo/build-tests-release --output-on-failure
```

Tests are not run by `cmake --build`; run `ctest` explicitly.

## Use The Installed Companion Package

A separate project that uses the companion package should point
`CMAKE_PREFIX_PATH` at the shared install prefix:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="/path/to/QMOMKit/install"
```
