# Monte Carlo companion package

QMOMKit includes the optional `QMOMKitMonteCarlo` companion package under
`tools/monte_carlo/`. It provides finite ensemble sampling, sample-derived raw
moments, advection sample evolution, and moment-history comparison utilities.
It is useful when you want a sampled reference beside a QMOM or CQMOM
calculation.

The companion package is **not** needed for ordinary QMOM or CQMOM
reconstruction.

## When to use it

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

For direct comparisons, initialize both the moment-based solution and the Monte
Carlo solution from the same ensemble. That eliminates initial sampling fluctuations.

## Inject a local population birth

`inject_population(...)` adds a generic birth contribution to an existing
Monte Carlo ensemble. The caller supplies the integrated birth amount
$\Delta n$ for its chosen interval; the companion package does not own physical
rate laws, time stepping, or operator splitting.

For an ensemble of $N$ samples and positive number density $n_0$, every sample
has representative weight $w = n_0 / N$. The operation draws
$K \sim \mathrm{Poisson}(\Delta n / w)$ and invokes the coordinate callback once per born
sample. That callback represents the normalized birth distribution and may
ignore its supplied random engine for deterministic injection. The returned
number density is $w (N + K)$, so existing and new samples retain the same
representative weight.

The operation owns one locally seeded `std::mt19937_64` for the Poisson draw
and all coordinate samples. Its reproducibility limits are the same as the
other sampling APIs. A zero birth amount is a successful no-op. Positive birth
requires a non-empty ensemble with positive finite number density; birth into
an empty or zero-density population is not supported by this initial API.

## Install QMOMKit core

Build and install the core QMOMKit package from the repository root first:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build
```

This writes the core package to `install/` in the repository root.

## Build the Monte Carlo package

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

## Using the C++ API

The optional package can also be used directly as an installed C++ library,
without using the supplied drivers. It provides reusable finite-ensemble
infrastructure; an application supplies its own physical event models. For
example, generic breakup, collision, recoil, and other model-specific rules in
the Monte Carlo companion examples are example-local code built on top of this
package, not library features.

Consume the separately installed package with its exported target. It brings
the required core QMOMKit package dependency with it:

```cmake
find_package(QMOMKitMonteCarlo CONFIG REQUIRED)
target_link_libraries(my_application
    PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)
```

Public headers are included as `#include <monte_carlo/...>` and their types
live in `qmomkit::monte_carlo`. Build and installation steps remain the same as
in the preceding section and in [installation.md](installation.md).

| Task | Public header / API family | Purpose |
|---|---|---|
| Represent a finite population | `ensemble.hpp`: `Sample`, `Ensemble`, `SamplingOptions` | Store generic phase-space samples in QMOMKit phase-space order. |
| Sample analytical or discrete data | `sampling.hpp`, `explicit_distribution.hpp` | Sample supported distributions, mixtures, independent products, or an explicit weighted discrete distribution. |
| Initialize from one finite realization | `sampled_initialization.hpp`: `sample_initialization(...)` | Pair an ensemble with the raw moments estimated from that same sample. |
| Estimate raw moments | `moments.hpp`: `compute_ensemble_moments(...)`, `compute_requested_ensemble_moments(...)` | Estimate complete one-dimensional/tensor moments or only requested multi-indices, with standard errors. |
| Advect samples or explicit quadrature nodes | `evolution.hpp`: `evolve_ensemble(...)`, `evolve_quadrature_nodes(...)`, `evolve_explicit_quadrature(...)` | Advance coordinates through a caller-supplied `PhaseSpaceVelocity`. Each sample is evaluated as a single packed node, and every velocity entry must be written. |
| Add a stochastic birth | `birth.hpp`: `inject_population(...)` | Add samples for a caller-supplied integrated birth amount and normalized birth sampler. |
| Write or compare estimates | `output.hpp`, `comparison.hpp` | Write moment estimates and compare QMOM/CQMOM raw moments with Monte Carlo estimates. |

The public-header comments are the detailed API reference. Most operations
return a result object with `ok()`, `status`, and `message`; use its payload
only after `ok()` succeeds. `sample_initialization(...)` follows the same
pattern and exposes both its `sampling` and `moments` results.

### Sample and initialize from one realization

This small example samples a log-normal population and computes moments from
that exact realization. Those moments can initialize a corresponding QMOM
calculation, while `initialization.sampling.ensemble` remains available as the
stochastic companion population.

```cpp
#include <monte_carlo/sampled_initialization.hpp>
#include <qmomkit/distributions.hpp>

#include <iostream>

int main()
{
    const auto initialization = qmomkit::monte_carlo::sample_initialization(
        qmomkit::LogNormalDistribution{-0.2, 0.4},
        2.5,
        3,
        qmomkit::monte_carlo::SamplingOptions{10000, 1201});

    if (!initialization.ok())
    {
        std::cerr << "Sampling: " << initialization.sampling.message
                  << "\nMoment estimation: " << initialization.moments.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::Ensemble& ensemble = initialization.sampling.ensemble;
    const qmomkit::MomentSet& moments = initialization.moments.moments;
    std::cout << "samples=" << ensemble.size() << ", m0=" << moments[0] << '\n';
    return 0;
}
```

For an existing ensemble, `compute_ensemble_moments(...)` estimates a complete
one-dimensional `MomentSet` or multidimensional `RawMomentTensor` and their
standard errors. `compute_requested_ensemble_moments(...)` evaluates only the
multi-indices listed in a `qmomkit::RawMomentRequest`, in the order requested,
which is the cheaper choice when you need a handful of moments from a large
phase space. Either way, these raw moments are what you compare against QMOM or
CQMOM results.

### Evolve an ensemble

`evolve_ensemble(...)` independently advances every sample with the supplied
phase-space velocity. It performs coordinate advection only: moment sources,
boundary rules, and event physics remain application code.

```cpp
#include <monte_carlo/evolution.hpp>

#include <iostream>
#include <vector>

int main()
{
    const qmomkit::monte_carlo::Ensemble ensemble{{0.5}, {1.0}, {2.0}};
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity(
            [](double, const std::vector<double>& coordinate)
            {
                return std::vector<double>{-0.1 * coordinate[0]};
            });

    const auto evolved =
        qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, 0.0, 0.1);
    if (!evolved.ok())
    {
        std::cerr << evolved.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::Ensemble& next_ensemble = evolved.ensemble;
    std::cout << "first evolved coordinate=" << next_ensemble.front()[0] << '\n';
    return 0;
}
```

The companion examples show how reusable sampling, moment estimation, and
advection utilities can be combined with application-specific particle models.
They are stochastic particle-based companion/reference realizations, not exact
truth solvers, so their later moments need not match deterministic QMOM/CQMOM
results exactly. See the Monte Carlo entries in the [examples
guide](examples.md) for the current 1D and multidimensional workflows.

## Run the example driver

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

## Compare moment histories

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

## Estimate selected multidimensional moments

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
duplicate indices. The estimator normalizes by total number density exactly as the dense form
does, so the two agree on any moment both compute. Check `estimate.ok()` and
read `estimate.message` when it fails.

## Output format

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

## Run tests separately

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

## Use the installed companion package

A separate project that uses the companion package should point
`CMAKE_PREFIX_PATH` at the shared install prefix:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="/path/to/QMOMKit/install"
```
