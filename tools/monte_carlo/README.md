# QMOMKitMonteCarlo

`QMOMKitMonteCarlo` is an optional companion package that represents a
population by a finite sample of particles. It samples distributions,
estimates raw moments with their standard errors, moves the particles with the
same velocity callback as QMOMKit, and compares moments. QMOMKit itself does not
need it.

Its purpose is a reference for QMOM and CQMOM calculations. Starting QMOMKit
from the moments of the sample, instead of the analytical ones, removes the
initial sampling difference, so what remains comes from the quadrature and its
evolution. The reference carries sampling noise that shrinks as $N^{-1/2}$ with
the sample count $N$: it is a reference, not an exact solution. Physical models,
such as breakage or collisions, stay in the application.

## Build and install

Install QMOMKit, then the companion package into the same prefix:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

A project then finds the package, which brings QMOMKit with it:

```cmake
find_package(QMOMKitMonteCarlo CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)
```

## What it provides

Headers are included as `<monte_carlo/...>`, and everything lives in
`qmomkit::monte_carlo`. The header comments are the API reference.

| Header | Provides |
| --- | --- |
| `ensemble.hpp` | `Ensemble`, a vector of samples in phase-space order, and `SamplingOptions` (sample count, seed) |
| `sampling.hpp`, `explicit_distribution.hpp` | Samples of a distribution, mixture, independent product, or explicit weighted points |
| `sampled_initialization.hpp` | `sample_initialization(...)`: a sample together with the moments of that sample, all of them up to a maximum order or only those of a `RawMomentRequest` |
| `moments.hpp` | Raw, selected and generalized moments of an ensemble, with standard errors |
| `evolution.hpp` | `evolve_ensemble(...)`: moves every sample with a `PhaseSpaceVelocity` (RK23) |
| `birth.hpp` | `inject_population(...)`: adds a Poisson number of samples for a given integrated birth amount |
| `comparison.hpp`, `output.hpp` | `compare_moments(...)` and moment-history files |

Each operation returns a result with `ok()`, `status` and `message`; use its
payload only after `ok()`. A sample repeats for the same seed, compiler and
standard library.

## Comparing QMOMKit with a sample

```cpp
namespace mc = qmomkit::monte_carlo;

// One sample of 100000 particles, and the moments of exactly that sample.
const mc::OneDimensionalSampledInitialization initial =
    mc::sample_initialization(qmomkit::LogNormalDistribution{-0.2, 0.4}, 2.5, 3, mc::SamplingOptions{100000, 1201});

// The same velocity moves the QMOM nodes and the particles: dx/dt = -0.5 x.
const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
    [](double, const std::vector<double>& x) { return std::vector<double>{-0.5 * x[0]}; });
const qmomkit::QmomMomentSource no_source = [](double, const qmomkit::QmomSourceContext& context)
{ return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

const qmomkit::QmomState state = qmomkit::QmomCase{initial.moments.moments, 2}.initial_state();
const qmomkit::QmomPhaseSpaceEvolutionResult qmom = qmomkit::advance_qmom_phase_space(state, 0.0, 1.0, velocity, no_source);
const mc::EnsembleEvolutionResult particles = mc::evolve_ensemble(initial.sampling.ensemble, velocity, 0.0, 1.0);

const mc::OneDimensionalMomentEstimate sampled =
    mc::compute_ensemble_moments(particles.ensemble, initial.sampling.number_density, 3);
const mc::MomentComparisonResult comparison =
    mc::compare_moments(qmomkit::MomentSet(qmom.state.moment_values()), sampled.moments, sampled.standard_errors);
```

Check `ok()` on each result before using the next. `compare_moments(reference,
monte_carlo, standard_errors)` returns one row per moment:

| Column | Definition |
| --- | --- |
| `signed_error` | Monte Carlo − reference |
| `absolute_error` | $\lvert$`signed_error`$\rvert$ |
| `relative_error` | `absolute_error` / max($\lvert$reference$\rvert$, standard error, $10^{-300}$) |
| `z_score` | `signed_error` / standard error; ±∞ for a nonzero error with a zero standard error |

A large z-score is a diagnostic, not a verdict. For CQMOM, sample an
independent product with the request
`RawMomentRequest{layout.required_indices.to_multi_indices()}` of the case's
layout; the estimated `values` then come in layout order, and
`CqmomMomentData{layout, values}` holds them.

## Example driver and comparison script

The driver samples the distribution of a `.qmom` file once and writes the
QMOMKit and Monte Carlo moment histories on one time grid, plus the Monte Carlo
statistics:

```bash
./tools/monte_carlo/build/monte_carlo_sample_based_initialization_example \
    tools/monte_carlo/configuration/qmom_distribution_lognormal.qmom \
    tools/monte_carlo/build/sample_based_qmom
```

For CQMOM, use `configuration/cqmom_independent_product_2d.qmom`. The script
compares two such histories, the reference first and the Monte Carlo history
second, as `compare_moments(...)` does:

```bash
python3 tools/monte_carlo/scripts/compare_moments.py \
    tools/monte_carlo/build/sample_based_qmom/sample_based_qmom_moments.txt \
    tools/monte_carlo/build/sample_based_qmom/sample_based_mc_moments.txt \
    --mc-statistics tools/monte_carlo/build/sample_based_qmom/sample_based_mc_statistics.txt \
    --output-dir tools/monte_carlo/build/sample_based_qmom/comparison
```

It writes `moment_comparison.txt`, one row per time and moment, with the
columns above. The standard errors come from `--mc-statistics`; without it,
there is no `z_score`, and `relative_error` divides by $\lvert$reference$\rvert$
alone. `--epsilon` sets the floor, and `--plots` adds PDF plots when matplotlib
is available.

## Tests

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build-tests -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_PREFIX_PATH="$PWD/install" -DBUILD_TESTING=ON
cmake --build tools/monte_carlo/build-tests
ctest --test-dir tools/monte_carlo/build-tests --output-on-failure
```

The tests check without `assert()`, so any build type runs them in full. Their
conventions and the estimator behind the standard errors are in
[Testing and validation](../../docs/developer_new/testing_and_validation.md#monte-carlo-references).
