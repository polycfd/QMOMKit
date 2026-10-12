# QMOMKitMonteCarlo

`QMOMKitMonteCarlo` is an optional companion package that represents a
population by a finite sample of particles. It samples distributions,
estimates raw moments with their standard errors, and moves the particles with
the same velocity callback as QMOMKit. QMOMKit itself does not need it.

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
| `sampling.hpp` | Samples of a distribution, mixture or independent product |
| `sampled_initialization.hpp` | `sample_initialization(...)`: a sample together with the moments of that sample, all of them up to a maximum order or only those of a `RawMomentRequest` |
| `moments.hpp` | Raw and selected moments of an ensemble with their standard errors, and generalized moments |
| `evolution.hpp` | `evolve_ensemble(...)`: moves every sample with a `PhaseSpaceVelocity` (RK23), calling it for one sample at a time, so a velocity that couples nodes acts on each particle alone |
| `birth.hpp` | `inject_population(...)`: adds a Poisson number of samples for a given integrated birth amount into a non-empty population, and returns the number density to use from then on |
| `output.hpp` | Moment-history files written with QMOMKit's `MomentWriter` |

Sampling, moment, evolution and birth results have `ok()`, `status` and
`message`; an initialization result has `ok()` and its `sampling` and
`moments` results. Use a payload only after `ok()`. The output helpers throw
instead. A sample repeats for the same seed, compiler and standard library.

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
for (std::size_t order = 0; order < sampled.moments.size(); ++order)
{
    std::cout << "m" << order << ": QMOM " << qmom.state.moment_values()[order] << ", Monte Carlo " << sampled.moments[order]
              << " +/- " << sampled.standard_errors[order] << '\n';
}
```

Check `ok()` on each result before using the next. A difference of a few
standard errors is sampling noise, not a discrepancy; $m_0$ has none, since
every sample carries the same weight. For CQMOM, sample an
independent product with the request
`RawMomentRequest{layout.required_indices.to_multi_indices()}` of the case's
layout; the estimated `values` then come in layout order, and
`CqmomMomentData{layout, values}` holds them.

## Examples

The Monte Carlo variants in [`examples/`](../../examples/README.md), such as
`qmom_1d_breakage_advection_with-mc`, run complete workflows: they sample the
initial population, add their own particle physics, and write QMOMKit and Monte
Carlo moment histories side by side.

## Tests

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build-tests -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_PREFIX_PATH="$PWD/install" -DBUILD_TESTING=ON
cmake --build tools/monte_carlo/build-tests
ctest --test-dir tools/monte_carlo/build-tests --output-on-failure
```

The tests check without `assert()`, so any build type runs them in full. Their
conventions and the estimator behind the standard errors are in
[Testing and validation](../../docs/developer/testing_and_validation.md#monte-carlo-references).
