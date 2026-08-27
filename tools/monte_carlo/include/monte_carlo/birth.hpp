#pragma once

#include "monte_carlo/sampling.hpp"

#include <functional>
#include <random>
#include <string>

namespace qmomkit::monte_carlo
{

/**
 * Generate one coordinate realization from a normalized birth distribution.
 *
 * The operation owns the pseudo-random engine for one injection call and passes
 * it to this callback. A deterministic injection may ignore the engine and
 * return the same phase-space coordinate every time.
 */
using BirthCoordinateSampler = std::function<Sample(std::mt19937_64&)>;

/**
 * Result of a local population birth/injection operation.
 *
 * On success, `ensemble` contains the original and born samples, and
 * `number_density` is the absolute density represented by that augmented
 * equally weighted ensemble. On failure, `ensemble` is empty and its contents
 * must not be used.
 */
struct PopulationBirthResult
{
    Ensemble ensemble;
    double number_density = 0.0;
    unsigned long long seed = 0;
    SamplingStatus status = SamplingStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == SamplingStatus::success;
    }
};

/**
 * Add a stochastic local population birth to an equally weighted ensemble.
 *
 * `birth_number_density` is the caller-supplied integrated birth amount over a
 * caller-chosen interval. For a non-empty ensemble of size `N` with positive
 * number density `n0`, every existing sample has representative weight
 * `w = n0 / N`. The operation draws `K ~ Poisson(birth_number_density / w)`,
 * appends `K` samples from `sample_birth`, and returns
 * `n0_new = w * (N + K)`. The existing representative weight is therefore
 * preserved.
 *
 * A zero birth amount succeeds without validating or sampling the ensemble.
 * Positive birth requires a non-empty, finite, dimensionally consistent
 * ensemble with positive finite number density. Positive birth into an empty
 * or zero-density population is not supported by this API.
 */
[[nodiscard]] PopulationBirthResult inject_population(const Ensemble& ensemble, double number_density, double birth_number_density,
                                                      const BirthCoordinateSampler& sample_birth, unsigned long long seed);

} // namespace qmomkit::monte_carlo
