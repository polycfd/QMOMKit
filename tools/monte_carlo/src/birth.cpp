#include "monte_carlo/birth.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <random>
#include <string>
#include <utility>

namespace qmomkit::monte_carlo
{
namespace
{

[[nodiscard]] PopulationBirthResult birth_failure(const unsigned long long seed, const SamplingStatus status, std::string message)
{
    return PopulationBirthResult{{}, 0.0, seed, status, std::move(message)};
}

[[nodiscard]] PopulationBirthResult validate_positive_birth_inputs(const Ensemble& ensemble, const double number_density,
                                                                   const BirthCoordinateSampler& sample_birth, const unsigned long long seed,
                                                                   std::size_t& dimensions, double& representative_weight)
{
    if (ensemble.empty())
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Positive population birth requires a non-empty ensemble.");
    }

    dimensions = ensemble.front().size();
    if (dimensions == 0)
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Samples must contain at least one phase-space coordinate.");
    }
    for (const Sample& sample : ensemble)
    {
        if (sample.size() != dimensions)
        {
            return birth_failure(seed, SamplingStatus::invalid_input, "All samples must have the same phase-space dimension.");
        }
        for (const double coordinate : sample)
        {
            if (!std::isfinite(coordinate))
            {
                return birth_failure(seed, SamplingStatus::non_finite_sample, "Sample coordinate is non-finite.");
            }
        }
    }

    if (!std::isfinite(number_density) || number_density <= 0.0)
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Positive population birth requires finite positive number density.");
    }
    if (!sample_birth)
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Birth coordinate sampler must be provided.");
    }

    representative_weight = number_density / static_cast<double>(ensemble.size());
    if (!std::isfinite(representative_weight) || representative_weight <= 0.0)
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Representative sample weight must be finite and positive.");
    }
    return {};
}

} // namespace

PopulationBirthResult inject_population(const Ensemble& ensemble, const double number_density, const double birth_number_density,
                                        const BirthCoordinateSampler& sample_birth, const unsigned long long seed)
{
    if (!std::isfinite(birth_number_density) || birth_number_density < 0.0)
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Birth number density must be finite and non-negative.");
    }
    if (birth_number_density == 0.0)
    {
        return PopulationBirthResult{ensemble, number_density, seed, SamplingStatus::success, ""};
    }

    std::size_t dimensions = 0;
    double representative_weight = 0.0;
    PopulationBirthResult validation =
        validate_positive_birth_inputs(ensemble, number_density, sample_birth, seed, dimensions, representative_weight);
    if (!validation.ok())
    {
        return validation;
    }

    const double mean_birth_count = birth_number_density / representative_weight;
    if (!std::isfinite(mean_birth_count))
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Birth sample-count mean must be finite.");
    }

    std::mt19937_64 engine(seed);
    std::poisson_distribution<std::size_t> birth_count_sampler(mean_birth_count);
    const std::size_t birth_count = birth_count_sampler(engine);
    if (birth_count > std::numeric_limits<std::size_t>::max() - ensemble.size())
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Birth sample count exceeds ensemble capacity.");
    }

    Ensemble augmented = ensemble;
    augmented.reserve(ensemble.size() + birth_count);
    for (std::size_t birth = 0; birth < birth_count; ++birth)
    {
        Sample sample;
        try
        {
            sample = sample_birth(engine);
        }
        catch (const std::exception& error)
        {
            return birth_failure(seed, SamplingStatus::invalid_input, std::string("Birth coordinate sampler failed: ") + error.what());
        }
        catch (...)
        {
            return birth_failure(seed, SamplingStatus::invalid_input, "Birth coordinate sampler failed with an unknown exception.");
        }

        if (sample.size() != dimensions)
        {
            return birth_failure(seed, SamplingStatus::invalid_input, "Birth coordinate sampler returned the wrong phase-space dimension.");
        }
        for (const double coordinate : sample)
        {
            if (!std::isfinite(coordinate))
            {
                return birth_failure(seed, SamplingStatus::non_finite_sample, "Birth coordinate sampler returned a non-finite coordinate.");
            }
        }
        augmented.push_back(std::move(sample));
    }

    const double updated_number_density = representative_weight * static_cast<double>(augmented.size());
    if (!std::isfinite(updated_number_density))
    {
        return birth_failure(seed, SamplingStatus::invalid_input, "Updated number density became non-finite.");
    }
    return PopulationBirthResult{std::move(augmented), updated_number_density, seed, SamplingStatus::success, ""};
}

} // namespace qmomkit::monte_carlo
