#pragma once

#include "monte_carlo/ensemble.hpp"

#include <cmath>
#include <cstddef>

namespace qmomkit::monte_carlo::internal
{

enum class EnsembleProblem
{
    none,
    empty,
    no_coordinates,
    inconsistent_dimensions,
    non_finite_coordinate
};

struct EnsembleCheck
{
    EnsembleProblem problem = EnsembleProblem::none;
    // Coordinates of the first sample; zero for an empty ensemble.
    std::size_t dimensions = 0;
    // The sample the problem was found in; zero for an empty ensemble or one
    // without coordinates.
    std::size_t sample_index = 0;
};

/**
 * Check that an ensemble is non-empty, that its samples have coordinates, that
 * every sample has the dimension of the first, and that every coordinate is
 * finite, and report the first problem found.
 *
 * By default the samples are checked one at a time, so a sample with a
 * non-finite coordinate is reported before a later sample of the wrong
 * dimension. With `dimensions_first`, every dimension is checked before any
 * coordinate, the order `compute_requested_ensemble_moments` reports in.
 */
[[nodiscard]] inline EnsembleCheck check_ensemble(const Ensemble& ensemble, const bool dimensions_first = false)
{
    if (ensemble.empty())
    {
        return {EnsembleProblem::empty, 0, 0};
    }
    const std::size_t dimensions = ensemble.front().size();
    if (dimensions == 0)
    {
        return {EnsembleProblem::no_coordinates, 0, 0};
    }

    const auto inconsistent = [&](const std::size_t sample) { return ensemble[sample].size() != dimensions; };
    const auto non_finite = [&](const std::size_t sample)
    {
        for (const double coordinate : ensemble[sample])
        {
            if (!std::isfinite(coordinate))
            {
                return true;
            }
        }
        return false;
    };

    if (dimensions_first)
    {
        for (std::size_t sample = 0; sample < ensemble.size(); ++sample)
        {
            if (inconsistent(sample))
            {
                return {EnsembleProblem::inconsistent_dimensions, dimensions, sample};
            }
        }
    }
    for (std::size_t sample = 0; sample < ensemble.size(); ++sample)
    {
        if (inconsistent(sample))
        {
            return {EnsembleProblem::inconsistent_dimensions, dimensions, sample};
        }
        if (non_finite(sample))
        {
            return {EnsembleProblem::non_finite_coordinate, dimensions, sample};
        }
    }
    return {EnsembleProblem::none, dimensions, 0};
}

// The message every module reports for a problem; birth words the empty case
// itself.
[[nodiscard]] inline const char* ensemble_problem_message(const EnsembleProblem problem)
{
    switch (problem)
    {
    case EnsembleProblem::empty:
        return "Ensemble must contain at least one sample.";
    case EnsembleProblem::no_coordinates:
        return "Samples must contain at least one phase-space coordinate.";
    case EnsembleProblem::inconsistent_dimensions:
        return "All samples must have the same phase-space dimension.";
    case EnsembleProblem::non_finite_coordinate:
        return "Sample coordinate is non-finite.";
    case EnsembleProblem::none:
        break;
    }
    return "";
}

} // namespace qmomkit::monte_carlo::internal
