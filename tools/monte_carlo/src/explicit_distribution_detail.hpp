#pragma once

#include "monte_carlo/explicit_distribution.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace qmomkit::monte_carlo::detail
{

struct ExplicitDistributionValidation
{
    bool ok = false;
    double number_density = 0.0;
    std::size_t dimensions = 0;
    std::string message;
};

[[nodiscard]] inline bool all_finite_values(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

[[nodiscard]] inline bool all_finite_nodes(const std::vector<std::vector<double>>& nodes)
{
    return std::all_of(nodes.begin(), nodes.end(), [](const std::vector<double>& node) { return all_finite_values(node); });
}

[[nodiscard]] inline ExplicitDistributionValidation validate_explicit_distribution(const ExplicitDiscreteDistribution& distribution)
{
    ExplicitDistributionValidation validation;
    if (distribution.nodes.empty())
    {
        validation.message = "Explicit discrete distribution must contain at least one node.";
        return validation;
    }
    if (distribution.nodes.size() != distribution.weights.size())
    {
        validation.message = "Explicit discrete distribution nodes and weights must have the same size.";
        return validation;
    }

    validation.dimensions = distribution.nodes.front().size();
    if (validation.dimensions == 0)
    {
        validation.message = "Explicit discrete distribution nodes must contain at least one coordinate.";
        return validation;
    }
    for (const std::vector<double>& node : distribution.nodes)
    {
        if (node.size() != validation.dimensions)
        {
            validation.message = "Explicit discrete distribution nodes must have consistent dimensions.";
            return validation;
        }
    }
    if (!all_finite_nodes(distribution.nodes) || !all_finite_values(distribution.weights))
    {
        validation.message = "Explicit discrete distribution contains a non-finite coordinate or weight.";
        return validation;
    }

    for (const double weight : distribution.weights)
    {
        if (weight < 0.0)
        {
            validation.message = "Explicit discrete distribution weights must be non-negative.";
            return validation;
        }
        validation.number_density += weight;
        if (!std::isfinite(validation.number_density))
        {
            validation.message = "Explicit discrete distribution weight sum became non-finite.";
            return validation;
        }
    }
    if (validation.number_density <= 0.0)
    {
        validation.message = "Explicit discrete distribution must have positive number density.";
        return validation;
    }

    if (distribution.number_density.has_value())
    {
        // An explicit density is metadata for the same absolute quadrature
        // weights, not an independent normalization.
        const double supplied = *distribution.number_density;
        if (!std::isfinite(supplied) || supplied <= 0.0)
        {
            validation.message = "Supplied number density must be finite and positive.";
            return validation;
        }
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(validation.number_density));
        if (std::abs(supplied - validation.number_density) > tolerance)
        {
            validation.message = "Supplied number density must match the sum of absolute weights.";
            return validation;
        }
        validation.number_density = supplied;
    }

    validation.ok = true;
    validation.message.clear();
    return validation;
}

} // namespace qmomkit::monte_carlo::detail
