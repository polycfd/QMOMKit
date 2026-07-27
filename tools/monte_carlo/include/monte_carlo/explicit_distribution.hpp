#pragma once

#include <optional>
#include <vector>

namespace qmomkit::monte_carlo
{

/**
 * Explicit finite distribution used for reference-sampling checks.
 *
 * `nodes` store phase-space coordinates in QMOMKit dimension order. `weights`
 * are absolute quadrature weights, not probabilities; sampling normalizes them
 * internally as P(node j) = w_j / sum_j w_j. If `number_density` is omitted,
 * the total number density is inferred as the sum of the absolute weights.
 */
struct ExplicitDiscreteDistribution
{
    std::vector<std::vector<double>> nodes;
    std::vector<double> weights;
    std::optional<double> number_density;
};

} // namespace qmomkit::monte_carlo
