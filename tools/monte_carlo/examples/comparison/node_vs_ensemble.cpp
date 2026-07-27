#include "monte_carlo/comparison.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/explicit_distribution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"
#include "qmomkit/multi_index.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

namespace
{

double integer_power(const double base, const std::size_t exponent)
{
    double result = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        result *= base;
    }
    return result;
}

std::vector<std::size_t> extents_from_maximum_orders(const std::vector<std::size_t>& maximum_orders)
{
    std::vector<std::size_t> extents(maximum_orders.size(), 0);
    for (std::size_t dimension = 0; dimension < maximum_orders.size(); ++dimension)
    {
        extents[dimension] = maximum_orders[dimension] + 1;
    }
    return extents;
}

double number_density(const qmomkit::monte_carlo::ExplicitDiscreteDistribution& distribution)
{
    double value = 0.0;
    for (const double weight : distribution.weights)
    {
        value += weight;
    }
    return value;
}

double quadrature_moment(const qmomkit::monte_carlo::ExplicitDiscreteDistribution& distribution, const qmomkit::MultiIndex& order)
{
    double value = 0.0;
    for (std::size_t node = 0; node < distribution.nodes.size(); ++node)
    {
        double monomial = 1.0;
        for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
        {
            monomial *= integer_power(distribution.nodes[node][dimension], order[dimension]);
        }
        value += distribution.weights[node] * monomial;
    }
    return value;
}

// Reconstruct reference raw moments directly from the evolved weighted nodes
// so the comparison is against the explicit quadrature, not the sampler.
qmomkit::RawMomentTensor quadrature_moment_tensor(const qmomkit::monte_carlo::ExplicitDiscreteDistribution& distribution,
                                                  const std::vector<std::size_t>& maximum_orders)
{
    const std::vector<std::size_t> extents = extents_from_maximum_orders(maximum_orders);
    const std::vector<qmomkit::MultiIndex> orders = qmomkit::enumerate_multi_indices(extents);
    std::vector<double> values(orders.size(), 0.0);
    for (std::size_t flat = 0; flat < orders.size(); ++flat)
    {
        values[flat] = quadrature_moment(distribution, orders[flat]);
    }
    return qmomkit::RawMomentTensor{extents, std::move(values)};
}

} // namespace

int main()
{
    const qmomkit::monte_carlo::ExplicitDiscreteDistribution initial_distribution{{{-0.25, 0.1}, {0.15, -0.2}, {0.35, 0.25}}, {0.9, 0.8, 1.0}, {}};
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinates)
    { return std::vector<double>{0.05 + 0.1 * coordinates[0], -0.03 * coordinates[1]}; };
    qmomkit::AdaptiveRk23Options rk23_options;
    rk23_options.relative_tolerance = 1.0e-9;
    rk23_options.absolute_tolerance = 1.0e-12;
    rk23_options.max_step = 0.05;

    constexpr double t0 = 0.0;
    constexpr double t1 = 0.2;
    const std::vector<std::size_t> maximum_orders{2, 2};
    const qmomkit::monte_carlo::SamplingOptions sampling_options{12000, 9101};
    const double density = number_density(initial_distribution);

    const qmomkit::monte_carlo::ExplicitQuadratureEvolutionResult reference =
        qmomkit::monte_carlo::evolve_explicit_quadrature(initial_distribution, velocity, t0, t1, rk23_options);
    if (!reference.ok())
    {
        std::cerr << reference.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult sampling =
        qmomkit::monte_carlo::sample_explicit_discrete_distribution(initial_distribution, sampling_options);
    if (!sampling.ok())
    {
        std::cerr << sampling.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::EnsembleEvolutionResult evolved_ensemble =
        qmomkit::monte_carlo::evolve_ensemble(sampling.ensemble, velocity, t0, t1, rk23_options);
    if (!evolved_ensemble.ok())
    {
        std::cerr << evolved_ensemble.message << '\n';
        return 1;
    }

    const qmomkit::RawMomentTensor reference_moments = quadrature_moment_tensor(reference.evolved_distribution, maximum_orders);
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(evolved_ensemble.ensemble, density, maximum_orders);
    if (!estimate.ok())
    {
        std::cerr << estimate.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::MomentComparisonResult comparison =
        qmomkit::monte_carlo::compare_moments(reference_moments, estimate.moments, estimate.standard_errors);
    if (!comparison.ok())
    {
        std::cerr << comparison.message << '\n';
        return 1;
    }

    qmomkit::monte_carlo::write_comparison_text(
        std::cout, comparison,
        qmomkit::monte_carlo::ComparisonOutputMetadata{t1, sampling_options.sample_count, sampling_options.seed, density, rk23_options});
    return 0;
}
