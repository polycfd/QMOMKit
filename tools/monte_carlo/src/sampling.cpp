#include "monte_carlo/sampling.hpp"

#include "explicit_distribution_detail.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace qmomkit::monte_carlo
{
namespace
{

template <class... Ts> struct Overloaded : Ts...
{
    using Ts::operator()...;
};

template <class... Ts> Overloaded(Ts...) -> Overloaded<Ts...>;

[[nodiscard]] SamplingResult invalid_sampling_result(const SamplingOptions& options, std::string message)
{
    return SamplingResult{{}, 0.0, options.seed, SamplingStatus::invalid_input, std::move(message)};
}

[[nodiscard]] ExplicitDiscreteSamplingResult invalid_explicit_sampling_result(const SamplingOptions& options, std::string message)
{
    return ExplicitDiscreteSamplingResult{{}, 0.0, options.seed, ComparisonStatus::invalid_input, std::move(message)};
}

[[nodiscard]] bool is_positive_finite(const double value)
{
    return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] bool is_nonnegative_finite(const double value)
{
    return std::isfinite(value) && value >= 0.0;
}

[[nodiscard]] bool valid_distribution(const qmomkit::DistributionSpecification& distribution)
{
    return std::visit(
        Overloaded{
            [](const qmomkit::NormalDistribution& normal) { return std::isfinite(normal.mean) && is_positive_finite(normal.standard_deviation); },
            [](const qmomkit::LogNormalDistribution& lognormal)
            { return std::isfinite(lognormal.normal_mean) && is_positive_finite(lognormal.normal_standard_deviation); },
            [](const qmomkit::BetaDistribution& beta) { return is_positive_finite(beta.alpha) && is_positive_finite(beta.beta); },
            [](const qmomkit::GammaDistribution& gamma) { return is_positive_finite(gamma.shape) && is_positive_finite(gamma.scale); },
            [](const qmomkit::ChiDistribution& chi) { return is_positive_finite(chi.degrees_of_freedom); },
            [](const qmomkit::ChiSquaredDistribution& chi_squared) { return is_positive_finite(chi_squared.degrees_of_freedom); },
            [](const qmomkit::WeibullDistribution& weibull) { return is_positive_finite(weibull.shape) && is_positive_finite(weibull.scale); }},
        distribution);
}

[[nodiscard]] double sample_value(const qmomkit::DistributionSpecification& distribution, std::mt19937_64& engine)
{
    return std::visit(Overloaded{[&engine](const qmomkit::NormalDistribution& normal)
                                 {
                                     std::normal_distribution<double> sampler(normal.mean, normal.standard_deviation);
                                     return sampler(engine);
                                 },
                                 [&engine](const qmomkit::LogNormalDistribution& lognormal)
                                 {
                                     std::lognormal_distribution<double> sampler(lognormal.normal_mean, lognormal.normal_standard_deviation);
                                     return sampler(engine);
                                 },
                                 [&engine](const qmomkit::BetaDistribution& beta)
                                 {
                                     // A beta variate can be sampled as X / (X + Y)
                                     // for independent gamma variates with unit scale.
                                     std::gamma_distribution<double> left(beta.alpha, 1.0);
                                     std::gamma_distribution<double> right(beta.beta, 1.0);
                                     const double x = left(engine);
                                     const double y = right(engine);
                                     return x / (x + y);
                                 },
                                 [&engine](const qmomkit::GammaDistribution& gamma)
                                 {
                                     std::gamma_distribution<double> sampler(gamma.shape, gamma.scale);
                                     return sampler(engine);
                                 },
                                 [&engine](const qmomkit::ChiDistribution& chi)
                                 {
                                     std::chi_squared_distribution<double> sampler(chi.degrees_of_freedom);
                                     return std::sqrt(sampler(engine));
                                 },
                                 [&engine](const qmomkit::ChiSquaredDistribution& chi_squared)
                                 {
                                     std::chi_squared_distribution<double> sampler(chi_squared.degrees_of_freedom);
                                     return sampler(engine);
                                 },
                                 [&engine](const qmomkit::WeibullDistribution& weibull)
                                 {
                                     std::weibull_distribution<double> sampler(weibull.shape, weibull.scale);
                                     return sampler(engine);
                                 }},
                      distribution);
}

[[nodiscard]] SamplingResult check_finished_samples(Ensemble ensemble, const double number_density, const SamplingOptions& options)
{
    for (const Sample& sample : ensemble)
    {
        for (const double coordinate : sample)
        {
            if (!std::isfinite(coordinate))
            {
                return SamplingResult{std::move(ensemble), number_density, options.seed, SamplingStatus::non_finite_sample,
                                      "Sampling produced a non-finite coordinate."};
            }
        }
    }
    return SamplingResult{std::move(ensemble), number_density, options.seed, SamplingStatus::success, ""};
}

} // namespace

SamplingResult sample_distribution(const qmomkit::DistributionSpecification& distribution, const double number_density,
                                   const SamplingOptions& options)
{
    if (options.sample_count == 0)
    {
        return invalid_sampling_result(options, "Sample count must be positive.");
    }
    if (!is_nonnegative_finite(number_density))
    {
        return invalid_sampling_result(options, "Number density must be finite and non-negative.");
    }
    if (!valid_distribution(distribution))
    {
        return invalid_sampling_result(options, "Distribution parameters are invalid.");
    }

    // The seed defines the finite stochastic realization; density is kept
    // separate so every sample can later receive equal weight n0 / N.
    std::mt19937_64 engine(options.seed);
    Ensemble ensemble(options.sample_count);
    for (Sample& sample : ensemble)
    {
        sample = {sample_value(distribution, engine)};
    }

    return check_finished_samples(std::move(ensemble), number_density, options);
}

SamplingResult sample_mixture(const qmomkit::DistributionMixture& mixture, const SamplingOptions& options)
{
    if (options.sample_count == 0)
    {
        return invalid_sampling_result(options, "Sample count must be positive.");
    }
    if (mixture.components.empty())
    {
        return invalid_sampling_result(options, "Mixture must contain at least one component.");
    }

    std::vector<double> component_weights;
    std::vector<std::size_t> component_indices;
    double number_density = 0.0;
    for (std::size_t component = 0; component < mixture.components.size(); ++component)
    {
        const double density = mixture.components[component].number_density;
        if (!is_nonnegative_finite(density))
        {
            return invalid_sampling_result(options, "Mixture component density must be finite and non-negative.");
        }
        if (density == 0.0)
        {
            continue;
        }
        if (!valid_distribution(mixture.components[component].distribution))
        {
            return invalid_sampling_result(options, "Positive-density mixture component has invalid distribution parameters.");
        }
        // Component densities are absolute NDF contributions. They are used as
        // relative probabilities only for selecting the sampled component.
        component_indices.push_back(component);
        component_weights.push_back(density);
        number_density += density;
        if (!std::isfinite(number_density))
        {
            return invalid_sampling_result(options, "Mixture number density accumulation became non-finite.");
        }
    }
    if (component_indices.empty() || number_density <= 0.0)
    {
        return invalid_sampling_result(options, "Mixture must contain at least one positive-density component.");
    }

    std::mt19937_64 engine(options.seed);
    std::discrete_distribution<std::size_t> component_sampler(component_weights.begin(), component_weights.end());
    Ensemble ensemble(options.sample_count);
    for (Sample& sample : ensemble)
    {
        const std::size_t compact_component = component_sampler(engine);
        const std::size_t component = component_indices[compact_component];
        sample = {sample_value(mixture.components[component].distribution, engine)};
    }

    return check_finished_samples(std::move(ensemble), number_density, options);
}

SamplingResult sample_independent_product(const qmomkit::IndependentProductDistribution& distribution, const SamplingOptions& options)
{
    if (options.sample_count == 0)
    {
        return invalid_sampling_result(options, "Sample count must be positive.");
    }
    if (distribution.marginals.empty())
    {
        return invalid_sampling_result(options, "Independent product must contain at least one marginal.");
    }
    if (!is_nonnegative_finite(distribution.number_density))
    {
        return invalid_sampling_result(options, "Number density must be finite and non-negative.");
    }
    for (const qmomkit::DistributionSpecification& marginal : distribution.marginals)
    {
        if (!valid_distribution(marginal))
        {
            return invalid_sampling_result(options, "Independent-product marginal has invalid distribution parameters.");
        }
    }

    // Independent-product sampling draws each phase-space coordinate from its
    // marginal with the same engine stream. Correlations are intentionally not
    // introduced by this companion utility.
    std::mt19937_64 engine(options.seed);
    Ensemble ensemble(options.sample_count);
    for (Sample& sample : ensemble)
    {
        sample.resize(distribution.marginals.size());
        for (std::size_t dimension = 0; dimension < distribution.marginals.size(); ++dimension)
        {
            sample[dimension] = sample_value(distribution.marginals[dimension], engine);
        }
    }

    return check_finished_samples(std::move(ensemble), distribution.number_density, options);
}

ExplicitDiscreteSamplingResult sample_explicit_discrete_distribution(const ExplicitDiscreteDistribution& distribution, const SamplingOptions& options)
{
    if (options.sample_count == 0)
    {
        return invalid_explicit_sampling_result(options, "Sample count must be positive.");
    }

    const detail::ExplicitDistributionValidation validation = detail::validate_explicit_distribution(distribution);
    if (!validation.ok)
    {
        return invalid_explicit_sampling_result(options, validation.message);
    }

    // Absolute node weights define sampling probabilities after normalization;
    // their sum remains the number-density scale for moment estimates.
    std::mt19937_64 engine(options.seed);
    std::discrete_distribution<std::size_t> node_sampler(distribution.weights.begin(), distribution.weights.end());

    Ensemble ensemble(options.sample_count);
    for (Sample& sample : ensemble)
    {
        const std::size_t node = node_sampler(engine);
        sample = distribution.nodes[node];
    }

    return ExplicitDiscreteSamplingResult{std::move(ensemble), validation.number_density, options.seed, ComparisonStatus::success, ""};
}

} // namespace qmomkit::monte_carlo
