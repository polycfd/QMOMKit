#include "monte_carlo/sampled_initialization.hpp"

#include <string>
#include <utility>

namespace qmomkit::monte_carlo
{
namespace
{

// Preserve the paired result shape while making failed sampling impossible to
// mistake for a valid zero-moment realization.
[[nodiscard]] OneDimensionalMomentEstimate skipped_one_dimensional_moments(const SamplingResult& sampling)
{
    return OneDimensionalMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, "Sampling failed: " + sampling.message};
}

[[nodiscard]] TensorMomentEstimate skipped_tensor_moments(const SamplingResult& sampling)
{
    return TensorMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, "Sampling failed: " + sampling.message};
}

} // namespace

OneDimensionalSampledInitialization sample_initialization(const qmomkit::DistributionSpecification& distribution, const double number_density,
                                                          const std::size_t maximum_order, const SamplingOptions& options)
{
    SamplingResult sampling = sample_distribution(distribution, number_density, options);
    if (!sampling.ok())
    {
        return {sampling, skipped_one_dimensional_moments(sampling)};
    }

    OneDimensionalMomentEstimate moments = compute_ensemble_moments(sampling.ensemble, sampling.number_density, maximum_order);
    return {std::move(sampling), std::move(moments)};
}

OneDimensionalSampledInitialization sample_initialization(const qmomkit::DistributionMixture& mixture, const std::size_t maximum_order,
                                                          const SamplingOptions& options)
{
    SamplingResult sampling = sample_mixture(mixture, options);
    if (!sampling.ok())
    {
        return {sampling, skipped_one_dimensional_moments(sampling)};
    }

    OneDimensionalMomentEstimate moments = compute_ensemble_moments(sampling.ensemble, sampling.number_density, maximum_order);
    return {std::move(sampling), std::move(moments)};
}

TensorSampledInitialization sample_initialization(const qmomkit::IndependentProductDistribution& distribution,
                                                  const std::vector<std::size_t>& maximum_orders, const SamplingOptions& options)
{
    SamplingResult sampling = sample_independent_product(distribution, options);
    if (!sampling.ok())
    {
        return {sampling, skipped_tensor_moments(sampling)};
    }

    TensorMomentEstimate moments = compute_ensemble_moments(sampling.ensemble, sampling.number_density, maximum_orders);
    return {std::move(sampling), std::move(moments)};
}

} // namespace qmomkit::monte_carlo
