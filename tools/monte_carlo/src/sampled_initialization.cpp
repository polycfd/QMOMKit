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

[[nodiscard]] SelectedEnsembleMomentEstimate skipped_selected_moments(const SamplingResult& sampling)
{
    return SelectedEnsembleMomentEstimate{{}, {}, {}, EnsembleMomentStatus::invalid_input, "Sampling failed: " + sampling.message};
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

SelectedSampledInitialization sample_initialization(const qmomkit::IndependentProductDistribution& distribution,
                                                    const qmomkit::RawMomentRequest& request, const SamplingOptions& options)
{
    // Sampling is the same draw as the dense overload for the same options; only
    // the moment estimator differs, so the two agree on any shared index.
    SamplingResult sampling = sample_independent_product(distribution, options);
    if (!sampling.ok())
    {
        return {sampling, skipped_selected_moments(sampling)};
    }

    SelectedEnsembleMomentEstimate moments = compute_requested_ensemble_moments(sampling.ensemble, sampling.number_density, request);
    return {std::move(sampling), std::move(moments)};
}

} // namespace qmomkit::monte_carlo
