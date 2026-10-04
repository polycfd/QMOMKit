#pragma once

#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"
#include "qmomkit/distributions.hpp"

#include <cstddef>
#include <vector>

namespace qmomkit::monte_carlo
{

/**
 * Pair a finite sample realization with moments estimated from that same
 * realization.
 *
 * The helper exists for comparisons that need a common stochastic initial
 * state:
 *
 *     distribution
 *         -> sample ensemble
 *         -> sample moments
 *         -> QMOM/CQMOM initialization
 *
 * The sample ensemble can also be evolved directly as an independent Monte
 * Carlo reference. A physical application may interpret samples as particles,
 * but QMOMKitMonteCarlo treats them only as phase-space realizations.
 */
template <class MomentEstimate> struct SampledInitializationResult
{
    SamplingResult sampling;
    MomentEstimate moments;

    [[nodiscard]] bool ok() const noexcept
    {
        return sampling.ok() && moments.ok();
    }
};

using OneDimensionalSampledInitialization = SampledInitializationResult<OneDimensionalMomentEstimate>;
using TensorSampledInitialization = SampledInitializationResult<TensorMomentEstimate>;
using SelectedSampledInitialization = SampledInitializationResult<SelectedEnsembleMomentEstimate>;

/**
 * Sample a one-dimensional distribution and compute raw moments from the same
 * finite sample ensemble.
 *
 * Moments preserve the absolute number-density convention:
 *
 *     M[k] = n0 / N * sum_i xi_i^k
 */
[[nodiscard]] OneDimensionalSampledInitialization sample_initialization(const qmomkit::DistributionSpecification& distribution, double number_density,
                                                                        std::size_t maximum_order, const SamplingOptions& options);

/**
 * Sample a one-dimensional absolute-density mixture and compute raw moments
 * from the same finite sample ensemble.
 *
 * Mixture component densities choose the sampled component probabilities; the
 * returned moments are scaled by the summed absolute number density.
 */
[[nodiscard]] OneDimensionalSampledInitialization sample_initialization(const qmomkit::DistributionMixture& mixture, std::size_t maximum_order,
                                                                        const SamplingOptions& options);

/**
 * Sample an independent multidimensional product and compute the raw moment
 * tensor from the same finite sample ensemble.
 *
 * Tensor moments use phase-space dimension order and row-major storage.
 */
[[nodiscard]] TensorSampledInitialization sample_initialization(const qmomkit::IndependentProductDistribution& distribution,
                                                                const std::vector<std::size_t>& maximum_orders, const SamplingOptions& options);

/**
 * Sample an independent multidimensional product and estimate only the
 * requested raw moments from the same finite sample ensemble.
 *
 * The request is evaluated exactly as `compute_requested_ensemble_moments`
 * does: request order and duplicate entries are preserved, and no full moment
 * tensor is constructed. A caller that needs only the moments a CQMOM hierarchy
 * requires therefore evaluates that set rather than the full rectangle, which
 * for higher-dimensional hierarchies is most of the estimation work.
 *
 * For identical `SamplingOptions` this draws the same ensemble as the tensor
 * overload, so a requested moment equals the corresponding entry of the full
 * tensor estimated from that ensemble.
 */
[[nodiscard]] SelectedSampledInitialization sample_initialization(const qmomkit::IndependentProductDistribution& distribution,
                                                                  const qmomkit::RawMomentRequest& request, const SamplingOptions& options);

} // namespace qmomkit::monte_carlo
