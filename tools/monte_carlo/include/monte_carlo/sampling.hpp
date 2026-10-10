#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/distributions.hpp"

#include <string>

namespace qmomkit::monte_carlo
{

/**
 * Status for constructing a finite sample ensemble from a distribution.
 *
 * A successful result means the requested number of samples was drawn and all
 * coordinates are finite. It does not mean the finite ensemble exactly matches
 * analytical moments; sampling uncertainty is intrinsic to Monte Carlo
 * realization.
 */
enum class SamplingStatus
{
    success,
    invalid_input,
    non_finite_sample
};

struct SamplingResult
{
    // Finite stochastic realization of the requested distribution.
    Ensemble ensemble;
    // Absolute number density carried by the ensemble when moments are formed.
    double number_density = 0.0;
    // Seed used to initialize the pseudo-random engine.
    unsigned long long seed = 0;
    SamplingStatus status = SamplingStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == SamplingStatus::success;
    }
};

/**
 * Draw a one-dimensional sample ensemble from a supported analytical
 * distribution.
 *
 * Each returned sample has one coordinate. The input distribution describes a
 * normalized shape; `number_density` supplies the absolute density scale used
 * later by moment estimation. With `N = options.sample_count`, the standard
 * moment convention treats each sample as having statistical weight
 * `number_density / N`.
 *
 * Sampling uses a deterministic pseudo-random engine initialized from
 * `options.seed`. Reusing the same distribution, sample count, seed, compiler,
 * and standard library gives reproducible samples for this implementation, but
 * Monte Carlo samples remain a finite stochastic realization rather than exact
 * distribution moments.
 */
[[nodiscard]] SamplingResult sample_distribution(const qmomkit::DistributionSpecification& distribution, double number_density,
                                                 const SamplingOptions& options);

/**
 * Draw a one-dimensional sample ensemble from an absolute-density mixture.
 *
 * Component densities are interpreted as absolute contributions to total
 * number density. During sampling they are normalized only to choose the
 * component:
 *
 *     P(component r) = n_r / sum_s n_s
 *
 * The returned `number_density` is the sum of positive component densities.
 * Samples are then equally weighted by that total density during moment
 * estimation.
 */
[[nodiscard]] SamplingResult sample_mixture(const qmomkit::DistributionMixture& mixture, const SamplingOptions& options);

/**
 * Draw a multidimensional sample ensemble from independent marginals.
 *
 * Each sample has one coordinate per marginal, stored in QMOMKit phase-space
 * dimension order. The implementation assumes independence between marginals;
 * it does not model correlations, copulas, or application-specific physics.
 * The distribution's `number_density` supplies the absolute density scale for
 * later raw-moment estimation.
 */
[[nodiscard]] SamplingResult sample_independent_product(const qmomkit::IndependentProductDistribution& distribution, const SamplingOptions& options);

} // namespace qmomkit::monte_carlo
