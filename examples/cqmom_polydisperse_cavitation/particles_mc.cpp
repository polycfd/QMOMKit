#include "particles_mc.hpp"

#include "monte_carlo/moments.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cqmom_polydisperse_cavitation
{
namespace
{

void require_two_dimensional_sample(const qmomkit::monte_carlo::Sample& sample)
{
    if (sample.size() != 2)
    {
        throw std::runtime_error("The 2D MC particle example expects samples ordered as (R,U).");
    }
}

std::vector<std::size_t> maximum_orders_from_extents(const std::vector<std::size_t>& extents)
{
    std::vector<std::size_t> maximum_orders;
    maximum_orders.reserve(extents.size());
    for (const std::size_t extent : extents)
    {
        if (extent == 0)
        {
            throw std::invalid_argument("The 2D MC particle example requires non-empty moment extents.");
        }
        maximum_orders.push_back(extent - 1);
    }
    return maximum_orders;
}

} // namespace

qmomkit::RawMomentTensor compute_particle_moments(const qmomkit::monte_carlo::Ensemble& ensemble, const double particle_weight,
                                                  const std::vector<std::size_t>& moment_extents)
{
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        require_two_dimensional_sample(sample);
    }
    const double number_density = particle_weight * static_cast<double>(ensemble.size());
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, maximum_orders_from_extents(moment_extents));
    if (!estimate.ok())
    {
        throw std::runtime_error(estimate.message);
    }
    return estimate.moments;
}

} // namespace cqmom_polydisperse_cavitation
