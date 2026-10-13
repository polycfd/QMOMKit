#include "particles_mc.hpp"

#include "monte_carlo/moments.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace paper_2d_bubbles_cavitation
{
namespace
{

std::vector<std::size_t> maximum_orders_from_extents(const std::vector<std::size_t>& extents)
{
    std::vector<std::size_t> maximum_orders;
    maximum_orders.reserve(extents.size());
    for (const std::size_t extent : extents)
    {
        maximum_orders.push_back(extent - 1);
    }
    return maximum_orders;
}

} // namespace

qmomkit::RawMomentTensor compute_particle_moments(const qmomkit::monte_carlo::Ensemble& ensemble, const double particle_weight,
                                                  const std::vector<std::size_t>& moment_extents)
{
    const double number_density = particle_weight * static_cast<double>(ensemble.size());
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, maximum_orders_from_extents(moment_extents));
    if (!estimate.ok())
    {
        throw std::runtime_error(estimate.message);
    }
    return estimate.moments;
}

} // namespace paper_2d_bubbles_cavitation
