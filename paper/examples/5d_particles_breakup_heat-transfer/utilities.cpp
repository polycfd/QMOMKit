#include "utilities.hpp"

#include "monte_carlo/moments.hpp"

#include <stdexcept>

namespace paper_5d_particles_breakup_heat_transfer
{

std::vector<std::size_t> maximum_orders(const std::vector<std::size_t>& extents)
{
    std::vector<std::size_t> orders;
    orders.reserve(extents.size());
    for (const std::size_t extent : extents)
    {
        orders.push_back(extent - 1);
    }
    return orders;
}

qmomkit::RawMomentTensor ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, const double particle_weight,
                                          const std::vector<std::size_t>& extents)
{
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, particle_weight * static_cast<double>(ensemble.size()), maximum_orders(extents));
    if (!estimate.ok())
    {
        throw std::runtime_error(estimate.message);
    }
    return estimate.moments;
}

} // namespace paper_5d_particles_breakup_heat_transfer
