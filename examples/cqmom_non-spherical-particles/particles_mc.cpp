#include "particles_mc.hpp"

#include "qmomkit/multi_index.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace cqmom_non_spherical_particles
{
namespace
{

[[nodiscard]] double integer_power(double base, std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t count = 0; count < exponent; ++count)
    {
        value *= base;
    }
    return value;
}

void require_particle_sample(const qmomkit::monte_carlo::Sample& sample)
{
    if (sample.size() != 5)
    {
        throw std::runtime_error("The 5D MC particle example expects samples ordered as (d, chi, u, v, w).");
    }
}

} // namespace

qmomkit::RawMomentTensor compute_ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, const double particle_weight,
                                                  const std::vector<std::size_t>& moment_extents)
{
    if (moment_extents.size() != 5)
    {
        throw std::runtime_error("The 5D MC particle example expects moment extents ordered as (d, chi, u, v, w).");
    }

    const std::size_t moment_count = qmomkit::product_extents(moment_extents);
    std::vector<double> moments(moment_count, 0.0);
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        require_particle_sample(sample);
        for (std::size_t flat_index = 0; flat_index < moment_count; ++flat_index)
        {
            const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, moment_extents);
            double monomial = 1.0;
            for (std::size_t dimension = 0; dimension < sample.size(); ++dimension)
            {
                monomial *= integer_power(sample[dimension], order[dimension]);
            }
            moments[flat_index] += particle_weight * monomial;
        }
    }
    return qmomkit::RawMomentTensor{moment_extents, std::move(moments)};
}

} // namespace cqmom_non_spherical_particles
