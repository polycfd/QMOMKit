#include "particles_mc.hpp"

#include "qmomkit/multi_index.hpp"
#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace cqmom_2d_continuous_breakage
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

void require_two_dimensional_sample(const qmomkit::monte_carlo::Sample& sample)
{
    if (sample.size() != 2)
    {
        throw std::runtime_error("The 2D MC particle example expects samples ordered as (d,u).");
    }
}

} // namespace

void apply_binary_breakage(qmomkit::monte_carlo::Ensemble& ensemble, const double step_size, std::mt19937_64& generator)
{
    const double breakage_probability = 1.0 - std::exp(-fragmentation_rate * step_size);
    std::bernoulli_distribution breaks(breakage_probability);

    const std::size_t initial_particle_count = ensemble.size();
    ensemble.reserve(ensemble.size() + initial_particle_count);
    for (std::size_t particle_index = 0; particle_index < initial_particle_count; ++particle_index)
    {
        qmomkit::monte_carlo::Sample& parent = ensemble[particle_index];
        require_two_dimensional_sample(parent);

        if (!breaks(generator))
        {
            continue;
        }

        parent[0] *= 0.5;
        ensemble.push_back(qmomkit::monte_carlo::Sample{parent[0], parent[1]});
    }
}

qmomkit::RawMomentTensor compute_particle_moments(const qmomkit::monte_carlo::Ensemble& ensemble, const double particle_weight,
                                                    const std::vector<std::size_t>& moment_extents)
{
    const std::size_t moment_count = qmomkit::product_extents(moment_extents);
    std::vector<double> moments(moment_count, 0.0);
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        require_two_dimensional_sample(sample);
        for (std::size_t flat_index = 0; flat_index < moment_count; ++flat_index)
        {
            const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, moment_extents);
            moments[flat_index] += particle_weight * integer_power(sample[0], order[0]) * integer_power(sample[1], order[1]);
        }
    }
    return qmomkit::RawMomentTensor{moment_extents, std::move(moments)};
}

} // namespace cqmom_2d_continuous_breakage
