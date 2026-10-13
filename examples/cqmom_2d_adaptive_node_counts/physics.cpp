#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace cqmom_2d_adaptive_node_counts
{
namespace
{

constexpr double time_tolerance = 1.0e-12;

double binomial_coefficient(const std::size_t n, const std::size_t k)
{
    if (k > n)
    {
        return 0.0;
    }

    const std::size_t reduced_k = std::min(k, n - k);
    double coefficient = 1.0;
    for (std::size_t i = 1; i <= reduced_k; ++i)
    {
        coefficient *= static_cast<double>(n + 1 - i);
        coefficient /= static_cast<double>(i);
    }
    return coefficient;
}

} // namespace

bool source_pulse_active(const double time, const VarianceInjectionParameters& parameters)
{
    return time + time_tolerance >= parameters.pulse_start && time < parameters.pulse_end - time_tolerance;
}

qmomkit::PhaseSpaceVelocity advection_velocity(const VarianceInjectionParameters& parameters)
{
    // Velocities relax towards the equilibrium velocity; diameters do not move.
    return [parameters](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            const double u = coordinates[first + velocity_dimension];
            velocities[first] = 0.0;
            velocities[first + velocity_dimension] = -parameters.relaxation_rate * (u - parameters.equilibrium_velocity);
        }
    };
}

qmomkit::CqmomMomentSource variance_injection_source(const VarianceInjectionParameters& parameters)
{
    return [parameters](const double time, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> derivative(moments.size(), 0.0);
        if (!source_pulse_active(time, parameters))
        {
            return qmomkit::CqmomMomentRate{std::move(derivative)};
        }

        for (std::size_t entry = 0; entry < moments.size(); ++entry)
        {
            const qmomkit::MultiIndexView order = moments.index(entry);
            const std::size_t velocity_order = order[1];
            double redistribution_change = 0.0;
            for (std::size_t displacement_order = 2; displacement_order <= velocity_order; displacement_order += 2)
            {
                qmomkit::MultiIndex lower_order = order.to_multi_index();
                lower_order[1] = velocity_order - displacement_order;
                redistribution_change += binomial_coefficient(velocity_order, displacement_order) *
                                         std::pow(parameters.velocity_displacement, static_cast<double>(displacement_order)) *
                                         moments.moment(lower_order);
            }

            derivative[entry] = parameters.source_rate * redistribution_change;
        }

        return qmomkit::CqmomMomentRate{std::move(derivative)};
    };
}

} // namespace cqmom_2d_adaptive_node_counts
