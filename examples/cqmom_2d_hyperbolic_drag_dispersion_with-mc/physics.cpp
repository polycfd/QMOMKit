#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace cqmom_2d_hyperbolic_drag_dispersion
{

double response_time(const double diameter)
{
    return particle_density * diameter * diameter / (18.0 * fluid_viscosity);
}

double dispersion_coefficient(const double diameter)
{
    return turbulent_velocity_scale * turbulent_velocity_scale / response_time(diameter);
}

qmomkit::PhaseSpaceVelocity drag_velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            const double u = coordinates[first];
            const double d = coordinates[first + 1];
            velocities[first] = (flow_velocity - u) / response_time(d);
            velocities[first + 1] = 0.0;
        }
    };
}

qmomkit::CqmomMomentSource dispersion_source()
{
    return [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        moments.for_each_terminal_node(
            [&](const double weight, const std::vector<double>& coordinates)
            {
                const double u = coordinates[0];
                const double d = coordinates[1];
                const double diffusion = weight * dispersion_coefficient(d);
                for (std::size_t entry = 0; entry < moments.size(); ++entry)
                {
                    const qmomkit::MultiIndexView order = moments.index(entry);
                    const std::size_t k = order[0];
                    if (k >= 2)
                    {
                        rate[entry] += static_cast<double>(k * (k - 1)) * diffusion * std::pow(u, static_cast<int>(k - 2)) *
                                       std::pow(d, static_cast<int>(order[1]));
                    }
                }
            });
        return qmomkit::CqmomMomentRate{std::move(rate)};
    };
}

} // namespace cqmom_2d_hyperbolic_drag_dispersion
