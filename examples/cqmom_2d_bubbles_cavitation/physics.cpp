#include "physics.hpp"

#include <cmath>
#include <vector>

namespace cqmom_2d_bubbles_cavitation
{
qmomkit::PhaseSpaceVelocity advection_velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            const double radius = coordinates[first];
            const double velocity = coordinates[first + 1];
            const double radius_ratio = reference_radius / radius;
            const double compressed_ratio = std::pow(radius_ratio, 3.0 * heat_capacity_ratio);
            const double acceleration =
                (compressed_ratio - 1.0 / pressure_coefficient - (2.0 / (weber_number * reference_radius)) * (radius_ratio - compressed_ratio) -
                 1.5 * velocity * velocity - (4.0 / reynolds_number) * velocity / radius) /
                radius;
            velocities[first] = velocity;
            velocities[first + 1] = acceleration;
        }
    };
}

} // namespace cqmom_2d_bubbles_cavitation
