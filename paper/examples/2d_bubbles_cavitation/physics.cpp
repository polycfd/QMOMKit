#include "physics.hpp"

#include <cmath>
#include <vector>

namespace paper_2d_bubbles_cavitation
{
// Dimensionless Rayleigh-Plesset dynamics of a gas bubble, after Bryngelson et
// al. (2023), https://doi.org/10.1016/j.jcp.2023.111917: dR/dt = U and
// dU/dt = [(R_o/R)^(3 gamma) - 1/C_p - 2/(We R_o) (R_o/R - (R_o/R)^(3 gamma))
//          - 3/2 U^2 - 4/Re U/R] / R,
// with the radius R, the radial velocity U and the time all dimensionless.
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

} // namespace paper_2d_bubbles_cavitation
