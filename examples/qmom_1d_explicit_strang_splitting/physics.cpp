#include "physics.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace qmom_1d_explicit_strang_splitting
{

qmomkit::PhaseSpaceVelocity advection_velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, std::size_t)
    {
        constexpr double relaxation_rate = 500.0;
        constexpr double equilibrium_velocity = 0.5;
        for (std::size_t node = 0; node < coordinates.size(); ++node)
        {
            velocities[node] = -relaxation_rate * (coordinates[node] - equilibrium_velocity);
        }
    };
}

qmomkit::QmomMomentSource variance_injection_source()
{
    return [](const double time, const qmomkit::QmomSourceContext& context)
    {
        constexpr double source_start = 0.16;
        constexpr double source_end = 0.25;
        constexpr double variance_injection_rate = 12.0;

        std::vector<double> derivative(context.size(), 0.0);
        if (time >= source_start && time < source_end)
        {
            derivative[2] = variance_injection_rate;
        }
        return qmomkit::QmomMomentRate{std::move(derivative)};
    };
}

} // namespace qmom_1d_explicit_strang_splitting
