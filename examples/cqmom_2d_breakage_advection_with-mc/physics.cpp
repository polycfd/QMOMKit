#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace cqmom_2d_advection_breakage
{
qmomkit::PhaseSpaceVelocity advection_velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            const double volume = coordinates[first];
            const double u = coordinates[first + 1];
            velocities[first] = 1.0 - volume;
            velocities[first + 1] = -u;
        }
    };
}

qmomkit::CqmomMomentSource binary_fragmentation_source()
{
    return [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> source(moments.size(), 0.0);
        for (std::size_t index = 0; index < moments.size(); ++index)
        {
            const qmomkit::MultiIndexView order = moments.index(index);
            const double volume_order = static_cast<double>(order[0]);
            // Binary fragmentation V -> V/2 leaves u unchanged:
            // dM_ij/dt = lambda * (2^(1-i) - 1) * M_ij.
            const double daughter_factor = std::pow(2.0, 1.0 - volume_order);
            source[index] = fragmentation_rate * (daughter_factor - 1.0) * moments.value(index);
        }
        return qmomkit::CqmomMomentRate{std::move(source)};
    };
}

} // namespace cqmom_2d_advection_breakage
