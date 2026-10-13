#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace qmom_1d_breakage_advection_with_mc
{
qmomkit::PhaseSpaceVelocity advection_velocity()
{
    // Diameters relax towards d = 1: smaller particles grow, larger ones shrink.
    return [](double, const std::vector<double>& diameters, std::vector<double>& velocities, std::size_t)
    {
        for (std::size_t node = 0; node < diameters.size(); ++node)
        {
            velocities[node] = 1.0 - diameters[node];
        }
    };
}

qmomkit::QmomMomentSource binary_fragmentation_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> source(context.size(), 0.0);
        for (std::size_t order = 0; order < context.size(); ++order)
        {
            // Binary fragmentation d -> d/2 gives dM_k/dt = lambda * (2^(1-k) - 1) * M_k.
            const double daughter_factor = std::pow(2.0, 1.0 - static_cast<double>(order));
            source[order] = fragmentation_rate * (daughter_factor - 1.0) * context.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(source)};
    };
}

} // namespace qmom_1d_breakage_advection_with_mc
