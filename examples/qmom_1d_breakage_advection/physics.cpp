#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace qmom_1d_breakage_advection
{
namespace
{

constexpr double fragmentation_rate = 0.5;

} // namespace

qmomkit::PhaseSpaceVelocity advection_velocity()
{
    return [](double, const std::vector<double>& coordinates)
    {
        if (coordinates.size() != 1)
        {
            throw std::invalid_argument("The 1D QMOM example expects one phase-space coordinate.");
        }

        const double diameter = coordinates[0];
        return std::vector<double>{1.0 - diameter};
    };
}

qmomkit::MomentSource binary_fragmentation_source()
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
        return source;
    };
}

} // namespace qmom_1d_breakage_advection
