#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace qmom_1d_minimal
{

qmomkit::PhaseSpaceVelocity advection_velocity()
{
    return [](double, const std::vector<double>& coordinates)
    {
        if (coordinates.size() != 1)
        {
            throw std::invalid_argument("The minimal QMOM example expects one phase-space coordinate.");
        }

        const double diameter = coordinates[0];
        return std::vector<double>{1.0 - diameter};
    };
}

qmomkit::MomentSource fragmentation_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        constexpr double fragmentation_rate = 0.5;

        std::vector<double> derivative(context.size(), 0.0);
        for (std::size_t order = 0; order < context.size(); ++order)
        {
            const double daughter_factor = std::pow(2.0, 1.0 - static_cast<double>(order));
            derivative[order] = fragmentation_rate * (daughter_factor - 1.0) * context.value(order);
        }
        return derivative;
    };
}

} // namespace qmom_1d_minimal
