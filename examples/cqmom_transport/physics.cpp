#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cqmom_transport
{
namespace
{

constexpr double fragmentation_rate = 0.5;
constexpr double injection_rate = 0.1;
constexpr double injection_size = 0.1;

} // namespace

qmomkit::PhaseSpaceVelocity continuous_velocity()
{
    return [](double, const std::vector<double>& coordinates)
    {
        if (coordinates.empty())
        {
            throw std::invalid_argument("The CQMOM transport example expects at least one phase-space coordinate.");
        }

        std::vector<double> velocity(coordinates.size(), 0.0);
        velocity[0] = 1.0 - coordinates[0];
        for (std::size_t dimension = 1; dimension < coordinates.size(); ++dimension)
        {
            velocity[dimension] = -coordinates[dimension];
        }
        return velocity;
    };
}

qmomkit::CqmomMomentSource binary_fragmentation_source()
{
    return [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> derivative(context.size(), 0.0);
        for (std::size_t index = 0; index < context.size(); ++index)
        {
            const qmomkit::MultiIndex& order = context.index(index);
            if (order.empty())
            {
                throw std::invalid_argument("The CQMOM fragmentation source received an empty moment index.");
            }

            const double breakage_coordinate_order = static_cast<double>(order[0]);
            // Binary fragmentation x0 -> x0/2 leaves all other coordinates unchanged.
            const double daughter_factor = std::pow(2.0, 1.0 - breakage_coordinate_order);
            derivative[index] = fragmentation_rate * (daughter_factor - 1.0) * context.value(index);
        }
        return qmomkit::CqmomMomentRate{std::move(derivative)};
    };
}

qmomkit::CqmomMomentBoundaryContribution injection_boundary_contribution()
{
    return [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);

        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);

            if (order.empty())
            {
                throw std::invalid_argument("The CQMOM injection contribution received an empty moment index.");
            }

            // Birth of new entities in the first phase-space dimension only.
            // The contribution to M_alpha is J * d_b^(alpha_0).
            rate[entry] = injection_rate * std::pow(injection_size, static_cast<double>(order[0]));
        }

        return qmomkit::CqmomMomentRate{std::move(rate)};
    };
}

} // namespace cqmom_transport
