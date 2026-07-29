#include "physics.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace cqmom_2d_vanishing_variance
{

qmomkit::PhaseSpaceVelocity continuous_velocity()
{
    return [](double, const std::vector<double>& coordinates)
    {
        const double u = coordinates[1];
        return std::vector<double>{0.0, -5.0 * u};
    };
}

qmomkit::CqmomMomentSource zero_moment_source()
{
    return [](double, const qmomkit::CqmomSourceContext& moments) { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
}

} // namespace cqmom_2d_vanishing_variance
