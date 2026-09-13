#include "physics.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace qmom_1d_explicit_strang_splitting
{

qmomkit::PhaseSpaceVelocity advection_velocity()
{
    return qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& coordinates)
        {
            constexpr double relaxation_rate = 500.0;
            if (coordinates.size() != 1)
            {
                throw std::invalid_argument("The explicit Strang example expects one phase-space coordinate.");
            }
            return std::vector<double>{-relaxation_rate * coordinates[0]};
        });
}

qmomkit::MomentSource variance_injection_source()
{
    return [](const double time, const qmomkit::QmomSourceContext& context)
    {
        constexpr double source_start = 0.16;
        constexpr double source_end = 0.25;
        constexpr double variance_injection_rate = 12.0;

        std::vector<double> derivative(context.size(), 0.0);
        if (time >= source_start && time < source_end && derivative.size() > 2)
        {
            derivative[2] = variance_injection_rate;
        }
        return derivative;
    };
}

} // namespace qmom_1d_explicit_strang_splitting
