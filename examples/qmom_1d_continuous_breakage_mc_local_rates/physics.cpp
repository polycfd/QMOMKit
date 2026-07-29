#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace qmom_1d_continuous_breakage_mc_local_rates
{
double breakage_event_rate() noexcept
{
    return fragmentation_rate;
}

double daughter_diameter(const double parent_diameter) noexcept
{
    return 0.5 * parent_diameter;
}

double breakage_probability(const double step_size) noexcept
{
    return 1.0 - std::exp(-breakage_event_rate() * step_size);
}

double nucleation_rate() noexcept
{
    return constant_nucleation_rate;
}

double nucleation_diameter() noexcept
{
    return fixed_nucleation_diameter;
}

double nucleation_moment_rate(const std::size_t order) noexcept
{
    return nucleation_rate() * std::pow(nucleation_diameter(), static_cast<double>(order));
}

qmomkit::PhaseSpaceVelocity continuous_velocity()
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
    return [](double, const std::vector<double>& moments)
    {
        std::vector<double> source(moments.size(), 0.0);
        for (std::size_t order = 0; order < moments.size(); ++order)
        {
            // Binary fragmentation d -> d/2 gives dM_k/dt = lambda * (2^(1-k) - 1) * M_k.
            const double daughter_factor = 2.0 * std::pow(daughter_diameter(1.0), static_cast<double>(order));
            source[order] = breakage_event_rate() * (daughter_factor - 1.0) * moments[order];
        }
        return source;
    };
}

qmomkit::MomentBoundaryContribution boundary_contribution()
{
    return [](double, const std::vector<double>& moments)
    {
        std::vector<double> boundary_rates(moments.size(), 0.0);
        for (std::size_t order = 0; order < boundary_rates.size(); ++order)
        {
            boundary_rates[order] = nucleation_moment_rate(order);
        }
        return qmomkit::MomentBoundaryRate{boundary_rates};
    };
}

} // namespace qmom_1d_continuous_breakage_mc_local_rates
