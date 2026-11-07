#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace hyqmom_1d_drag_bgk_collisions
{
namespace
{

// Glass beads of radius 1 mm in air. The drag coefficient is nearly constant
// for particle Reynolds numbers above about 10^3, here above about 7.5 m/s, and
// held constant for slower particles too.
constexpr double drag_coefficient = 0.44;   // dimensionless
constexpr double gas_density = 1.2;         // kg/m^3
constexpr double particle_density = 2500.0; // kg/m^3
constexpr double particle_radius = 1.0e-3;  // m
constexpr double relaxation_time = 0.1;     // s

// du/dt = -c |u| u, with c = 3 C_D rho_g / (8 rho_p r_p), in 1/m.
constexpr double drag_constant = 3.0 * drag_coefficient * gas_density / (8.0 * particle_density * particle_radius);

} // namespace

qmomkit::PhaseSpaceVelocity drag_velocity()
{
    return [](double, const std::vector<double>& velocities, std::vector<double>& accelerations, std::size_t)
    {
        for (std::size_t node = 0; node < velocities.size(); ++node)
        {
            accelerations[node] = -drag_constant * std::abs(velocities[node]) * velocities[node];
        }
    };
}

qmomkit::QmomMomentSource bgk_collision_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        // BGK relaxation dm_k/dt = (m_k^eq - m_k) / theta, where m^eq are the
        // moments of the normal distribution with the current m0, mean and
        // variance, m_k^eq = mean m_{k-1}^eq + (k - 1) variance m_{k-2}^eq. It
        // leaves m0, m1 and m2, number, momentum and energy, unchanged.
        const double number_density = context.value(0);
        const double mean = context.value(1) / number_density;
        const double variance = context.value(2) / number_density - mean * mean;

        std::vector<double> rate(context.size(), 0.0);
        double previous_equilibrium_moment = 0.0;
        double equilibrium_moment = number_density;
        for (std::size_t order = 0; order < context.size(); ++order)
        {
            rate[order] = (equilibrium_moment - context.value(order)) / relaxation_time;
            const double next_equilibrium_moment = mean * equilibrium_moment + static_cast<double>(order) * variance * previous_equilibrium_moment;
            previous_equilibrium_moment = equilibrium_moment;
            equilibrium_moment = next_equilibrium_moment;
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
    };
}

} // namespace hyqmom_1d_drag_bgk_collisions
