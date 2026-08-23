#include "physics.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cqmom_spherical_breakup_moments_5d
{
namespace
{

constexpr double pi = 3.14159265358979323846;
constexpr double gas_density = 1.2;
constexpr double gas_viscosity = 1.8e-5;
constexpr double gas_thermal_conductivity = 0.026;
constexpr double gas_prandtl_number = 0.71;
constexpr double gas_temperature = 300.0;
constexpr double particle_density = 2000.0;
constexpr double particle_heat_capacity = 1000.0;
constexpr double gas_u = 2.0;
constexpr double gas_v = 1.5;
constexpr double gas_w = 0.8;
constexpr double breakup_rate = 20.0;

double diameter_from_volume(const double volume)
{
    if (!(volume > 0.0))
    {
        throw std::invalid_argument("Particle volume must remain positive.");
    }
    return std::cbrt(6.0 * volume / pi);
}

} // namespace

qmomkit::PhaseSpaceVelocity particle_velocity()
{
    return [](double, const std::vector<double>& coordinates)
    {
        if (coordinates.size() != 5)
        {
            throw std::invalid_argument("Thermal particle coordinates must be ordered as (V,T,u,v,w).");
        }
        const double d = diameter_from_volume(coordinates[0]);
        const double du = gas_u - coordinates[2];
        const double dv = gas_v - coordinates[3];
        const double dw = gas_w - coordinates[4];
        const double slip = std::sqrt(du * du + dv * dv + dw * dw);
        const double reynolds = gas_density * d * slip / gas_viscosity;
        const double drag = reynolds > 1.0e-14 ? 3.0 * gas_density / (4.0 * particle_density * d) *
                                                     (24.0 / reynolds * (1.0 + 0.15 * std::pow(reynolds, 0.687))) * slip
                                               : 18.0 * gas_viscosity / (particle_density * d * d);
        const double nusselt = 2.0 + 0.6 * std::sqrt(reynolds) * std::cbrt(gas_prandtl_number);
        const double temperature_rate =
            6.0 * gas_thermal_conductivity * nusselt / (particle_density * particle_heat_capacity * d * d) * (gas_temperature - coordinates[1]);
        return std::vector<double>{0.0, temperature_rate, drag * du, drag * dv, drag * dw};
    };
}

qmomkit::CqmomMomentSource breakup_source()
{
    return [](double, const qmomkit::CqmomSourceContext& context)
    {
        if (context.dimension_count() != 5)
        {
            throw std::invalid_argument("The thermal breakup source requires phase space (V,T,u,v,w).");
        }
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);
            const double p = static_cast<double>(order[0]);
            const double daughter_factor = std::pow(2.0, 1.0 - p) - 1.0;
            if (daughter_factor == 0.0)
            {
                continue;
            }
            rates[entry] = breakup_rate * daughter_factor * context.value(entry);
        }
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };
}

void apply_breakup(qmomkit::monte_carlo::Ensemble& ensemble, const double interval, std::mt19937_64& generator)
{
    qmomkit::monte_carlo::Ensemble completed;
    completed.reserve(2 * ensemble.size());
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        std::vector<std::pair<qmomkit::monte_carlo::Sample, double>> pending{{sample, interval}};
        while (!pending.empty())
        {
            auto [parent, remaining] = std::move(pending.back());
            pending.pop_back();
            if (parent.size() != 5)
            {
                throw std::runtime_error("MC particles must be ordered as (V,T,u,v,w).");
            }
            std::exponential_distribution<double> wait(breakup_rate);
            const double event_time = wait(generator);
            if (event_time >= remaining)
            {
                completed.push_back(std::move(parent));
                continue;
            }
            parent[0] *= 0.5;
            const double remaining_after_event = remaining - event_time;
            pending.push_back({parent, remaining_after_event});
            pending.push_back({std::move(parent), remaining_after_event});
        }
    }
    ensemble = std::move(completed);
}

} // namespace cqmom_spherical_breakup_moments_5d
