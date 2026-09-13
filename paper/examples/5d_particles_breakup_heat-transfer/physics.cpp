#include "physics.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace paper_5d_particles_breakup_heat_transfer
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
// These constants define the deliberately simple breakup temperature mapping
// used to demonstrate generalized-moment source closure.
constexpr double breakup_rate = 10.0;
constexpr double reference_diameter = 1.0e-4;
constexpr double breakup_temperature_relaxation = 1.0e-2;
constexpr double reference_volume = pi * reference_diameter * reference_diameter * reference_diameter / 6.0;

double diameter_from_volume(const double volume)
{
    if (!(volume > 0.0))
    {
        throw std::invalid_argument("Particle volume must remain positive.");
    }
    return std::cbrt(6.0 * volume / pi);
}

double binomial_coefficient(const std::size_t n, const std::size_t k)
{
    double coefficient = 1.0;
    for (std::size_t factor = 1; factor <= k; ++factor)
    {
        coefficient *= static_cast<double>(n - k + factor) / static_cast<double>(factor);
    }
    return coefficient;
}

double daughter_temperature(const double daughter_volume, const double parent_temperature)
{
    // The daughter differs from its parent through a small, surface-related
    // temperature redistribution. For a sphere, V/S scales as V^(1/3).
    const double relaxation = breakup_temperature_relaxation * std::cbrt(daughter_volume / reference_volume);
    return parent_temperature + relaxation * (gas_temperature - parent_temperature);
}

} // namespace

qmomkit::PhaseSpaceVelocity particle_velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        if (dimensions != 5)
        {
            throw std::invalid_argument("Thermal particle coordinates must be ordered as (V,T,u,v,w).");
        }
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            const double d = diameter_from_volume(coordinates[first]);
            const double du = gas_u - coordinates[first + 2];
            const double dv = gas_v - coordinates[first + 3];
            const double dw = gas_w - coordinates[first + 4];
            const double slip = std::sqrt(du * du + dv * dv + dw * dw);
            const double reynolds = gas_density * d * slip / gas_viscosity;
            const double drag =
                3.0 * gas_density / (4.0 * particle_density * d) * (24.0 / reynolds * (1.0 + 0.15 * std::pow(reynolds, 0.687))) * slip;
            const double nusselt = 2.0 + 0.6 * std::sqrt(reynolds) * std::cbrt(gas_prandtl_number);
            const double temperature_rate = 6.0 * gas_thermal_conductivity * nusselt / (particle_density * particle_heat_capacity * d * d) *
                                            (gas_temperature - coordinates[first + 1]);
            velocities[first] = 0.0;
            velocities[first + 1] = temperature_rate;
            velocities[first + 2] = drag * du;
            velocities[first + 3] = drag * dv;
            velocities[first + 4] = drag * dw;
        }
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
        // Expanding the daughter-temperature mapping leaves integer volume
        // terms available through context.moment(...) and exposes fractional
        // V^(1/3) terms that require generalized_moment(...).
        const double daughter_temperature_factor = breakup_temperature_relaxation * std::cbrt(0.5 / reference_volume);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndexView order = context.index(entry);
            const std::size_t p = order[0];
            const std::size_t q = order[1];
            double daughter_temperature_moment = 0.0;
            for (std::size_t temperature_difference_order = 0; temperature_difference_order <= q; ++temperature_difference_order)
            {
                const double volume_factor = std::pow(daughter_temperature_factor, static_cast<double>(temperature_difference_order));
                for (std::size_t parent_temperature_order = 0; parent_temperature_order <= temperature_difference_order; ++parent_temperature_order)
                {
                    const std::size_t temperature_order = q - temperature_difference_order + parent_temperature_order;
                    const double coefficient =
                        binomial_coefficient(q, temperature_difference_order) *
                        binomial_coefficient(temperature_difference_order, parent_temperature_order) * volume_factor *
                        std::pow(gas_temperature, static_cast<double>(temperature_difference_order - parent_temperature_order)) *
                        (parent_temperature_order % 2 == 0 ? 1.0 : -1.0);

                    if (temperature_difference_order % 3 == 0)
                    {
                        // Integer volume powers are ordinary raw moments.
                        daughter_temperature_moment +=
                            coefficient * context.moment({p + temperature_difference_order / 3, temperature_order, order[2], order[3], order[4]});
                    }
                    else
                    {
                        // Fractional volume powers are derived from the
                        // source-stage quadrature without changing the
                        // transported integer moment layout.
                        daughter_temperature_moment +=
                            coefficient *
                            context.generalized_moment({static_cast<double>(p) + static_cast<double>(temperature_difference_order) / 3.0,
                                                        static_cast<double>(temperature_order), static_cast<double>(order[2]),
                                                        static_cast<double>(order[3]), static_cast<double>(order[4])});
                    }
                }
            }
            rates[entry] = breakup_rate * (std::pow(2.0, 1.0 - static_cast<double>(p)) * daughter_temperature_moment - context.moment(order));
        }
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };
}

void apply_breakup(qmomkit::monte_carlo::Ensemble& ensemble, const double interval, std::mt19937_64& generator)
{
    // The MC reference samples the same binary daughter mapping as the CQMOM
    // source; only the population representation differs.
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
            parent[1] = daughter_temperature(parent[0], parent[1]);
            const double remaining_after_event = remaining - event_time;
            pending.push_back({parent, remaining_after_event});
            pending.push_back({std::move(parent), remaining_after_event});
        }
    }
    ensemble = std::move(completed);
}

} // namespace paper_5d_particles_breakup_heat_transfer
