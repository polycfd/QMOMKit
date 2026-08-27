#include "droplet_model.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace cqmom_3d_droplets_condensation_nucleation
{
namespace
{
constexpr double kelvin_offset = 273.15;
constexpr double minimum_liquid_water_temperature = 273.15;
constexpr double low_reynolds_number = 1.0e-8;

struct DropletState
{
    double radius = 0.0;
    double vertical_velocity = 0.0;
    double temperature = 0.0;
};

struct DropletRates
{
    double radius_rate = 0.0;
    double vertical_velocity_rate = 0.0;
    double temperature_rate = 0.0;
};

struct AmbientConditions
{
    double temperature = 313.15;
    double relative_humidity = 1.0;
    double vertical_air_velocity = 0.0;
    double gravity = 9.81;
};

struct DropletModelParameters
{
    double water_density = 997.0;
    double water_heat_capacity = 4181.0;
    double surface_tension = 0.072;
    double water_molar_mass = 18.0153e-3;
    double gas_constant = 8.3145;
    double water_critical_temperature = 647.096;
};

struct MoistAirProperties
{
    double density = 0.0;
    double dynamic_viscosity = 0.0;
    double thermal_conductivity = 0.0;
    double heat_capacity = 0.0;
};

[[nodiscard]] bool is_positive_finite(const double value)
{
    return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] bool is_nonnegative_finite(const double value)
{
    return std::isfinite(value) && value >= 0.0;
}

void require_positive_finite(const double value, const char* name)
{
    if (!is_positive_finite(value))
    {
        throw std::invalid_argument(std::string("Droplet model value '") + name + "' must be supplied as a positive finite value.");
    }
}

void require_nonnegative_finite(const double value, const char* name)
{
    if (!is_nonnegative_finite(value))
    {
        throw std::invalid_argument(std::string("Droplet model value '") + name + "' must be supplied as a non-negative finite value.");
    }
}

void require_finite(const double value, const char* name)
{
    if (!std::isfinite(value))
    {
        throw std::invalid_argument(std::string("Droplet model value '") + name + "' must be supplied as a finite value.");
    }
}

[[nodiscard]] DropletModelParameters droplet_model_parameters()
{
    return {};
}

[[nodiscard]] double integer_power(const double base, const std::size_t exponent) noexcept
{
    double result = 1.0;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        result *= base;
    }
    return result;
}

void validate_temperature_for_air_properties(const double temperature, const char* name)
{
    require_positive_finite(temperature, name);
    if (temperature < minimum_liquid_water_temperature || temperature > 373.15)
    {
        throw std::invalid_argument(std::string("Droplet model value '") + name +
                                    "' must be in the liquid-water air-property correlation range [273.15 K, 373.15 K].");
    }
}

void validate_temperature_for_water_properties(const double temperature, const DropletModelParameters& parameters, const char* name)
{
    require_positive_finite(temperature, name);
    if (temperature < minimum_liquid_water_temperature || temperature >= parameters.water_critical_temperature)
    {
        throw std::invalid_argument(std::string("Droplet model value '") + name +
                                    "' must be in the liquid-water property correlation range [273.15 K, T_c).");
    }
}

[[nodiscard]] double polynomial_air_property(const double theta_celsius, const double relative_humidity, const double a1, const double b1,
                                             const double b2, const double c1, const double c2)
{
    return a1 + theta_celsius * (b1 + b2 * relative_humidity) + theta_celsius * theta_celsius * (c1 + c2 * relative_humidity);
}

[[nodiscard]] MoistAirProperties moist_air_properties(const AmbientConditions& ambient)
{
    const double theta = ambient.temperature - kelvin_offset;
    const double phi = ambient.relative_humidity;
    MoistAirProperties properties;
    properties.density = polynomial_air_property(theta, phi, 1.2932, -4.6688e-3, -8.6965e-4, 1.6165e-5, 2.2437e-5);
    properties.dynamic_viscosity = polynomial_air_property(theta, phi, 1.7192e-5, 4.9467e-8, -2.2425e-9, -4.2412e-11, -3.2389e-10);
    properties.thermal_conductivity = polynomial_air_property(theta, phi, 2.4031e-2, 7.5904e-5, -3.1199e-6, -4.8653e-8, -1.3015e-7);
    properties.heat_capacity = polynomial_air_property(theta, phi, 1006.4, 3.9654e-2, 2.0101, 4.15208e-4, -1.673e-1);
    assert(properties.density > 0.0);
    assert(properties.dynamic_viscosity > 0.0);
    assert(properties.thermal_conductivity > 0.0);
    assert(properties.heat_capacity > 0.0);
    return properties;
}

[[nodiscard]] double saturation_pressure(const double temperature)
{
    return std::exp(34.494 - 4924.99 / (temperature - 36.05)) * std::pow(temperature - 168.15, -1.57);
}

[[nodiscard]] double water_vapor_diffusivity(const double ambient_temperature)
{
    // Eq. (31) is used at the example's fixed ambient pressure p_inf = 101325 Pa.
    return 2.178e-5 * std::pow(ambient_temperature / kelvin_offset, 1.81);
}

[[nodiscard]] double latent_heat(const double droplet_temperature, const DropletModelParameters& parameters)
{
    const double theta_d = (parameters.water_critical_temperature - droplet_temperature) / droplet_temperature;
    const double theta_c = (parameters.water_critical_temperature - droplet_temperature) / parameters.water_critical_temperature;
    const double theta_c_0835 = std::pow(theta_c, 0.835);
    return 6.8477 * (std::pow(theta_d, 0.325) + 1.7383 * theta_c_0835) / (1.0 + 1.3913 * theta_c_0835) * parameters.gas_constant *
           parameters.water_critical_temperature / parameters.water_molar_mass;
}

[[nodiscard]] double reynolds_number(const DropletState& state, const AmbientConditions& ambient, const MoistAirProperties& air)
{
    return 2.0 * state.radius * air.density * std::abs(ambient.vertical_air_velocity - state.vertical_velocity) / air.dynamic_viscosity;
}

[[nodiscard]] double weber_number(const DropletState& state, const AmbientConditions& ambient, const MoistAirProperties& air,
                                  const DropletModelParameters& parameters)
{
    const double relative_velocity = ambient.vertical_air_velocity - state.vertical_velocity;
    return 2.0 * state.radius * air.density * relative_velocity * relative_velocity / parameters.surface_tension;
}

[[nodiscard]] double transfer_number(const double reynolds, const double transport_property_ratio)
{
    const double criterion = std::sqrt(reynolds) * std::cbrt(transport_property_ratio);
    if (criterion < 1.4)
    {
        return 2.0 + 0.316 * reynolds * std::pow(transport_property_ratio, 2.0 / 3.0);
    }
    return 1.56 + 0.616 * std::sqrt(reynolds) * std::cbrt(transport_property_ratio);
}

[[nodiscard]] double sherwood_number(const double reynolds, const MoistAirProperties& air, const double diffusivity)
{
    return transfer_number(reynolds, air.dynamic_viscosity / (air.density * diffusivity));
}

[[nodiscard]] double base_nusselt_number(const double reynolds, const MoistAirProperties& air)
{
    return transfer_number(reynolds, air.dynamic_viscosity * air.heat_capacity / air.thermal_conductivity);
}

[[nodiscard]] double radius_rate_at_temperature(const DropletState& state, const AmbientConditions& ambient, const DropletModelParameters& parameters,
                                                const double sherwood, const double droplet_temperature)
{
    const double vapor_drive = ambient.relative_humidity * saturation_pressure(ambient.temperature) / ambient.temperature -
                               saturation_pressure(droplet_temperature) / droplet_temperature;
    const double diffusivity = water_vapor_diffusivity(ambient.temperature);
    return sherwood * diffusivity * parameters.water_molar_mass * vapor_drive /
           (2.0 * state.radius * parameters.water_density * parameters.gas_constant);
}

[[nodiscard]] double nusselt_number(const MoistAirProperties& air, const double reynolds)
{
    // This condensation-only example uses the base heat-transfer correlation
    // from Eq. (23) and deliberately omits the evaporation-oriented wet-bulb
    // correction from Eqs. (22), (25), and (26).
    return base_nusselt_number(reynolds, air);
}

[[nodiscard]] double vertical_drag_acceleration(const DropletState& state, const AmbientConditions& ambient, const MoistAirProperties& air,
                                                const DropletModelParameters& parameters)
{
    const double relative_velocity = ambient.vertical_air_velocity - state.vertical_velocity;
    const double relative_speed = std::abs(relative_velocity);
    if (relative_speed == 0.0)
    {
        return 0.0;
    }

    const double reynolds = reynolds_number(state, ambient, air);
    if (reynolds < low_reynolds_number)
    {
        return 4.5 * air.dynamic_viscosity * relative_velocity / (parameters.water_density * state.radius * state.radius);
    }

    const double weber = weber_number(state, ambient, air, parameters);
    const double deformation = 1.0 + 2.632 * (1.0 - 1.0 / ((1.0 + 0.19 * std::sqrt(weber)) * (1.0 + 0.19 * std::sqrt(weber))));
    const double cd0 = 24.0 / reynolds * (1.0 + 0.15 * std::pow(reynolds, 0.687) + 0.0175 * reynolds / (1.0 + 42500.0 * std::pow(reynolds, -1.16)));
    const double drag_coefficient = cd0 * deformation;
    return 3.0 * air.density * drag_coefficient * relative_speed * relative_velocity / (8.0 * parameters.water_density * state.radius);
}

void validate_state(const DropletState& state, const DropletModelParameters& parameters)
{
    require_positive_finite(state.radius, "radius");
    require_finite(state.vertical_velocity, "vertical_velocity");
    validate_temperature_for_water_properties(state.temperature, parameters, "temperature");
}

void validate_ambient_conditions(const AmbientConditions& ambient)
{
    validate_temperature_for_air_properties(ambient.temperature, "ambient.temperature");
    require_finite(ambient.vertical_air_velocity, "ambient.vertical_air_velocity");
    require_nonnegative_finite(ambient.gravity, "ambient.gravity");
    if (!std::isfinite(ambient.relative_humidity) || ambient.relative_humidity < 0.0 || ambient.relative_humidity > 1.0)
    {
        throw std::invalid_argument("Droplet model value 'ambient.relative_humidity' must be supplied as a finite value in [0, 1].");
    }
}

void validate_droplet_model_parameters(const DropletModelParameters& parameters)
{
    require_positive_finite(parameters.water_density, "water_density");
    require_positive_finite(parameters.water_heat_capacity, "water_heat_capacity");
    require_positive_finite(parameters.surface_tension, "surface_tension");
    require_positive_finite(parameters.water_molar_mass, "water_molar_mass");
    require_positive_finite(parameters.gas_constant, "gas_constant");
    require_positive_finite(parameters.water_critical_temperature, "water_critical_temperature");
}

DropletRates evaluate_droplet_rates(const DropletState& state, const AmbientConditions& ambient, const DropletModelParameters& parameters)
{
    validate_droplet_model_parameters(parameters);
    validate_ambient_conditions(ambient);
    validate_state(state, parameters);

    const MoistAirProperties air = moist_air_properties(ambient);
    const double reynolds = reynolds_number(state, ambient, air);
    const double diffusivity = water_vapor_diffusivity(ambient.temperature);
    const double sherwood = sherwood_number(reynolds, air, diffusivity);
    const double radius_rate = radius_rate_at_temperature(state, ambient, parameters, sherwood, state.temperature);
    const double nu = nusselt_number(air, reynolds);

    DropletRates rates;
    rates.radius_rate = radius_rate;
    rates.vertical_velocity_rate =
        -ambient.gravity * (1.0 - air.density / parameters.water_density) + vertical_drag_acceleration(state, ambient, air, parameters);
    rates.temperature_rate =
        3.0 / (state.radius * parameters.water_heat_capacity) *
        (nu * air.thermal_conductivity * (ambient.temperature - state.temperature) / (2.0 * state.radius * parameters.water_density) +
         latent_heat(state.temperature, parameters) * radius_rate);

    if (!std::isfinite(rates.radius_rate) || !std::isfinite(rates.vertical_velocity_rate) || !std::isfinite(rates.temperature_rate))
    {
        throw std::runtime_error("Droplet model returned a non-finite rate.");
    }

    return rates;
}

} // namespace

qmomkit::PhaseSpaceVelocity droplet_phase_space_velocity()
{
    const AmbientConditions ambient;
    const DropletModelParameters parameters = droplet_model_parameters();
    validate_ambient_conditions(ambient);
    validate_droplet_model_parameters(parameters);

    return [ambient, parameters](double, const std::vector<double>& coordinates)
    {
        if (coordinates.size() != 3)
        {
            throw std::invalid_argument("The condensing-droplet example expects phase-space coordinates ordered as (r,w,T).");
        }

        const DropletRates rates = evaluate_droplet_rates({coordinates[0], coordinates[1], coordinates[2]}, ambient, parameters);
        return std::vector<double>{rates.radius_rate, rates.vertical_velocity_rate, rates.temperature_rate};
    };
}

NucleationModel make_nucleation_model(const double initial_number_density)
{
    require_positive_finite(initial_number_density, "initial_number_density");
    require_positive_finite(nucleation_parameters.relative_population_increase, "nucleation.relative_population_increase");
    const double duration = nucleation_parameters.end_time - nucleation_parameters.start_time;
    require_positive_finite(duration, "nucleation duration");

    // A constant birth flux adds the configured relative population increase
    // over the active nucleation interval.
    return {nucleation_parameters.relative_population_increase * initial_number_density / duration};
}

void validate_nucleation_interval(const double initial_time, const double final_time)
{
    require_finite(initial_time, "initial_time");
    require_finite(final_time, "final_time");
    if (nucleation_parameters.start_time < initial_time || nucleation_parameters.end_time <= nucleation_parameters.start_time ||
        nucleation_parameters.end_time > final_time)
    {
        throw std::invalid_argument("The nucleation interval must lie within the driver time interval.");
    }
}

double limit_step_to_nucleation_interval(const double time, const double maximum_step)
{
    require_finite(time, "time");
    require_positive_finite(maximum_step, "maximum_step");
    if (time < nucleation_parameters.start_time)
    {
        return std::min(maximum_step, nucleation_parameters.start_time - time);
    }
    if (time < nucleation_parameters.end_time)
    {
        return std::min(maximum_step, nucleation_parameters.end_time - time);
    }
    return maximum_step;
}

bool nucleation_is_active(const double time) noexcept
{
    return time >= nucleation_parameters.start_time && time < nucleation_parameters.end_time;
}

qmomkit::CqmomMomentBoundaryContribution nucleation_boundary_contribution(const NucleationModel& model)
{
    require_positive_finite(model.flux, "nucleation flux");

    return [flux = model.flux](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        // A birth at (r_n, w_n, T_n) contributes J r_n^a w_n^b T_n^c to
        // dM[a,b,c]/dt. This is a phase-space boundary birth rather than node motion.
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);
            if (order.size() != 3)
            {
                throw std::runtime_error("The nucleation contribution requires phase-space order (r,w,T).");
            }
            rates[entry] = flux * integer_power(nucleation_parameters.radius, order[0]) *
                           integer_power(nucleation_parameters.vertical_velocity, order[1]) *
                           integer_power(nucleation_parameters.temperature, order[2]);
        }
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };
}

qmomkit::monte_carlo::PopulationBirthResult inject_nucleated_droplets(const NucleationModel& model, const qmomkit::monte_carlo::Ensemble& ensemble,
                                                                      const double number_density, const double time, const double time_step,
                                                                      const unsigned long long random_seed)
{
    require_positive_finite(model.flux, "nucleation flux");
    require_positive_finite(time_step, "nucleation time_step");
    if (!nucleation_is_active(time))
    {
        throw std::invalid_argument("Nucleation injection requested outside the active interval.");
    }

    // The Monte Carlo birth sampler uses the same deterministic birth state and
    // integrated flux J Delta t as the CQMOM boundary contribution.
    const qmomkit::monte_carlo::BirthCoordinateSampler sample_nucleated_droplet = [](std::mt19937_64&)
    {
        return qmomkit::monte_carlo::Sample{nucleation_parameters.radius, nucleation_parameters.vertical_velocity, nucleation_parameters.temperature};
    };
    return qmomkit::monte_carlo::inject_population(ensemble, number_density, model.flux * time_step, sample_nucleated_droplet, random_seed);
}

} // namespace cqmom_3d_droplets_condensation_nucleation
