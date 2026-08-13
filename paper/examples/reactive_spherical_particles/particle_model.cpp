#include "particle_model.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cqmom_reactive_spherical_particles
{
namespace
{
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
        throw std::invalid_argument(std::string("Particle model parameter '") + name + "' must be supplied as a positive finite value.");
    }
}

void require_nonnegative_finite(const double value, const char* name)
{
    if (!is_nonnegative_finite(value))
    {
        throw std::invalid_argument(std::string("Particle model parameter '") + name + "' must be supplied as a non-negative finite value.");
    }
}

void require_finite(const double value, const char* name)
{
    if (!std::isfinite(value))
    {
        throw std::invalid_argument(std::string("Particle model parameter '") + name + "' must be supplied as a finite value.");
    }
}

} // namespace

ParticleModelParameters particle_model_parameters()
{
    ParticleModelParameters parameters;

    return parameters;
}

void validate_particle_model_parameters(const ParticleModelParameters& parameters)
{
    require_positive_finite(parameters.fluid_density, "rho_f / fluid_density");
    require_positive_finite(parameters.fluid_dynamic_viscosity, "mu_f / fluid_dynamic_viscosity");
    require_positive_finite(parameters.fluid_thermal_conductivity, "k_f / fluid_thermal_conductivity");
    require_positive_finite(parameters.fluid_prandtl_number, "Pr / fluid_prandtl_number");
    require_positive_finite(parameters.particle_density, "rho_p / particle_density");
    require_positive_finite(parameters.particle_specific_heat, "c_p / particle_specific_heat");
    require_positive_finite(parameters.particle_thermal_conductivity, "k_p / particle_thermal_conductivity");
    require_nonnegative_finite(parameters.gravity, "g / gravity");
    require_nonnegative_finite(parameters.attrition_rate, "ATTRITION_RATE_CONST / attrition_rate");
    require_positive_finite(parameters.attrition_transition_diameter, "attrition_transition_diameter");
    require_positive_finite(parameters.attrition_velocity_scale, "attrition_velocity_scale");
    require_nonnegative_finite(parameters.surface_rate, "surface_rate");
    require_positive_finite(parameters.surface_temperature_transition_width, "surface_temperature_transition_width");
    require_positive_finite(parameters.surface_diameter_transition, "surface_diameter_transition");
    require_finite(parameters.carrier_temperature_center, "carrier_temperature_center");
    require_nonnegative_finite(parameters.carrier_temperature_amplitude, "carrier_temperature_amplitude");
    require_finite(parameters.carrier_temperature_angular_frequency, "carrier_temperature_angular_frequency");
    require_positive_finite(parameters.diameter_regularization, "diameter_regularization");
}

qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters)
{
    validate_particle_model_parameters(parameters);

    return [parameters](const double time, const std::vector<double>& state)
    {
        if (state.size() != 5)
        {
            throw std::invalid_argument("The reactive spherical-particle model expects coordinates ordered as (d, u, v, w, T).");
        }

        // Ordering defined in configuration.qmom: (d, u, v, w, T), with T in kK.
        const double d = state[0];
        const double u = state[1];
        const double v = state[2];
        const double w = state[3];
        const double temperature = state[4];
        const double d_regularized = std::sqrt(d * d + parameters.diameter_regularization * parameters.diameter_regularization);

        const double fluid_velocity_x = 2.0 * std::cos(time);
        const double fluid_velocity_y = 2.0 * std::sin(time);
        constexpr double fluid_velocity_z = 0.0;

        const double relative_velocity_x = fluid_velocity_x - u;
        const double relative_velocity_y = fluid_velocity_y - v;
        const double relative_velocity_z = fluid_velocity_z - w;

        const double relative_velocity_magnitude = std::sqrt(relative_velocity_x * relative_velocity_x + relative_velocity_y * relative_velocity_y +
                                                             relative_velocity_z * relative_velocity_z);

        const double reynolds_number = parameters.fluid_density * d_regularized * relative_velocity_magnitude / parameters.fluid_dynamic_viscosity;
        const double drag_correction = 1.0 + 0.15 * std::pow(reynolds_number, 0.687);
        const double drag_rate =
            18.0 * parameters.fluid_dynamic_viscosity * drag_correction / (parameters.particle_density * d_regularized * d_regularized);
        const double nusselt_number = 2.0 + 0.6 * std::sqrt(reynolds_number) * std::cbrt(parameters.fluid_prandtl_number);
        const double heat_transfer_rate = 6.0 * nusselt_number * parameters.fluid_thermal_conductivity /
                                          (parameters.particle_density * parameters.particle_specific_heat * d_regularized * d_regularized);
        const double carrier_temperature =
            parameters.carrier_temperature_center +
            parameters.carrier_temperature_amplitude * std::sin(parameters.carrier_temperature_angular_frequency * time);

        std::vector<double> rate(5);

        const double relative_speed_squared = relative_velocity_magnitude * relative_velocity_magnitude;
        const double attrition_speed_squared =
            relative_speed_squared / (1.0 + relative_speed_squared / (parameters.attrition_velocity_scale * parameters.attrition_velocity_scale));
        const double attrition_taper = std::tanh((d - parameters.minimum_attrition_diameter) / parameters.attrition_transition_diameter);
        const double surface_rate =
            -parameters.surface_rate *
            std::tanh((temperature - parameters.carrier_temperature_center) / parameters.surface_temperature_transition_width) *
            std::tanh(d / parameters.surface_diameter_transition);
        rate[0] = surface_rate - parameters.attrition_rate * attrition_speed_squared * attrition_taper;
        rate[1] = drag_rate * relative_velocity_x;
        rate[2] = drag_rate * relative_velocity_y;
        rate[3] = drag_rate * relative_velocity_z - parameters.gravity;
        rate[4] = heat_transfer_rate * (carrier_temperature - temperature);

        return rate;
    };
}

} // namespace cqmom_reactive_spherical_particles
