#include "particle_model.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{
namespace
{
constexpr double spherical_limit_tolerance = 0.001;

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

[[nodiscard]] double oblate_spheroid_shape_factor(const double aspect_ratio, const double reynolds_number)
{
    // Apply the non-spherical drag correction, using the spherical-limit form
    // when the aspect ratio is close enough to one to avoid the oblate formula.
    const double e = std::sqrt(std::abs(1.0 - aspect_ratio * aspect_ratio));
    double k0 = 1.0;
    if (std::abs(1.0 - aspect_ratio) >= spherical_limit_tolerance)
    {
        const double one_minus_chi_squared = 1.0 - aspect_ratio * aspect_ratio;
        const double atan_value = aspect_ratio < 1.0 ? std::atan(e / aspect_ratio) : 0.0;
        const double term1 = 8.0 / 3.0 / std::pow(aspect_ratio, 1.0 / 3.0);
        const double term2 = 2.0 * aspect_ratio / one_minus_chi_squared +
                             2.0 * (1.0 - 2.0 * aspect_ratio * aspect_ratio) / std::pow(one_minus_chi_squared, 1.5) * atan_value;
        k0 = term1 / term2;
    }

    return k0 + 0.15 * std::pow(aspect_ratio, 95.91) * std::pow(reynolds_number, 0.687) +
           0.2927 * std::pow(1.0 - aspect_ratio, 0.4374) * std::pow(reynolds_number, 0.7512);
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
    require_positive_finite(parameters.particle_density, "rho_p / particle_density");
    require_nonnegative_finite(parameters.gravity, "g / gravity");
    require_finite(parameters.fluid_velocity_amplitude_x, "U_AMPLITUDE / fluid_velocity_amplitude_x");
    require_finite(parameters.fluid_velocity_amplitude_y, "V_AMPLITUDE / fluid_velocity_amplitude_y");
    require_finite(parameters.fluid_velocity_amplitude_z, "W_AMPLITUDE / fluid_velocity_amplitude_z");
    require_finite(parameters.oscillation_angular_frequency, "OMEGA / oscillation_angular_frequency");
    require_nonnegative_finite(parameters.attrition_rate, "ATTRITION_RATE_CONST / attrition_rate");
    require_positive_finite(parameters.aspect_ratio_relaxation_time, "TAU_RELAX / aspect_ratio_relaxation_time");
}

qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters)
{
    validate_particle_model_parameters(parameters);

    return [parameters](const double time, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        if (dimensions != 5)
        {
            throw std::invalid_argument("The non-spherical particle model expects coordinates ordered "
                                        "as (d, chi, u, v, w).");
        }

        // The fluid state depends only on time, so it is evaluated once for the
        // whole node set rather than once per node.
        const double fluid_velocity_x = parameters.fluid_velocity_amplitude_x * std::cos(parameters.oscillation_angular_frequency * time);

        const double fluid_velocity_y = parameters.fluid_velocity_amplitude_y * std::sin(parameters.oscillation_angular_frequency * time);

        const double fluid_velocity_z = parameters.fluid_velocity_amplitude_z * std::cos(parameters.oscillation_angular_frequency * time);

        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            // Ordering defined in configuration.qmom: (d, chi, u, v, w).
            const std::size_t first = node * dimensions;
            const double d = coordinates[first];
            const double chi = coordinates[first + 1];
            const double u = coordinates[first + 2];
            const double v = coordinates[first + 3];
            const double w = coordinates[first + 4];

            const double relative_velocity_x = fluid_velocity_x - u;
            const double relative_velocity_y = fluid_velocity_y - v;
            const double relative_velocity_z = fluid_velocity_z - w;

            const double relative_velocity_magnitude = std::sqrt(
                relative_velocity_x * relative_velocity_x + relative_velocity_y * relative_velocity_y + relative_velocity_z * relative_velocity_z);

            const double reynolds_number = parameters.fluid_density * relative_velocity_magnitude * d / parameters.fluid_dynamic_viscosity;

            const double shape_factor = oblate_spheroid_shape_factor(chi, reynolds_number);

            /*
             * This is algebraically equivalent to evaluating the drag force
             *
             *   F_d = 1/2 rho_f A (24/Re) K Re |U_rel| U_rel
             *
             * and dividing by the spherical particle mass. The resulting
             * acceleration coefficient is 18 mu_f K / (rho_p d^2).
             */
            const double drag_rate = 18.0 * parameters.fluid_dynamic_viscosity * shape_factor / (parameters.particle_density * d * d);

            // Diameter attrition uses a saturated velocity response and a tanh taper
            // near the configured minimum diameter; no callback clipping is applied.
            const double relative_speed_squared = relative_velocity_magnitude * relative_velocity_magnitude;
            const double attrition_speed_squared =
                relative_speed_squared / (1.0 + relative_speed_squared / (parameters.attrition_velocity_scale * parameters.attrition_velocity_scale));
            const double attrition_taper = std::tanh((d - parameters.minimum_attrition_diameter) / parameters.attrition_transition_diameter);
            velocities[first] = -parameters.attrition_rate * attrition_speed_squared * attrition_taper;

            velocities[first + 1] = (1.0 - chi) / parameters.aspect_ratio_relaxation_time;

            velocities[first + 2] = drag_rate * relative_velocity_x;
            velocities[first + 3] = drag_rate * relative_velocity_y;
            velocities[first + 4] = drag_rate * relative_velocity_z - parameters.gravity;
        }
    };
}

} // namespace cqmom_5d_non_spherical_particles
