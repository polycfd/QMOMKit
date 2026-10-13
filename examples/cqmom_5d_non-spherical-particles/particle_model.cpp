#include "particle_model.hpp"

#include <cmath>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{
namespace
{
constexpr double spherical_limit_tolerance = 0.001;

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

qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters, const CarrierFlowSamples& carrier_flow)
{
    // Capturing the reference parameter by reference refers to the driver's
    // samples, not to the parameter, so the capture stays valid after this
    // function returns.
    return [parameters, &carrier_flow](const double time, const std::vector<double>& coordinates, std::vector<double>& velocities,
                                       const std::size_t dimensions)
    {
        // The carrier flow is host data rather than part of the particle model. It
        // is reconstructed once for the whole node set, at the stage time QMOMKit
        // supplies, from the two host levels the driver published.
        const BackgroundFluidVelocity fluid = carrier_flow_at(carrier_flow, time);
        const double fluid_velocity_x = fluid.u;
        const double fluid_velocity_y = fluid.v;
        const double fluid_velocity_z = fluid.w;

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
