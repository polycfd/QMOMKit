#pragma once

#include "qmomkit/node_advection.hpp"

namespace cqmom_non_spherical_particles
{

struct ParticleModelParameters
{
    double fluid_density = 1.2;
    double fluid_dynamic_viscosity = 1.8e-5;
    double particle_density = 2000.0;
    double gravity = 9.81;
    double fluid_velocity_amplitude_x = 2.0;
    double fluid_velocity_amplitude_y = 2.0;
    double fluid_velocity_amplitude_z = 0.0;
    double oscillation_angular_frequency = 1.0;
    double attrition_rate = 0.0; // 1.0e-3;
    double aspect_ratio_relaxation_time = 2.0;
};

[[nodiscard]] ParticleModelParameters particle_model_parameters();

void validate_particle_model_parameters(const ParticleModelParameters& parameters);

[[nodiscard]] qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters);

} // namespace cqmom_non_spherical_particles
