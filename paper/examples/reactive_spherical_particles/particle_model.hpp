#pragma once

#include "qmomkit/node_advection.hpp"

namespace cqmom_reactive_spherical_particles
{

struct ParticleModelParameters
{
    double fluid_density = 1.2;
    double fluid_dynamic_viscosity = 1.8e-5;
    double fluid_thermal_conductivity = 0.03;
    double fluid_prandtl_number = 0.71;
    double particle_density = 1800.0;
    double particle_specific_heat = 900.0;
    double particle_thermal_conductivity = 1.5;
    double gravity = 9.81;
    double attrition_rate = 1.0e-4;
    double minimum_attrition_diameter = 0.0;
    double attrition_transition_diameter = 1.0e-5;
    double attrition_velocity_scale = 5.0;
    double surface_rate = 1.2e-4;
    double surface_temperature_transition_width = 0.010;
    double surface_diameter_transition = 1.0e-5;
    double carrier_temperature_center = 0.900;
    double carrier_temperature_amplitude = 0.020;
    double carrier_temperature_angular_frequency = 12.566370614359172;
    double diameter_regularization = 1.0e-8;
};

[[nodiscard]] ParticleModelParameters particle_model_parameters();

void validate_particle_model_parameters(const ParticleModelParameters& parameters);

// Returns the phase-space velocity for states ordered as
// (d, u, v, w, T), with T in kK, capturing the validated parameters by value.
[[nodiscard]] qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters);

} // namespace cqmom_reactive_spherical_particles
