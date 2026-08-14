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
    // A 250 micrometre particle loses 12.5% of its diameter over one second.
    double surface_recession_rate = 3.125e-5;
    double carrier_temperature = 0.950;
    double diameter_regularization = 1.0e-8;
};

[[nodiscard]] ParticleModelParameters particle_model_parameters();

void validate_particle_model_parameters(const ParticleModelParameters& parameters);

// Returns the phase-space velocity for states ordered as
// (d, u, v, w, T), with T in kK, capturing the validated parameters by value.
[[nodiscard]] qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters);

} // namespace cqmom_reactive_spherical_particles
