#pragma once

#include "fluid_field.hpp"
#include "qmomkit/node_advection.hpp"

namespace cqmom_5d_non_spherical_particles
{

struct ParticleModelParameters
{
    double fluid_density = 1.2;
    double fluid_dynamic_viscosity = 1.8e-5;
    double particle_density = 2000.0;
    double gravity = 9.81;
    double attrition_rate = 1.0e-6;
    double minimum_attrition_diameter = 0.0;
    double attrition_transition_diameter = 1.0e-5;
    double attrition_velocity_scale = 5.0;
    double aspect_ratio_relaxation_time = 2.0;
};

[[nodiscard]] ParticleModelParameters particle_model_parameters();

// Returns the phase-space velocity for states ordered as (d, chi, u, v, w),
// capturing the parameters by value and the carrier-flow samples by reference.
// The samples must outlive the returned callback. The callback reads them on
// every invocation, so the driver can publish new host levels between steps
// without rebuilding the callback.
[[nodiscard]] qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters,
                                                                        const CarrierFlowSamples& carrier_flow);

// The callback keeps a reference to the samples, so a temporary would dangle.
// Binding one is rejected at compile time.
qmomkit::PhaseSpaceVelocity particle_phase_space_velocity(const ParticleModelParameters& parameters, CarrierFlowSamples&& carrier_flow) = delete;

} // namespace cqmom_5d_non_spherical_particles
