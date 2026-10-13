#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

#include <cstddef>

namespace cqmom_2d_adaptive_node_counts
{

// The velocity u is phase-space coordinate 1 of (d, u).
constexpr std::size_t velocity_dimension = 1;

struct VarianceInjectionParameters
{
    double relaxation_rate = 5.0;
    double equilibrium_velocity = 0.5;
    double pulse_start = 4.0;
    double pulse_end = 4.05;
    double source_rate = 20.0;
    double velocity_displacement = 0.10;
};

[[nodiscard]] bool source_pulse_active(double time, const VarianceInjectionParameters& parameters);

[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity(const VarianceInjectionParameters& parameters);

[[nodiscard]] qmomkit::CqmomMomentSource variance_injection_source(const VarianceInjectionParameters& parameters);

} // namespace cqmom_2d_adaptive_node_counts
