#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

#include <cstddef>
#include <string>

namespace cqmom_2d_variance_reinjection
{

struct VarianceInjectionParameters
{
    double relaxation_rate = 5.0;
    double equilibrium_velocity = 0.0;
    double pulse_start = 4.0;
    double pulse_end = 4.05;
    double source_rate = 20.0;
    double velocity_displacement = 0.10;
};

[[nodiscard]] VarianceInjectionParameters default_variance_injection_parameters();

[[nodiscard]] std::size_t velocity_phase_space_dimension();

[[nodiscard]] bool source_pulse_active(double time, const VarianceInjectionParameters& parameters);

[[nodiscard]] qmomkit::PhaseSpaceVelocity continuous_velocity(const VarianceInjectionParameters& parameters);

[[nodiscard]] qmomkit::CqmomMomentSource variance_injection_source(const VarianceInjectionParameters& parameters);

[[nodiscard]] bool validate_variance_injection_source_rates(const VarianceInjectionParameters& parameters, std::string& message);

} // namespace cqmom_2d_variance_reinjection
