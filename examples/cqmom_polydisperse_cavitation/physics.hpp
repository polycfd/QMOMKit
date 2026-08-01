#pragma once

#include "qmomkit/node_advection.hpp"

namespace cqmom_polydisperse_cavitation
{

constexpr double reference_radius = 1.0;
constexpr double pressure_coefficient = 0.3;
constexpr double heat_capacity_ratio = 1.4;
constexpr double reynolds_number = 10.0;
constexpr double weber_number = 1.39;

[[nodiscard]] qmomkit::PhaseSpaceVelocity continuous_velocity();

} // namespace cqmom_polydisperse_cavitation
