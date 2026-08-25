#pragma once

#include "qmomkit/node_advection.hpp"

namespace paper_2d_bubbles_cavitation
{

constexpr double reference_radius = 1.0;
constexpr double pressure_coefficient = 0.3;
constexpr double heat_capacity_ratio = 1.4;
constexpr double reynolds_number = 10.0;
constexpr double weber_number = 1.39;

[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

} // namespace paper_2d_bubbles_cavitation
