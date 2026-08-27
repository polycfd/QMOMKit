#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace cqmom_2d_advection_breakage
{

constexpr double fragmentation_rate = 0.5;

[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

[[nodiscard]] qmomkit::CqmomMomentSource binary_fragmentation_source();

} // namespace cqmom_2d_advection_breakage
