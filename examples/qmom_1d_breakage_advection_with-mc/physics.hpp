#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace qmom_1d_breakage_advection_with_mc
{

constexpr double fragmentation_rate = 0.5;

/**
 * User-editable phase-space velocity for xi = V.
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

/**
 * User-editable moment source for binary fragmentation V -> V/2.
 */
[[nodiscard]] qmomkit::QmomMomentSource binary_fragmentation_source();

} // namespace qmom_1d_breakage_advection_with_mc
