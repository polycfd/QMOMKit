#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace qmom_1d_breakage_advection
{

/**
 * User-editable phase-space velocity for xi = d.
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

/**
 * User-editable moment source for binary fragmentation d -> d/2.
 */
[[nodiscard]] qmomkit::QmomMomentSource binary_fragmentation_source();

} // namespace qmom_1d_breakage_advection
