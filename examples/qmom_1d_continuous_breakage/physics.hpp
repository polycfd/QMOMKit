#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace qmom_1d_continuous_breakage
{

/**
 * User-editable continuous phase-space velocity for xi = d.
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity continuous_velocity();

/**
 * User-editable discontinuous source for binary fragmentation d -> d/2.
 */
[[nodiscard]] qmomkit::MomentSource binary_fragmentation_source();

} // namespace qmom_1d_continuous_breakage
