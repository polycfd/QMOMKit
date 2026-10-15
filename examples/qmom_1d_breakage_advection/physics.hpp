#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace qmom_1d_breakage_advection
{

/**
 * User-editable phase-space velocity for xi = V, the particle volume.
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

/**
 * User-editable moment source for size-dependent binary fragmentation
 * V -> V/2, written over the nodes of the source-stage quadrature.
 */
[[nodiscard]] qmomkit::QmomMomentSource binary_fragmentation_source();

} // namespace qmom_1d_breakage_advection
