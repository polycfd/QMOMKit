#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace cqmom_transport
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity continuous_velocity();

[[nodiscard]] qmomkit::CqmomMomentSource binary_fragmentation_source();

[[nodiscard]] qmomkit::CqmomMomentBoundaryContribution injection_boundary_contribution();

} // namespace cqmom_transport
