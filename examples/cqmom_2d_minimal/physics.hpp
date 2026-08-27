#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace cqmom_2d_minimal
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

[[nodiscard]] qmomkit::CqmomMomentSource binary_fragmentation_source();

[[nodiscard]] qmomkit::CqmomMomentBoundaryContribution injection_boundary_contribution();

} // namespace cqmom_2d_minimal
