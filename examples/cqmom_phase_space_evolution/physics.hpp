#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace cqmom_phase_space_evolution
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

[[nodiscard]] qmomkit::CqmomMomentSource binary_fragmentation_source();

[[nodiscard]] qmomkit::CqmomMomentBoundaryContribution injection_boundary_contribution();

} // namespace cqmom_phase_space_evolution
