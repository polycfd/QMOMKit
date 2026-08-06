#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

#include <cstddef>

namespace qmom_1d_breakage_advection_mc_local_rates
{

constexpr double fragmentation_rate = 0.5;
constexpr double constant_nucleation_rate = 0.05;
constexpr double fixed_nucleation_diameter = 0.2;

/**
 * Physical breakage event rate for binary fragmentation.
 */
[[nodiscard]] double breakage_event_rate() noexcept;

/**
 * Daughter state for one binary fragmentation event.
 */
[[nodiscard]] double daughter_diameter(double parent_diameter) noexcept;

/**
 * Event probability over one particle-level breakage step.
 */
[[nodiscard]] double breakage_probability(double step_size) noexcept;

/**
 * Constant physical number-rate for lower-boundary nucleation.
 */
[[nodiscard]] double nucleation_rate() noexcept;

/**
 * Diameter assigned to newly nucleated particles.
 */
[[nodiscard]] double nucleation_diameter() noexcept;

/**
 * Raw-moment contribution from the nucleation boundary model.
 */
[[nodiscard]] double nucleation_moment_rate(std::size_t order) noexcept;

/**
 * User-editable phase-space velocity for xi = d.
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();

/**
 * User-editable moment source for binary fragmentation d -> d/2.
 */
[[nodiscard]] qmomkit::MomentSource binary_fragmentation_source();

/**
 * Boundary contribution used by this local model.
 *
 * Nucleation injects particles at the lower-size boundary with fixed diameter
 * and constant physical number-rate.
 */
[[nodiscard]] qmomkit::MomentBoundaryContribution boundary_contribution();

} // namespace qmom_1d_breakage_advection_mc_local_rates
