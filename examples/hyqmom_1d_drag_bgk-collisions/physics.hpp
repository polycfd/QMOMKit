#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace hyqmom_1d_drag_bgk_collisions
{

/**
 * Deceleration du/dt of particles of velocity u (m/s) by the drag of a
 * still gas, with a constant drag coefficient.
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity drag_velocity();

/**
 * BGK model of collisions: relaxation of the velocity distribution toward the
 * normal distribution of the same number density, mean and variance.
 */
[[nodiscard]] qmomkit::QmomMomentSource bgk_collision_source();

} // namespace hyqmom_1d_drag_bgk_collisions
