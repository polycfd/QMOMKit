#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace cqmom_2d_vanishing_variance
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity continuous_velocity();

[[nodiscard]] qmomkit::CqmomMomentSource zero_moment_source();

} // namespace cqmom_2d_vanishing_variance
