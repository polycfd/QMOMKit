#pragma once

#include "qmomkit/qmomkit.hpp"

namespace qmom_1d_minimal
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();
[[nodiscard]] qmomkit::MomentSource fragmentation_source();

} // namespace qmom_1d_minimal
