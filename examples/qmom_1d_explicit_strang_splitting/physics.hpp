#pragma once

#include "qmomkit/qmomkit.hpp"

namespace qmom_1d_explicit_strang_splitting
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity advection_velocity();
[[nodiscard]] qmomkit::QmomMomentSource variance_injection_source();

} // namespace qmom_1d_explicit_strang_splitting
