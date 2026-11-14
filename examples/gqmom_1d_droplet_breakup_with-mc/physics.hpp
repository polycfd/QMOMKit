#pragma once

#include "qmomkit/moment_source.hpp"

namespace gqmom_1d_droplet_breakup_with_mc
{

// The droplet volume V, time t and both parameters are dimensionless. Droplets
// larger than the critical volume break up, at a frequency that rises with
// their excess volume.
constexpr double critical_volume = 1.0;
constexpr double breakup_coefficient = 1.0;

/**
 * Breakup frequency of a droplet of volume V, breakup_coefficient times
 * (V - critical_volume) above the critical volume and zero below it.
 */
[[nodiscard]] double breakup_frequency(double volume);

/**
 * Moment source of binary breakup V -> V/2 + V/2, summed over the nodes of the
 * stage quadrature.
 */
[[nodiscard]] qmomkit::QmomMomentSource binary_breakup_source();

} // namespace gqmom_1d_droplet_breakup_with_mc
