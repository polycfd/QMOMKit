#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/moment_source.hpp"

#include <random>

namespace gqmom_1d_droplet_breakup_with_mc
{

// The droplet volume V, time t and both parameters are dimensionless. Droplets
// larger than the critical volume break up, at a frequency that rises with
// their excess volume.
constexpr double critical_volume = 1.0;
constexpr double breakup_coefficient = 1.0;
constexpr unsigned long long breakup_seed = 67890;

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

/**
 * Breaks each droplet of the Monte Carlo ensemble within a step of the given
 * size with probability 1 - exp(-g(V) step), where g is breakup_frequency. A
 * broken droplet of volume V becomes two droplets of volume V/2.
 */
void apply_binary_breakup(qmomkit::monte_carlo::Ensemble& droplets, double step_size, std::mt19937_64& generator);

} // namespace gqmom_1d_droplet_breakup_with_mc
