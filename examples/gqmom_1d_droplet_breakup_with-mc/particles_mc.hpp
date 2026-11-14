#pragma once

#include "monte_carlo/ensemble.hpp"

#include <random>

namespace gqmom_1d_droplet_breakup_with_mc
{

constexpr unsigned long long breakup_seed = 67890;

/**
 * Breaks each droplet of the ensemble within a step of the given size with
 * probability 1 - exp(-g(V) step), where g is breakup_frequency. A broken
 * droplet of volume V becomes two droplets of volume V/2.
 */
void apply_binary_breakup(qmomkit::monte_carlo::Ensemble& droplets, double step_size, std::mt19937_64& generator);

} // namespace gqmom_1d_droplet_breakup_with_mc
