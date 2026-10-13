#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/distributions.hpp"

#include <cstddef>
#include <vector>

namespace paper_5d_particles_breakup_heat_transfer
{

[[nodiscard]] std::vector<std::size_t> maximum_orders(const std::vector<std::size_t>& extents);

/**
 * The raw moments of the ensemble within the given extents, with each sample
 * carrying `particle_weight`. Throws when the estimate fails.
 */
[[nodiscard]] qmomkit::RawMomentTensor ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                        const std::vector<std::size_t>& extents);

} // namespace paper_5d_particles_breakup_heat_transfer
