#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/conditional_moments.hpp"

#include <vector>

namespace cqmom_2d_bubbles_cavitation
{

[[nodiscard]] qmomkit::RawMomentTensor compute_particle_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                                const std::vector<std::size_t>& moment_extents);

} // namespace cqmom_2d_bubbles_cavitation
