#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/conditional_moments.hpp"

#include <random>
#include <vector>

namespace cqmom_2d_advection_breakage
{

constexpr unsigned long long breakage_seed = 67890;

void apply_binary_breakage(qmomkit::monte_carlo::Ensemble& ensemble, double step_size, std::mt19937_64& generator);

[[nodiscard]] qmomkit::RawMomentTensor compute_particle_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                                const std::vector<std::size_t>& moment_extents);

} // namespace cqmom_2d_advection_breakage
