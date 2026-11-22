#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/raw_moment_tensor.hpp"

#include <random>
#include <vector>

namespace cqmom_2d_hyperbolic_drag_dispersion
{

constexpr unsigned long long dispersion_seed = 67890;

/**
 * Turbulent dispersion over a step of `step_size` in s: each particle's
 * velocity receives the kick sqrt(2 D(d) dt) xi with xi standard normal, which
 * is exact for diffusion with a constant coefficient over the step.
 */
void apply_dispersion(qmomkit::monte_carlo::Ensemble& ensemble, double step_size, std::mt19937_64& generator);

[[nodiscard]] qmomkit::RawMomentTensor compute_particle_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                                const std::vector<std::size_t>& moment_extents);

} // namespace cqmom_2d_hyperbolic_drag_dispersion
