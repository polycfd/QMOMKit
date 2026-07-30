#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/conditional_moments.hpp"

#include <vector>

namespace cqmom_non_spherical_particles
{

// Computes weighted, unnormalized raw moments for samples ordered as
// (d, chi, u, v, w), using QMOMKit's dense multi-index flattening convention.
[[nodiscard]] qmomkit::RawMomentTensor compute_ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                                const std::vector<std::size_t>& moment_extents);

} // namespace cqmom_non_spherical_particles
