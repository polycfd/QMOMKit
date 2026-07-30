#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/moment_source.hpp"

#include <random>

namespace cqmom_non_spherical_particles
{

constexpr unsigned long long breakage_seed = 67890;

[[nodiscard]] qmomkit::CqmomMomentSource binary_breakage_source();

void apply_binary_breakage(qmomkit::monte_carlo::Ensemble& ensemble, double step_size, std::mt19937_64& generator);

} // namespace cqmom_non_spherical_particles
