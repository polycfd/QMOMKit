#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

#include <random>

namespace cqmom_4d_particles_breakage_recoil
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity particle_velocity();
[[nodiscard]] qmomkit::CqmomMomentSource breakup_source();
void apply_breakup(qmomkit::monte_carlo::Ensemble& ensemble, double interval, std::mt19937_64& generator);

} // namespace cqmom_4d_particles_breakage_recoil
