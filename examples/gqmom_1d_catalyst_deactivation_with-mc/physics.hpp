#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

#include <random>

namespace gqmom_1d_catalyst_deactivation_with_mc
{

// The coked fraction theta of a particle's surface lies in [0, 1]. Time t and
// all parameters are dimensionless. Coke builds up on every particle, and
// particles whose coked fraction exceeds the critical value are withdrawn for
// regeneration, at a frequency that rises with their excess coked fraction.
constexpr double coking_rate = 0.1;
constexpr double critical_fraction = 0.8;
constexpr double removal_coefficient = 20.0;
constexpr unsigned long long removal_seed = 67890;

/**
 * Coking velocity d(theta)/dt = coking_rate (1 - theta), which slows as the
 * free surface runs out and keeps every particle in [0, 1].
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity coking_velocity();

/**
 * Removal frequency of a particle of coked fraction theta, removal_coefficient
 * times (theta - critical_fraction) above the critical fraction and zero below
 * it.
 */
[[nodiscard]] double removal_frequency(double coked_fraction);

/**
 * Moment source of the removal, summed over the nodes of the stage quadrature.
 */
[[nodiscard]] qmomkit::QmomMomentSource removal_source();

/**
 * Removes each particle of the Monte Carlo ensemble within a step of the given
 * size with probability 1 - exp(-g(theta) step), where g is removal_frequency.
 */
void apply_removal(qmomkit::monte_carlo::Ensemble& particles, double step_size, std::mt19937_64& generator);

} // namespace gqmom_1d_catalyst_deactivation_with_mc
